/// @file vertex_order_fit_dp.C
/// @brief Order clusters so that CLUSTER 0 IS THE VERTEX END, independent of track direction,
///        then fit the LONGEST prefix that still passes a chi2/ndf threshold.
///
/// STEP 1 -- which end is the vertex.
/// The vertex sits on the beam axis, so |xy| is the natural signal, but it is NOT sufficient on
/// its own: a spiralling track curls back and BOTH ends hug the axis (run_0020 entry 16939 ends
/// at |xy| = 24.6 and 27.2 mm, and the physics says the FARTHER one is the vertex). So three
/// signals vote, all of them physical for a particle that loses energy:
///     |xy|      SMALLER at the vertex   (the vertex is on the beam axis)
///     curvature LARGER  at the vertex   (it is fastest there; the spiral tightens as it slows)
///     charge    LOWER   at the vertex   (dE/dx rises toward the stopping end -- Bragg)
/// Measured on 118 a1975 (d,p) backward tracks, curvature and charge agree with each other on
/// 73%. Where all three agree the answer is unambiguous; the vote records its own margin so a
/// 2-1 call can be told from a 3-0 one downstream.
///
/// The order ALONG the track is the clusteriser's own -- only its DIRECTION is set here. Do not
/// let AtGenfitter re-sort by z afterwards (SetUseClusterOrder(kTRUE)); its z sort is what
/// scrambles shallow-pitch tracks in the first place.
///
/// STEP 2 -- longest passing prefix.
/// Grow from the vertex and keep the LONGEST prefix with chi2/ndf < c2Max, not the BEST chi2:
/// chi2 falls monotonically with fewer points (0.00 at 5% of run_0020/9789 against 0.07 at 80%),
/// so "best chi2" selects the fewest clusters and returns a 13%-biased KE. Longest-passing uses
/// the most information that is still consistent.
///
///   root -b -q 'vertex_order_fit_dp.C("run_0020_multifit",9789,0)'
///   root -b -q 'vertex_order_fit_dp.C("",0,0,"list.txt","out.csv")'
#include <algorithm>
#include <fstream>
#include <numeric>
#include <vector>

struct VoOut { double ke, th, c2n, ndf, rms; int nUsed; bool ok; double vxy{-1}, vz{0}; };

/// circle fit that also hands back the centre -- the guiding centre of the spiral
static double voCircleRC(const std::vector<AtHitCluster> &P, int i0, int i1, double &cxo, double &cyo)
{
   const int m = i1 - i0;
   if (m < 5) return -1;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
   for (int i=i0;i<i1;++i){ auto p=P[i].GetPosition(); double x=p.X(), y=p.Y();
      Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y; }
   double C=m*Sxx-Sx*Sx, D=m*Sxy-Sx*Sy, E=m*Sxxx+m*Sxyy-(Sxx+Syy)*Sx;
   double G=m*Syy-Sy*Sy, H=m*Sxxy+m*Syyy-(Sxx+Syy)*Sy;
   double den=2*(C*G-D*D);
   if (std::fabs(den)<1e-9) return -1;
   cxo=(E*G-D*H)/den; cyo=(C*H-D*E)/den;
   double R=0;
   for (int i=i0;i<i1;++i){ auto p=P[i].GetPosition(); R+=std::hypot(p.X()-cxo,p.Y()-cyo); }
   return R/m;
}

static double voCircleR(const std::vector<AtHitCluster> &P, int i0, int i1)
{
   const int m = i1 - i0;
   if (m < 5) return -1;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
   for (int i = i0; i < i1; ++i) { auto p = P[i].GetPosition(); double x=p.X(), y=p.Y();
      Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y; }
   double C=m*Sxx-Sx*Sx, D=m*Sxy-Sx*Sy, E=m*Sxxx+m*Sxyy-(Sxx+Syy)*Sx;
   double G=m*Syy-Sy*Sy, H=m*Sxxy+m*Syyy-(Sxx+Syy)*Sy;
   double den=2*(C*G-D*D);
   if (std::fabs(den)<1e-9) return -1;
   double cx=(E*G-D*H)/den, cy=(C*H-D*E)/den, R=0;
   for (int i=i0;i<i1;++i){ auto p=P[i].GetPosition(); R+=std::hypot(p.X()-cx,p.Y()-cy); }
   return R/m;
}

/// Returns true if the stored order must be REVERSED to put the vertex at index 0.
/// `votes` receives the three individual verdicts (+1 = "reverse", -1 = "keep", 0 = no opinion).
static bool voNeedsReverse(const std::vector<AtHitCluster> &hc, int votes[3])
{
   const int n = hc.size();
   const int k = std::max(5, n / 3);
   auto meanXY = [&](int i0, int i1) { double s = 0; int m = 0;
      for (int i = i0; i < i1 && i < n; ++i) { auto p = hc[i].GetPosition(); s += std::hypot(p.X(), p.Y()); ++m; }
      return m ? s / m : 0.0; };
   auto meanQ = [&](int i0, int i1) { double s = 0; int m = 0;
      for (int i = i0; i < i1 && i < n; ++i) { s += hc[i].GetCharge(); ++m; }
      return m ? s / m : 0.0; };
   const double xyA = meanXY(0, k), xyB = meanXY(n - k, n);
   const double qA = meanQ(0, k), qB = meanQ(n - k, n);
   const double rA = voCircleR(hc, 0, k), rB = voCircleR(hc, n - k, n);
   // each votes +1 if the END (B) looks more like the vertex than the START (A)
   votes[0] = (std::fabs(xyA - xyB) < 1e-6) ? 0 : (xyB < xyA ? +1 : -1);   // vertex is nearer the axis
   votes[1] = (rA < 0 || rB < 0) ? 0 : (rB > rA ? +1 : -1);                // vertex is the faster end
   votes[2] = (qA <= 0 || qB <= 0) ? 0 : (qB < qA ? +1 : -1);              // vertex has the lower dE/dx
   // THE STOPPING END IS INSIDE THE SPIRAL: the radius shrinks about a guiding centre that does
   // not move, so the end farther OUT from that centre is the vertex. Stronger than |xy|, which
   // measures distance from the BEAM AXIS and can be smallest at the stopping end on a curled
   // track (run_0016 entry 13708 votes |xy| = -1, i.e. wrong, while radius and charge are right).
   {
      double cx = 0, cy = 0;
      if (voCircleRC(hc, 0, n, cx, cy) > 0) {
         double dA = 0, dB = 0; int mA = 0, mB = 0;
         for (int i = 0; i < k && i < n; ++i) { auto p = hc[i].GetPosition(); dA += std::hypot(p.X()-cx, p.Y()-cy); ++mA; }
         for (int i = std::max(0, n-k); i < n; ++i) { auto p = hc[i].GetPosition(); dB += std::hypot(p.X()-cx, p.Y()-cy); ++mB; }
         if (mA && mB) votes[1] = ((dB/mB) > (dA/mA)) ? +1 : -1;
      }
   }
   // THE RADIUS DECIDES. Measured on 118 tracks: the stopping end has the smaller radius on 81%,
   // and the 19% that seem to break it have median chi2/ndf 17.1 against 0.02 -- there it is the
   // FITTED VERTEX used as the reference that is wrong, not the radius. Keying on the radius does
   // not require the fit to have succeeded, which is the point: it works on the broken tracks.
   if (votes[1] != 0) return votes[1] > 0;
   if (votes[2] != 0) return votes[2] > 0;
   return votes[0] > 0;
}

static std::unique_ptr<EventFit::AtGenfitter> voFitter(TString eloss)
{
   auto f = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 1.00782503207, 1, eloss.Data(), kFALSE, 2, 5);
   f->SetCatimaMaterial(kTRUE, kTRUE); f->SetCatimaELoss(kTRUE, kFALSE);
   f->SetELossHybrid(kTRUE, 6.61e-5);
   f->SetZPadPlane(1000.0); f->SetMeasSigma(4.0); f->SetSeedFromSpyral(kFALSE);
   f->SetMatEffectsFallback(kFALSE); f->SetThetaWindow(10.0, 170.0);
   f->SetBackExtrapToAxis(kTRUE);
   // *** THE DIRECTION IS DECIDED ONCE, AND IT IS THE CLUSTER ORDER ***
   // SetUseClusterOrder stops the z sort; SetBackwardSeedFix(kFALSE) stops the OTHER reversal.
   // Both are required. Leaving the seed fix on applies `if (backwardSeed) reverse(addSeq)` ON TOP
   // of a vertex-first order, i.e. a DOUBLE reversal that puts the vertex back at the stopping
   // end. Measured: run_0020 entry 16939 went from a good production fit (KE 0.989, chi2/ndf 0.01)
   // to KE 0.083 -- the reordering "broke" a track it had actually ordered correctly.
   // With the seed fix off, iVtx = order.front() = cluster 0 = the vertex, and the seed direction
   // comes from (order[1] - order[0]), which points along the travel direction by construction.
   // So forward and backward tracks need no special case at all.
   f->SetUseClusterOrder(kTRUE);
   f->SetBackwardSeedFix(kFALSE);
   f->Init();
   return f;
}

static VoOut voFit(const AtTrack &proto, const std::vector<AtHitCluster> &use, TString eloss)
{
   VoOut o{0,0,0,0,0,(int)use.size(),false};
   if (use.size() < 6) return o;
   auto fitter = voFitter(eloss);
   AtTrack track = proto;
   *track.GetHitClusterArray() = use;
   AtPatternEvent pe; pe.AddTrack(track);
   AtTrackingEvent tev;
   fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
   if (tev.GetFittedTracks().empty()) return o;
   auto *ft = tev.GetFittedTracks().front().get();
   auto &k = ft->GetKinematicsXtr();
   double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
   auto &sp = ft->GetSmoothedPositions();
   o.ke=k.kineticEnergy; o.th=k.theta*TMath::RadToDeg();
   o.c2n = ndf>0 ? c2/ndf : 1e9; o.ndf=ndf; o.ok=true;
   double s2=0;
   for (auto &c : use) { auto q=c.GetPosition(); double bd=1e30;
      for (auto &m : sp) { double dz=q.Z()-(1000.0-m.Z());
         double d=(q.X()-m.X())*(q.X()-m.X())+(q.Y()-m.Y())*(q.Y()-m.Y())+dz*dz;
         if (d<bd) bd=d; }
      s2+=bd; }
   o.rms=std::sqrt(s2/use.size());
   // WHERE THE VERTEX LANDS. genfit back-extrapolates the fitted state to the beam axis, so a
   // correct fit must put it ON that axis. Fitting from the WRONG end extrapolates away from the
   // beam and the vertex comes out far off it. This is the one internal check that is physical
   // rather than statistical -- chi2 only says a fit is self-consistent, and a fit run backwards
   // down the track is perfectly self-consistent while being wrong.
   auto v = ft->GetVertex();
   o.vxy = std::hypot(v.X(), v.Y());
   o.vz  = v.Z();
   return o;
}

/// dirSource 2 = DON'T DECIDE. Run the longest-passing prefix scan in BOTH directions -- from
///                cluster 0 up, and from the last cluster down -- and keep whichever gives the
///                longer passing prefix (lower chi2 breaks a tie). One of the two IS the vertex
///                end, so this sidesteps the direction determination altogether.
///                Justified here because the arc-length theta that the framework's direction
///                comes from is itself unreliable: AtPRA unwraps the azimuth in HIT-ARRAY order
///                and AtTrackFinderHDBSCAN hands it an unordered array, which moves theta by
///                >10 deg on 40% of tracks and across the 90 deg boundary on 32%.
///                A 2-way choice is also far less bias-prone than picking the best of many
///                subsets: the alternatives are a physical dichotomy, not a free parameter.
/// dirSource 0 = the FRAMEWORK's direction (PRA GeoTheta > 90 deg means the vertex is the
///                high-z_lab end, the same rule backwardSeedFix uses) -- consistent by
///                construction with the rest of the chain.
///           1 = the physics vote (|xy| + curvature + charge), as a cross-check.
/// Either way the two verdicts are BOTH recorded, so a disagreement is visible rather than
/// silently resolved: those are the tracks where a convention needs checking, not the fit.
void vertex_order_fit_dp(TString runTag="run_0020_multifit", Long64_t entry=9789, int tid=0,
                         TString listFile="", TString csvOut="", int dirSource=0,
                         double c2Max=5.0, int stepPct=5, int minPct=25,
                         /// Leave GeoTheta < 90 tracks exactly as production. The forward
                         /// analysis is good; nothing here should touch it.
                         bool forwardUntouched=true,
                         TString gfDir="/mnt/f/a1975/gf_dp_cateloss/",
                         TString geoName="ATTPC_D300torr_v2_geomanager.root",
                         TString eLossTable="proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf=TFile::Open(dir+"/geometry/"+geoName); gf->Get("FAIRGeom"); }
   TString eloss = dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable;

   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   if (listFile.Length()) { std::ifstream in(listFile.Data()); std::string r; Long64_t e; int t;
      while (in>>r>>e>>t) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,t}); }
   else jobs.push_back({runTag,entry,tid});

   std::ofstream csv;
   if (csvOut.Length()) { csv.open(csvOut.Data());
      csv << "run,entry,tid,ncl,reversed,dir_agree,geotheta,vote_xy,vote_R,vote_q,"
             "ke_prod,c2_prod,ke_new,th_new,c2_new,n_new,pct,dirboth,pctfwd,pctbwd,"
             "vxyf,vxyb,kef,keb,c2f,c2b\n"; }
   else printf("\n  %-24s %5s %5s %-5s %-9s | %9s %9s | %9s %8s %9s %6s %5s\n","track","ncl","rev","dir","votes",
               "KE_prod","chi2_prod","KE_new","th_new","chi2_new","nclus","%");

   TString open; TFile *f=nullptr; TTree *tr=nullptr; TClonesArray *te=nullptr; long done=0;
   for (auto &j : jobs) {
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if (rt!=open) { if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} tr=(TTree*)f->Get("cbmsim"); te=nullptr;
         tr->SetBranchAddress("AtTrackingEvent",&te); open=rt; }
      tr->GetEntry(en);
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      auto trks=ev->GetTrackArray(); AtTrack *src=nullptr;
      for (auto &t2:trks) if (t2.GetTrackID()==id) src=&t2;
      if (!src) continue;
      std::vector<AtHitCluster> all = *src->GetHitClusterArray();
      if (all.size() < 12) continue;

      // --- STEP 1: cluster 0 = the vertex end -------------------------------
      int votes[3]={0,0,0};
      const bool revPhys = voNeedsReverse(all, votes);
      // The FRAMEWORK's view: PRA tags backward when GeoTheta > 90, and for those the vertex is
      // the LAST cluster of the ascending-z order, which is what its own reverse() encodes.
      const double geoTh = src->GetGeoTheta() * TMath::RadToDeg();
      const bool taggedBack = std::isfinite(geoTh) && geoTh >= 90.0;
      const bool revFw = taggedBack;
      // dirSource 2 does not reverse up front: the scan tries both ends itself
      // dirSource 2 AND 3 must NOT pre-reverse: both scan the array from each end themselves,
      // so any reversal here would make "forward"/"backward" relative to an already-flipped array.
      // That bug made every fit-free predictor score BELOW chance (31-43%) when compared against
      // the resulting labels -- the labels were inverted for the tracks the physics vote flipped,
      // not the predictors.
      const bool rev = (dirSource >= 2) ? false : ((dirSource == 0) ? revFw : revPhys);
      int dirBoth = 0, pctFwd = 0, pctBwd = 0;
      VoOut best{}; best.ok = false; int bestPct = 0;
      double vxyF=-1, vxyB=-1, keF=-1, keB=-1, c2F=-1, c2B=-1;
      // production result, computed first so the forward guard below can return it directly
      VoOut prod = [&]{ auto fitter=std::make_unique<EventFit::AtGenfitter>(-2.85,2212,1.00782503207,1,
                          eloss.Data(),kFALSE,2,5);
         fitter->SetCatimaMaterial(kTRUE,kTRUE); fitter->SetCatimaELoss(kTRUE,kFALSE);
         fitter->SetELossHybrid(kTRUE,6.61e-5); fitter->SetZPadPlane(1000.0); fitter->SetMeasSigma(4.0);
         fitter->SetMatEffectsFallback(kFALSE); fitter->SetThetaWindow(10.0,170.0);
         fitter->SetBackwardSeedFix(kTRUE); fitter->SetBackExtrapToAxis(kTRUE);
         fitter->SetUseClusterOrder(kFALSE); fitter->Init();
         AtTrack t3=*src; AtPatternEvent pe; pe.AddTrack(t3); AtTrackingEvent tev;
         fitter->FitEvent(&tev,&pe,nullptr,nullptr,nullptr);
         VoOut o{}; o.ok=false;
         if (tev.GetFittedTracks().empty()) return o;
         auto *ft=tev.GetFittedTracks().front().get(); auto &k=ft->GetKinematicsXtr();
         double ndf=ft->GetTrackMetadata()->GetNdf(), c2=ft->GetTrackMetadata()->GetChi2();
         o.ke=k.kineticEnergy; o.th=k.theta*TMath::RadToDeg(); o.c2n=ndf>0?c2/ndf:1e9; o.ok=true;
         o.nUsed=all.size();
         auto v=ft->GetVertex(); o.vxy=std::hypot(v.X(),v.Y()); o.vz=v.Z();
         return o; }();
      const bool dirAgree = (revFw == revPhys);
      if (rev) std::reverse(all.begin(), all.end());

      // *** DO NOT TOUCH FORWARD TRACKS ***
      // The forward analysis is good and must stay bit-for-bit as production. AtGenfitter already
      // guarantees this internally -- backwardSeed is only ever true for GeoTheta > 90, and with it
      // false addSeq == order, so nothing is reversed. Any reordering scheme here has to honour the
      // same boundary, or it "fixes" tracks that were never broken. This is a HARD guard, not a
      // statistical one: a forward track returns the production result untouched.
      if (forwardUntouched && std::isfinite(geoTh) && geoTh < 90.0) {
         best = prod; bestPct = 100; dirBoth = 0; pctFwd = 100; pctBwd = 0;
         if (csvOut.Length()) {
            csv << rt<<","<<en<<","<<id<<","<<all.size()<<",0,"<<(dirAgree?1:0)<<","<<geoTh<<","
                << votes[0]<<","<<votes[1]<<","<<votes[2]<<","
                << prod.ke<<","<<prod.c2n<<","<<prod.ke<<","<<prod.th<<","<<prod.c2n<<","
                << all.size()<<",100,0,100,0,-1,-1,-1,-1,-1,-1\n";
            if (++done%20==0) printf("  %ld / %zu\n", done, jobs.size());
         } else
            printf("  %-24s %5zu %5s %-5s %+d %+d %+d | %9.3f %9.2f | %9.3f %8.2f %9.2f %6zu %5d  FORWARD: untouched\n",
                   (rt+Form(" e%lld",en)).Data(), all.size(), "no", "fwd", votes[0],votes[1],votes[2],
                   prod.ke, prod.c2n, prod.ke, prod.th, prod.c2n, all.size(), 100);
         continue;
      }

      // --- STEP 2: longest prefix that still passes -------------------------
      // scan(rev=false) grows from cluster 0; scan(rev=true) grows from the last cluster.
      auto scan = [&](bool fromEnd, VoOut &out) -> int {
         out = VoOut{}; out.ok = false;
         for (int pct=100; pct>=minPct; pct-=stepPct) {
            int m = (int)std::lround(all.size()*pct/100.0);
            if (m < 6) break;
            std::vector<AtHitCluster> use;
            if (!fromEnd) use.assign(all.begin(), all.begin()+m);
            else { use.assign(all.end()-m, all.end()); std::reverse(use.begin(), use.end()); }
            VoOut o = voFit(*src, use, eloss);
            if (o.ok && o.c2n < c2Max) { out=o; return pct; }   // first from the top = longest
         }
         return 0;
      };
      if (dirSource == 2 || dirSource == 3) {
         VoOut fwd{}, bwd{};
         const int pf = scan(false, fwd), pb = scan(true, bwd);
         bool takeB;
         if (dirSource == 3) {
            // CHOOSE BY WHERE THE VERTEX LANDS: the beam axis is at x = y = 0, so the direction
            // whose fitted vertex sits nearer it is the one that extrapolated the right way.
            // Falls back to the longer prefix when only one direction produced a fit.
            if (fwd.ok && bwd.ok)      takeB = (bwd.vxy < fwd.vxy);
            else if (bwd.ok && !fwd.ok) takeB = true;
            else                        takeB = false;
         } else {
            takeB = (pb > pf) || (pb == pf && pb > 0 && bwd.ok && (!fwd.ok || bwd.c2n < fwd.c2n));
         }
         best = takeB ? bwd : fwd; bestPct = takeB ? pb : pf;
         dirBoth = takeB ? -1 : +1;
         pctFwd = pf; pctBwd = pb;
         vxyF = fwd.ok ? fwd.vxy : -1; vxyB = bwd.ok ? bwd.vxy : -1;
         keF  = fwd.ok ? fwd.ke  : -1; keB  = bwd.ok ? bwd.ke  : -1;
         c2F  = fwd.ok ? fwd.c2n : -1; c2B  = bwd.ok ? bwd.c2n : -1;
      } else {
         bestPct = scan(false, best);
      }
      VoOut prodLate = [&]{ auto fitter=std::make_unique<EventFit::AtGenfitter>(-2.85,2212,1.00782503207,1,
                          eloss.Data(),kFALSE,2,5);
         fitter->SetCatimaMaterial(kTRUE,kTRUE); fitter->SetCatimaELoss(kTRUE,kFALSE);
         fitter->SetELossHybrid(kTRUE,6.61e-5); fitter->SetZPadPlane(1000.0); fitter->SetMeasSigma(4.0);
         fitter->SetMatEffectsFallback(kFALSE); fitter->SetThetaWindow(10.0,170.0);
         fitter->SetBackwardSeedFix(kTRUE); fitter->SetBackExtrapToAxis(kTRUE);
         fitter->SetUseClusterOrder(kFALSE); fitter->Init();
         AtTrack t3=*src; AtPatternEvent pe; pe.AddTrack(t3); AtTrackingEvent tev;
         fitter->FitEvent(&tev,&pe,nullptr,nullptr,nullptr);
         VoOut o{0,0,0,0,0,0,false};
         if (tev.GetFittedTracks().empty()) return o;
         auto *ft=tev.GetFittedTracks().front().get(); auto &k=ft->GetKinematicsXtr();
         double ndf=ft->GetTrackMetadata()->GetNdf(), c2=ft->GetTrackMetadata()->GetChi2();
         o.ke=k.kineticEnergy; o.th=k.theta*TMath::RadToDeg(); o.c2n=ndf>0?c2/ndf:1e9; o.ok=true;
         return o; }();

      if (csvOut.Length()) {
         csv << rt<<","<<en<<","<<id<<","<<all.size()<<","<<(rev?1:0)<<","
             << (dirAgree?1:0)<<","<<geoTh<<","
             << votes[0]<<","<<votes[1]<<","<<votes[2]<<","
             << prod.ke<<","<<prod.c2n<<","<<best.ke<<","<<best.th<<","<<best.c2n<<","
             << best.nUsed<<","<<bestPct<<","<<dirBoth<<","<<pctFwd<<","<<pctBwd<<","
             << vxyF<<","<<vxyB<<","<<keF<<","<<keB<<","<<c2F<<","<<c2B<<"\n";
         if (++done%20==0) printf("  %ld / %zu\n", done, jobs.size());
      } else {
         printf("  %-24s %5zu %5s %-5s %+d %+d %+d | %9.3f %9.2f | %9.3f %8.2f %9.2f %6d %5d\n",
                (rt+Form(" e%lld",en)).Data(), all.size(), rev?"YES":"no",
                dirAgree?"ok":"DIFFER", votes[0],votes[1],votes[2], prod.ke, prod.c2n,
                best.ke, best.th, best.c2n, best.nUsed, bestPct);
      }
   }
   if (f) f->Close();
   if (csvOut.Length()) { csv.close(); printf("  wrote %s (%ld tracks)\n", csvOut.Data(), done); }
}
