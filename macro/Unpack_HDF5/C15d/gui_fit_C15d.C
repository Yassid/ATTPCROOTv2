/// @file gui_fit_C15d.C
/// @brief LIVE single-track genfit refitting with a ROOT GUI: change a fit parameter, press
///        Refit, see the trajectory and the numbers change. About a second per fit.
///
/// This is what the browser viewer cannot be. The HTML page can only ever show fits precomputed
/// on disk; here AtGenfitter is called for real, so ANY parameter can be varied -- ordering,
/// measSigma, iteration counts, seed source, backward-seed fix, material effects, or a subset of
/// the clusters -- against the same track.
///
/// ⚠ PORTED FROM a1975 16C(d,p) 2026-09-07, AND NOT YET ANCHORED. The a1975 original checked its
/// fitter against a KNOWN track -- run_0020 entry 9789 must return KE 0.385, theta 96.92,
/// chi2/ndf 143.66 -- because omitting SetELossHybrid ALONE moved that to KE 4.517. One missing
/// setter makes this a different fitter and then nothing it shows means anything.
/// THERE IS NO SUCH REFERENCE FOR C15d YET. What is reproduced here is fit_batch.sh's call into
/// fitGenfit_C15d.C, matched setter by setter; what is NOT yet done is fitting one track both ways
/// and checking the numbers agree. Do that before trusting a number off this window.
///
/// The C15d differences from the a1975 original, all of them:
///   gas density  6.61e-5 -> 6.5643e-5   (media.geo TargetD2_300, what the production uses)
///   fit file     <run>_genfitter_p.root -> <run>_genfit_<species>.root
///   eloss table  a copy lives in this folder; the a1975 tree is never referenced
/// Geometry (ATTPC_D300torr_v2_geomanager.root), B (-2.85 T) and ZPadPlane (1000 mm) are already
/// identical between the two experiments -- checked against ATTPC.C15d_a2091_D2.par.
///
///   root -l 'gui_fit_C15d.C("run_0013",0,0)'
///   setsid bash -c 'tail -f /dev/null | root -l "gui_fit_C15d.C(\"run_0013\",0,0)"' &
///
/// EXIT: Quit, or the window's [X], calls gSystem->Exit(0) -- a HARD exit. It is NOT
/// gApplication->Terminate(), and this launcher never calls gApplication->Run(): that pairing is
/// what made the gate drawer respawn a new window every time one was closed, ~2 per second,
/// unstoppable from inside the session. See reference-gate-drawers.

#include <algorithm>
#include <numeric>
#include <fstream>
#include <vector>

static std::vector<int> guiArcWalk(const std::vector<AtHitCluster> &hc)
{
   const int n = hc.size();
   std::vector<ROOT::Math::XYZPoint> p(n);
   for (int i = 0; i < n; ++i) p[i] = hc[i].GetPosition();
   double cx=0, cy=0, cz=0;
   for (auto &q : p) { cx+=q.X(); cy+=q.Y(); cz+=q.Z(); }
   cx/=n; cy/=n; cz/=n;
   auto d2=[&](int a,int b){ return (p[a].X()-p[b].X())*(p[a].X()-p[b].X())
        +(p[a].Y()-p[b].Y())*(p[a].Y()-p[b].Y())+(p[a].Z()-p[b].Z())*(p[a].Z()-p[b].Z()); };
   int far=0; double best=-1;
   for(int i=0;i<n;++i){ double d=(p[i].X()-cx)*(p[i].X()-cx)+(p[i].Y()-cy)*(p[i].Y()-cy)+(p[i].Z()-cz)*(p[i].Z()-cz);
      if(d>best){best=d;far=i;} }
   int start=far; best=-1;
   for(int i=0;i<n;++i) if(d2(far,i)>best){best=d2(far,i);start=i;}
   std::vector<char> used(n,0); std::vector<int> out; int cur=start;
   used[cur]=1; out.push_back(cur);
   for(int k=1;k<n;++k){ int nx=-1; double bd=1e30;
      for(int i=0;i<n;++i) if(!used[i]&&d2(cur,i)<bd){bd=d2(cur,i);nx=i;}
      if(nx<0) break; used[nx]=1; out.push_back(nx); cur=nx; }
   return out;
}

static double guiCircleRC(const std::vector<AtHitCluster> &P, int i0, int i1, double &cxOut, double &cyOut)
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
   cxOut=(E*G-D*H)/den; cyOut=(C*H-D*E)/den;
   double R=0;
   for (int i=i0;i<i1;++i){ auto p=P[i].GetPosition(); R+=std::hypot(p.X()-cxOut,p.Y()-cyOut); }
   return R/m;
}

static double guiCircleR(const std::vector<AtHitCluster> &P, int i0, int i1)
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
   double cx=(E*G-D*H)/den, cy=(C*H-D*E)/den, R=0;
   for (int i=i0;i<i1;++i){ auto p=P[i].GetPosition(); R+=std::hypot(p.X()-cx,p.Y()-cy); }
   return R/m;
}

/// Should the stored order be reversed so CLUSTER 0 IS THE VERTEX? Three physical signals vote:
/// the vertex is nearer the beam axis, is the FASTER end (larger curvature radius), and has the
/// LOWER dE/dx (the Bragg rise sits at the stopping end). |xy| only breaks ties -- a spiralling
/// track curls back and both ends can hug the axis (run_0020/16939 ends at 24.6 and 27.2 mm).
static bool guiNeedsReverse(const std::vector<AtHitCluster> &hc, int votes[3])
{
   const int n = hc.size(), k = std::max(5, n/3);
   auto mXY=[&](int a,int b){ double s=0; int m=0;
      for(int i=a;i<b&&i<n;++i){auto p=hc[i].GetPosition(); s+=std::hypot(p.X(),p.Y()); ++m;} return m?s/m:0.0; };
   auto mQ =[&](int a,int b){ double s=0; int m=0;
      for(int i=a;i<b&&i<n;++i){ s+=hc[i].GetCharge(); ++m;} return m?s/m:0.0; };
   double xyA=mXY(0,k), xyB=mXY(n-k,n), qA=mQ(0,k), qB=mQ(n-k,n);
   double rA=guiCircleR(hc,0,k), rB=guiCircleR(hc,n-k,n);
   votes[0] = (std::fabs(xyA-xyB)<1e-6)?0:(xyB<xyA?+1:-1);
   votes[1] = (rA<0||rB<0)?0:(rB>rA?+1:-1);
   votes[2] = (qA<=0||qB<=0)?0:(qB<qA?+1:-1);
   // THE STOPPING END IS INSIDE THE SPIRAL. A particle losing energy curls tighter about a
   // guiding centre that does not move (uniform B), so the far end of the track sits OUTSIDE and
   // the stopping end spirals INWARD. Measuring each end's distance from the guiding centre of
   // the WHOLE track says which is which, and it is stronger than |xy|: |xy| is the distance from
   // the BEAM AXIS, which for a curled-back track can be small at both ends (run_0020/16939 ends
   // at 24.6 and 27.2 mm) and can even be smallest at the stopping end.
   int voteIn = 0;
   {
      double cx=0, cy=0, Rall=guiCircleRC(hc, 0, n, cx, cy);
      if (Rall > 0) {
         double dA=0, dB=0; int mA=0, mB=0;
         for (int i=0;i<k&&i<n;++i){ auto p=hc[i].GetPosition(); dA+=std::hypot(p.X()-cx,p.Y()-cy); ++mA; }
         for (int i=std::max(0,n-k);i<n;++i){ auto p=hc[i].GetPosition(); dB+=std::hypot(p.X()-cx,p.Y()-cy); ++mB; }
         if (mA && mB) { dA/=mA; dB/=mB;
            // the END is the vertex when the END is the one farther OUT
            voteIn = (dB > dA) ? +1 : -1; }
      }
   }
   votes[1] = voteIn ? voteIn : votes[1];   // supersedes the per-third radius comparison
   // THE LOCAL CURVATURE RADIUS DECIDES. Measured on 118 a1975 (d,p) backward tracks: the
   // stopping end has the smaller radius on 81% of them, median 24.2 mm against 40.1 mm at the
   // vertex end -- and the 19% that appear to break the rule have median chi2/ndf 17.1 against
   // 0.02 for the rest, i.e. it is the FITTED VERTEX used as the reference that is wrong there,
   // not the radius. That is the whole argument for keying on the radius: it does not require
   // the fit to have succeeded, so it still works on exactly the tracks that need fixing.
   // Charge and |xy| are consulted only when the radius cannot be measured at both ends.
   if (votes[1] != 0) return votes[1] > 0;
   if (votes[2] != 0) return votes[2] > 0;   // Bragg: lower dE/dx at the vertex
   return votes[0] > 0;                      // |xy|: weakest, a curled track hugs the axis twice
}


#include "pra_theta.C"   // PraTheta / praThetaHits -- AtPRA's arc-vs-z, replayed

class DpFitGui {
public:
   DpFitGui(TString runTag, Long64_t entry, int tid, TString gfDir, TString listFile, TString species)
      : fRunTag(runTag), fEntry(entry), fTid(tid), fGfDir(gfDir), fSpecies(species)
   {
      TString dir = gSystem->Getenv("VMCWORKDIR");
      { TFile *gf = TFile::Open(dir + "/geometry/ATTPC_D300torr_v2_geomanager.root"); gf->Get("FAIRGeom"); }
      fEloss = dir + "/macro/Unpack_HDF5/C15d/proton_D2_300torr.txt";
      if (listFile.Length() && !gSystem->AccessPathName(listFile)) {
         std::ifstream in(listFile.Data()); std::string r; Long64_t e; int t;
         while (in >> r >> e >> t) fList.push_back({TString(r.c_str()), e, t});
         printf("  track list: %zu entries\n", fList.size());
      }
      MakeGui();
      Load();
      Refit();
   }

   // ---- data ----
   void Load()
   {
      fClusters.clear();
      TString fn = fGfDir + fRunTag + "_genfit_" + fSpecies + ".root";
      if (fFile) { fFile->Close(); fFile = nullptr; }
      fFile = TFile::Open(fn);
      if (!fFile || fFile->IsZombie()) { fLog->AddLine(Form("cannot open %s", fn.Data())); return; }
      auto *t = (TTree *)fFile->Get("cbmsim");
      TClonesArray *te = nullptr; t->SetBranchAddress("AtTrackingEvent", &te);
      t->GetEntry(fEntry);
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) { fLog->AddLine("no AtTrackingEvent"); return; }
      auto trks = ev->GetTrackArray();
      for (auto &tr : trks)
         if (tr.GetTrackID() == fTid) { fTrack = tr; fClusters = *tr.GetHitClusterArray(); }
      // --- the STORED fit: what the production (and FIND) actually produced for this track ---
      fStoredPts.clear(); fStoredKE = fStoredTh = fStoredC2 = fStoredVxy = -1; fStoredNpts = 0;
      for (auto &sft : ev->GetFittedTracks()) {
         if (!sft || sft->GetTrackID() != fTid) continue;
         auto &sk = sft->GetKinematicsXtr();
         double sndf = sft->GetTrackMetadata()->GetNdf(), sc2 = sft->GetTrackMetadata()->GetChi2();
         auto sv = sft->GetVertex();
         fStoredKE = sk.kineticEnergy; fStoredTh = sk.theta * TMath::RadToDeg();
         fStoredC2 = (sndf > 0) ? sc2 / sndf : -1;
         fStoredVxy = std::hypot(sv.X(), sv.Y());
         for (auto &m : sft->GetSmoothedPositions())
            fStoredPts.push_back({m.X(), m.Y(), 1000.0 - m.Z()});   // -> cluster frame
         fStoredNpts = (int)fStoredPts.size();
         break;
      }
      if (fClusters.empty()) { fLog->AddLine("track not found"); return; }
      if (fStoredKE > 0)
         fLog->AddLine(Form("  STORED fit (%s): KE %.3f  theta %.2f  chi2/ndf %.3f  vtx|xy| %.1f  pts %d",
                            gSystem->BaseName(fGfDir.Data()), fStoredKE, fStoredTh, fStoredC2,
                            fStoredVxy, fStoredNpts));
      else
         fLog->AddLine("  STORED fit: none in this production for this track");

      // Replay AtPRA's theta determination on the UNTOUCHED hit array, before any reordering
      // below -- pane 4 must show what the production PRA saw, not what this GUI rearranged.
      fPra = praThetaHits(fTrack);
      if (fPra.ok)
         fLog->AddLine(Form("  PRA theta plane: %d/%d hits on the RANSAC line, slope sign %+d"
                            " -> theta %.1f deg (stored GeoTheta %.1f)",
                            (int)std::count(fPra.inlier.begin(), fPra.inlier.end(), 1),
                            (int)fPra.arc.size(), fPra.sign, fPra.angleDeg,
                            fTrack.GetGeoTheta() * TMath::RadToDeg()));
      else
         fLog->AddLine("  PRA theta plane: no RANSAC line (pane 4 empty)");

      // *** CLUSTER 0 IS THE VERTEX, AND IT IS FIXED HERE -- ONCE, FOR EVERYTHING ***
      // This used to be done inside Refit, AFTER the index window had already been applied, so
      // "use clusters 0..30" still meant the first 30 of the STORED order -- the STOPPING end on
      // ~90% of tracks (run_0016 entry 13708: first third R 30.8 mm and charge 10856, i.e. the
      // stopping end, while the vertex is the last cluster). Reordering the array itself means the
      // index window, the box selection, the markers and the fit all speak the same order, and
      // index 0 means the vertex everywhere.
      // *** THE ARRAY IS PUT INTO AtGenfitter's OWN MEASUREMENT ORDER (addSeq) ***
      // It used to be the STORED cluster order with an optional reversal. That was wrong for an
      // index window, because the stored order is AtPRA's greedy nearest-neighbour walk, which
      // must consume every cluster and therefore APPENDS the ones it skipped at the end. Measured
      // on 797 long tracks: z total-variation / z span is 1.06 below 100 clusters but 2.20 at
      // 100-300, 3.64 at 300-600 -- i.e. above ~300 clusters the walk crosses the z range three or
      // four times, and 58 % of long tracks U-turn. Reversing such an array puts the LEFTOVERS
      // first: on run_0016 e31316 the first 100 clusters after reversal were two disconnected
      // clumps at z ~ 530 and z ~ 280, not the vertex region at all. A window of "clusters 0 to N"
      // then means nothing, and cannot be compared with what FIND truncates.
      //
      // AtGenfitter ignores the cluster order entirely (unless SetUseClusterOrder): it sorts by
      // z_lab = 1000 - z_cluster ascending, and with backwardSeedFix a PRA-backward track is fitted
      // from the far end of that sort. Reproducing exactly that here makes index 0 genfit's FIRST
      // MEASUREMENT POINT, so an index window is the same prefix FIND would take.
      if (fVertexFirst && fVertexFirst->IsOn()) {
         const double gt = fTrack.GetGeoTheta() * TMath::RadToDeg();
         const bool backwardSeed = std::isfinite(gt) && gt >= 90.0;
         std::vector<int> ord(fClusters.size());
         std::iota(ord.begin(), ord.end(), 0);
         // z_lab ascending == z_cluster DESCENDING
         std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            return fClusters[a].GetPosition().Z() > fClusters[b].GetPosition().Z(); });
         if (backwardSeed) std::reverse(ord.begin(), ord.end());
         std::vector<AtHitCluster> re; re.reserve(ord.size());
         for (int i : ord) re.push_back(fClusters[i]);
         fClusters = std::move(re);
         // how badly was the stored order broken? same metric as zuturn_dp.C, for context.
         double tv = 0, zlo = 1e30, zhi = -1e30;
         for (size_t i = 0; i < fClusters.size(); ++i) {
            double z = fClusters[i].GetPosition().Z();
            zlo = std::min(zlo, z); zhi = std::max(zhi, z);
         }
         int v[3] = {0, 0, 0};
         const bool revPhys = guiNeedsReverse(fClusters, v);
         fLog->AddLine(Form("  cluster order: PRODUCTION (z sort%s) -- index 0 is genfit's first "
                            "measurement.  GeoTheta %.1f -> backwardSeed %s ; physics vote says %s",
                            backwardSeed ? " + backwardSeedFix reversal" : "", gt,
                            backwardSeed ? "ON" : "off", revPhys ? "reverse" : "keep"));
         fLog->AddLine(Form("  z now runs %.1f -> %.1f mm monotonically (the stored PRA walk order "
                            "is discarded; above ~300 clusters it U-turns on ~85%% of tracks)",
                            fClusters.front().GetPosition().Z(), fClusters.back().GetPosition().Z()));
      } else {
         fLog->AddLine("  cluster order: RAW STORED (AtPRA greedy walk) -- an index window here is "
                       "NOT a contiguous piece of track on long tracks");
      }
      fI0->SetNumber(0);
      fI1->SetNumber(fClusters.size() - 1);
      fKeep.assign(fClusters.size(), 1);   // a fresh track starts with everything selected
   }

   /// The cluster list Refit() would use right now -- the index window ANDed with the box mask.
   /// Kept as its own method so a view-only redraw cannot drift from what the fit actually used.
   std::vector<AtHitCluster> CurrentSelection() const
   {
      std::vector<AtHitCluster> use;
      const int i0 = std::max(0, (int)fI0->GetNumber());
      const int i1 = std::min((int)fClusters.size() - 1, (int)fI1->GetNumber());
      for (int i = i0; i <= i1; ++i)
         if (i < (int)fKeep.size() && fKeep[i]) use.push_back(fClusters[i]);
      return use;
   }

   void Refit()
   {
      if (fClusters.empty()) { fLog->AddLine("no clusters loaded"); return; }
      const int mode = fOrder->GetSelected();
      const int i0 = std::max(0, (int)fI0->GetNumber());
      const int i1 = std::min((int)fClusters.size() - 1, (int)fI1->GetNumber());
      // The index window and the graphical mask are ANDed: the window is the coarse cut, the
      // boxes are the surgical one. Which clusters go into the fit is the whole point of the
      // tool, so the fit and the drawing use ONE list -- there is no way for the picture to
      // show a selection the fit did not actually use.
      std::vector<AtHitCluster> use;
      fUsedIdx.clear();
      for (int i = i0; i <= i1; ++i)
         if (i < (int)fKeep.size() && fKeep[i]) { use.push_back(fClusters[i]); fUsedIdx.push_back(i); }
      if (mode == 2) {
         auto ord = guiArcWalk(use);
         std::vector<AtHitCluster> re; re.reserve(ord.size());
         for (int i : ord) re.push_back(use[i]);
         use = re;
      }
      // mode 3 needs no reordering here: fClusters was already put vertex-first in Load(), so
      // `use` is a vertex-first prefix by construction. Mode 3 only tells the fitter to KEEP that
      // order (no z sort) and to skip its own reversal -- doing both would be a double reversal.
      if ((int)use.size() < 4) { fLog->AddLine("fewer than 4 clusters selected"); return; }

      const bool matEff = fMatFx->IsOn();
      auto fitter = std::make_unique<EventFit::AtGenfitter>(
         -2.85, 2212, 1.00782503207, 1, fEloss.Data(), !matEff,
         (int)fMinIt->GetNumber(), (int)fMaxIt->GetNumber());
      if (fCatima->IsOn()) { fitter->SetCatimaMaterial(kTRUE, kTRUE); fitter->SetCatimaELoss(kTRUE, kFALSE); }
      fitter->SetELossHybrid(kTRUE, 6.5643e-5);    // C15d D2 300 torr; fitGenfit_C15d.C uses the same
      fitter->SetZPadPlane(1000.0);
      fitter->SetMeasSigma(fSigma->GetNumber());
      fitter->SetSeedFromSpyral(fSpySeed->IsOn());
      fitter->SetMatEffectsFallback(kFALSE);
      fitter->SetThetaWindow(10.0, 170.0);
      // mode 3 supplies the direction through the ORDER, so the framework reversal is forced
      // off: the two must never both be active.
      fitter->SetBackwardSeedFix(mode == 3 ? kFALSE : fBackFix->IsOn());
      fitter->SetBackExtrapToAxis(kTRUE);
      fitter->SetUseClusterOrder(mode != 0);
      fitter->Init();

      AtTrack track = fTrack;
      *track.GetHitClusterArray() = use;
      AtPatternEvent pe; pe.AddTrack(track);
      AtTrackingEvent tev;
      fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);

      fFitPts.clear();
      fLastC2 = -1;
      if (tev.GetFittedTracks().empty()) {
         fLog->AddLine(Form("%s e%lld t%d  mode=%d  -> NO FIT RETURNED", fRunTag.Data(), fEntry, fTid, mode));
         Draw(use);
         return;
      }
      auto *ft = tev.GetFittedTracks().front().get();
      auto &k = ft->GetKinematicsXtr();
      auto v = ft->GetVertex();
      double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
      auto &sp = ft->GetSmoothedPositions();
      for (auto &m : sp) fFitPts.push_back({m.X(), m.Y(), 1000.0 - m.Z()}); // -> cluster frame
      double s2 = 0;
      for (auto &c : use) {
         auto q = c.GetPosition(); double bd = 1e30;
         for (auto &m : fFitPts) {
            double d = (q.X()-m[0])*(q.X()-m[0]) + (q.Y()-m[1])*(q.Y()-m[1]) + (q.Z()-m[2])*(q.Z()-m[2]);
            if (d < bd) bd = d;
         }
         s2 += bd;
      }
      double rms = std::sqrt(s2 / use.size());
      fLastC2 = ndf > 0 ? c2 / ndf : -1;
      const char *mn[] = {"z-sort", "cluster", "arcwalk"};
      {
         const int iF = GenfitFirstIndex(use, mode), iL = GenfitLastIndex(use, mode);
         auto p0 = use[iF].GetPosition(); auto pN = use[iL].GetPosition();
         // THE INDEPENDENT CHECK: the fitted vertex should be nearer genfit's point 0 than its
         // last point. If it is not, the ordering IS reversed, whatever the picture suggests.
         double dF = std::hypot(std::hypot(p0.X()-v.X(), p0.Y()-v.Y()), p0.Z()-v.Z());
         double dL = std::hypot(std::hypot(pN.X()-v.X(), pN.Y()-v.Y()), pN.Z()-v.Z());
         fLog->AddLine(Form("  genfit pt0 = cluster %d (%.1f,%.1f,%.1f) |xy| %.1f q %.0f  d_vtx %.0f"
                            "   |  last = cluster %d |xy| %.1f q %.0f  d_vtx %.0f   %s",
                            iF, p0.X(), p0.Y(), p0.Z(), std::hypot(p0.X(),p0.Y()), use[iF].GetCharge(), dF,
                            iL, std::hypot(pN.X(),pN.Y()), use[iL].GetCharge(), dL,
                            dF < dL ? "OK (pt0 nearer the vertex)" : "<<< REVERSED (pt0 is the FAR end)"));
      }
      fLog->AddLine(Form("%s e%lld t%d | %-7s sig %.1f it %d-%d %s%s | KE %7.3f  th %6.2f  chi2/ndf %8.2f  ndf %.0f  pts %zu  rms %.1f",
                         fRunTag.Data(), fEntry, fTid, mn[mode], fSigma->GetNumber(),
                         (int)fMinIt->GetNumber(), (int)fMaxIt->GetNumber(),
                         matEff ? "matFX " : "noMat ", fBackFix->IsOn() ? "bwFix" : "     ",
                         k.kineticEnergy, k.theta * TMath::RadToDeg(), ndf > 0 ? c2 / ndf : -1, ndf,
                         sp.size(), rms));
      fLog->ShowBottom();
      if (fSelLbl) { fSelLbl->SetText(Form("selected: %d / %d", (int)use.size(), (int)fClusters.size()));
                     }
      Draw(use);
   }

   /// The index in `use` that AtGenfitter will actually treat as measurement point 0.
   /// This MIRRORS AtGenfitter::GetFittedTrack: with useClusterOrder off it sorts by ascending
   /// z_lab (= ZPadPlane - z), and with backwardSeed it reverses the whole sequence. Marking
   /// use[0] instead was WRONG for the z-sort modes -- the star sat on a cluster genfit does not
   /// treat as first, which reads as "the ordering is still reversed" when it is only the picture
   /// that is. Anything drawn as "cluster 0" has to come from the same rule the fitter applies.
   int GenfitFirstIndex(const std::vector<AtHitCluster> &use, int mode) const
   {
      const int n = use.size();
      if (n == 0) return 0;
      std::vector<int> order(n);
      for (int i = 0; i < n; ++i) order[i] = i;
      const bool useClusterOrder = (mode != 0);
      if (!useClusterOrder)
         std::sort(order.begin(), order.end(), [&](int a, int b) {
            return (1000.0 - use[a].GetPosition().Z()) < (1000.0 - use[b].GetPosition().Z()); });
      const bool backFixOn = (mode == 3) ? false : fBackFix->IsOn();
      const double gt = fTrack.GetGeoTheta();
      const bool backwardSeed = backFixOn && std::isfinite(gt) && gt > M_PI / 2.0;
      return backwardSeed ? order.back() : order.front();
   }
   int GenfitLastIndex(const std::vector<AtHitCluster> &use, int mode) const
   {
      const int n = use.size();
      if (n == 0) return 0;
      std::vector<int> order(n);
      for (int i = 0; i < n; ++i) order[i] = i;
      if (mode == 0)
         std::sort(order.begin(), order.end(), [&](int a, int b) {
            return (1000.0 - use[a].GetPosition().Z()) < (1000.0 - use[b].GetPosition().Z()); });
      const bool backFixOn = (mode == 3) ? false : fBackFix->IsOn();
      const double gt = fTrack.GetGeoTheta();
      const bool backwardSeed = backFixOn && std::isfinite(gt) && gt > M_PI / 2.0;
      return backwardSeed ? order.front() : order.back();
   }

   void Draw(const std::vector<AtHitCluster> &use)
   {
      const int iFirst = GenfitFirstIndex(use, fOrder->GetSelected());
      const int iLast  = GenfitLastIndex(use, fOrder->GetSelected());
      auto *c = fCanvas->GetCanvas();
      c->Clear(); c->Divide(2, 2);
      const int ax[3] = {2, 2, 0}, ay[3] = {0, 1, 1};
      const char *xl[3] = {"z [mm]", "z [mm]", "x [mm]"};
      const char *yl[3] = {"x [mm]", "y [mm]", "y [mm]"};
      for (int pane = 0; pane < 3; ++pane) {
         c->cd(pane + 1);
         fPads[pane] = gPad;
         gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12);
         auto get = [&](const ROOT::Math::XYZPoint &p, int a) { return a == 0 ? p.X() : (a == 1 ? p.Y() : p.Z()); };
         auto getv = [&](const std::array<double,3> &p, int a) { return p[a]; };
         double lo1 = 1e30, hi1 = -1e30, lo2 = 1e30, hi2 = -1e30;
         // range over EVERY cluster, so the view does not jump each time hits are dropped
         for (auto &h : fClusters) { auto p = h.GetPosition();
            lo1 = std::min(lo1, get(p, ax[pane])); hi1 = std::max(hi1, get(p, ax[pane]));
            lo2 = std::min(lo2, get(p, ay[pane])); hi2 = std::max(hi2, get(p, ay[pane])); }
         for (auto &q : fFitPts) {   // the fit is inside the frame: a trajectory that wanders
            lo1 = std::min(lo1, getv(q, ax[pane])); hi1 = std::max(hi1, getv(q, ax[pane]));   // off the
            lo2 = std::min(lo2, getv(q, ay[pane])); hi2 = std::max(hi2, getv(q, ay[pane])); } // cloud is the point
         double m1 = std::max(8.0, (hi1 - lo1) * 0.10), m2 = std::max(8.0, (hi2 - lo2) * 0.10);
         auto *fr = gPad->DrawFrame(lo1 - m1, lo2 - m2, hi1 + m1, hi2 + m2);
         fr->SetTitle(Form(";%s;%s", xl[pane], yl[pane]));
         double qmax = 1;
         for (auto &h : use) qmax = std::max(qmax, (double)h.GetCharge());
         // dropped clusters stay visible in grey: a selection you cannot see is a selection you
         // cannot check, and silently vanishing hits is how a cut gets forgotten.
         auto *gd = new TGraph();
         int nd = 0;
         for (size_t i = 0; i < fClusters.size(); ++i) {
            bool used = std::find(fUsedIdx.begin(), fUsedIdx.end(), (int)i) != fUsedIdx.end();
            if (used) continue;
            auto p = fClusters[i].GetPosition();
            gd->SetPoint(nd++, get(p, ax[pane]), get(p, ay[pane]));
         }
         if (nd) { gd->SetMarkerStyle(20); gd->SetMarkerSize(0.35); gd->SetMarkerColor(kGray + 1);
                   gd->Draw("P same"); }
         auto *g = new TGraph();
         for (size_t i = 0; i < use.size(); ++i) { auto p = use[i].GetPosition();
            g->SetPoint(i, get(p, ax[pane]), get(p, ay[pane])); }
         g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(kAzure + 1);
         g->Draw("P same");
         // FIRST and LAST cluster of the order actually handed to the fitter. This is the whole
         // point of the vertex-first mode -- without seeing where cluster 0 landed there is no way
         // to tell a correct ordering from a reversed one, and a reversed one still draws a
         // perfectly plausible track.
         if (!use.empty()) {
            auto p0 = use[iFirst].GetPosition();
            auto *m0 = new TMarker(get(p0, ax[pane]), get(p0, ay[pane]), 29); // filled star
            m0->SetMarkerColor(kGreen + 2); m0->SetMarkerSize(2.6); m0->Draw();
            auto *m0b = new TMarker(get(p0, ax[pane]), get(p0, ay[pane]), 30); // open star on top
            m0b->SetMarkerColor(kBlack); m0b->SetMarkerSize(2.6); m0b->Draw();
            if (use.size() > 1) {
               auto pN = use[iLast].GetPosition();
               auto *mN = new TMarker(get(pN, ax[pane]), get(pN, ay[pane]), 21); // filled square
               mN->SetMarkerColor(kOrange + 7); mN->SetMarkerSize(1.5); mN->Draw();
               auto *mNb = new TMarker(get(pN, ax[pane]), get(pN, ay[pane]), 25);
               mNb->SetMarkerColor(kBlack); mNb->SetMarkerSize(1.5); mNb->Draw();
            }
            if (pane == 0) {
               auto *lg = new TLegend(0.55, 0.78, 0.98, 0.98);
               lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
               auto *s0 = new TMarker(0,0,29); s0->SetMarkerColor(kGreen+2); s0->SetMarkerSize(2.0);
               auto *sN = new TMarker(0,0,21); sN->SetMarkerColor(kOrange+7); sN->SetMarkerSize(1.3);
               lg->AddEntry(s0, "genfit point 0  (vertex end)", "p");
               lg->AddEntry(sN, "genfit last point  (far end)", "p");
               auto *lS = new TGraph(2); lS->SetLineColor(kGreen+2); lS->SetLineWidth(3);
               auto *lL = new TGraph(2); lL->SetLineColor(kRed+1); lL->SetLineWidth(2); lL->SetLineStyle(2);
               lg->AddEntry(lS, "STORED production fit (what FIND chose)", "l");
               lg->AddEntry(lL, "live refit in this window (no FIND)", "l");
               lg->Draw();
            }
         }
         // ---- optional layer: the RAW HITS under the clusters ---------------------------
         // The clusters are charge-weighted centroids; the hits are what they were made from.
         // Seeing both tells a clustering artefact from a real feature of the track.
         if (fShowHits && fShowHits->IsOn()) {
            auto *gh = new TGraph();
            int nh = 0;
            for (auto &h : fTrack.GetHitArray()) {
               auto p = h->GetPosition();
               gh->SetPoint(nh++, get(p, ax[pane]), get(p, ay[pane]));
            }
            if (nh) { gh->SetMarkerStyle(1); gh->SetMarkerColor(kGray + 2); gh->Draw("P same"); }
         }
         // ---- optional layer: the STORED CLUSTER ORDER, drawn as a stroke ------------------
         // This is the thing to judge. A good ordering is one clean stroke along the track; a
         // broken one throws long chords across it (median largest jump measured at 124 mm for
         // the nearest-neighbour walk, 77 mm for a z sort, on 800 long tracks).
         if (fShowOrder && fShowOrder->IsOn() && fClusters.size() > 1) {
            auto *go = new TGraph();
            for (size_t i = 0; i < fClusters.size(); ++i) {
               auto p = fClusters[i].GetPosition();
               go->SetPoint(i, get(p, ax[pane]), get(p, ay[pane]));
            }
            go->SetLineColor(kMagenta + 1); go->SetLineWidth(1); go->Draw("L same");
            // mark where the order jumps, so a bad link is visible and not just implied
            double mean = 0; int nst = 0;
            for (size_t i = 1; i < fClusters.size(); ++i) {
               mean += (fClusters[i].GetPosition() - fClusters[i-1].GetPosition()).R(); ++nst; }
            if (nst) mean /= nst;
            for (size_t i = 1; i < fClusters.size(); ++i) {
               double d = (fClusters[i].GetPosition() - fClusters[i-1].GetPosition()).R();
               if (d <= 3.0 * mean) continue;
               auto p0 = fClusters[i-1].GetPosition(), p1 = fClusters[i].GetPosition();
               auto *jl = new TGraph(2);
               jl->SetPoint(0, get(p0, ax[pane]), get(p0, ay[pane]));
               jl->SetPoint(1, get(p1, ax[pane]), get(p1, ay[pane]));
               jl->SetLineColor(kRed); jl->SetLineWidth(2); jl->SetLineStyle(2); jl->Draw("L same");
            }
         }
         // THE STORED PRODUCTION FIT (green): what FIND actually chose and what the physics was
         // built from. Drawn UNDER the live refit so the two are visibly distinct -- if they
         // differ, the window is telling you the production did something the GUI would not.
         if (fStoredPts.size() > 1) {
            auto *gs = new TGraph();
            for (size_t i = 0; i < fStoredPts.size(); ++i)
               gs->SetPoint(i, getv(fStoredPts[i], ax[pane]), getv(fStoredPts[i], ay[pane]));
            gs->SetLineColor(kGreen + 2); gs->SetLineWidth(3); gs->SetLineStyle(1);
            gs->Draw("L same");
         }
         if (fFitPts.size() > 1) {
            auto *gf = new TGraph();
            for (size_t i = 0; i < fFitPts.size(); ++i)
               gf->SetPoint(i, getv(fFitPts[i], ax[pane]), getv(fFitPts[i], ay[pane]));
            gf->SetLineColor(kRed + 1); gf->SetLineWidth(2); gf->SetLineStyle(2); gf->Draw("L same");
         }
      }

      // ---- pane 4: AtPRA's (arc length, z) plane and the line RANSAC that sets GeoTheta -------
      // The whole forward/backward decision is the SIGN of this one line. Drawing it is the only
      // way to see whether that sign was taken from a real line or from a scrambled cloud.
      c->cd(4);
      fPads[3] = gPad;
      gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.12);
      if (!fPra.ok || fPra.arc.empty()) {
         auto *fr = gPad->DrawFrame(0, 0, 1, 1);
         fr->SetTitle("AtPRA theta: no RANSAC line;arc length [mm];z [mm]");
      } else {
         double alo = *std::min_element(fPra.arc.begin(), fPra.arc.end());
         double ahi = *std::max_element(fPra.arc.begin(), fPra.arc.end());
         double zlo = *std::min_element(fPra.z.begin(), fPra.z.end());
         double zhi = *std::max_element(fPra.z.begin(), fPra.z.end());
         double ma = std::max(5.0, (ahi - alo) * 0.08), mz = std::max(5.0, (zhi - zlo) * 0.08);
         auto *fr = gPad->DrawFrame(alo - ma, zlo - mz, ahi + ma, zhi + mz);
         // stored GeoTheta beside the recomputed one: they must agree, and if they do not the
         // track was reconstructed by a different code path than this pane is replaying.
         fr->SetTitle(Form("RANSAC theta %.1f#circ (sign %+d, %d/%d inliers) | stored GeoTheta %.1f#circ"
                           ";arc length [mm];z [mm]",
                           fPra.angleDeg, fPra.sign, (int)std::count(fPra.inlier.begin(), fPra.inlier.end(), 1),
                           (int)fPra.arc.size(), fTrack.GetGeoTheta() * TMath::RadToDeg()));
         // outliers grey, inliers blue -- RANSAC discards up to a third of the points and the
         // slope is taken from what survives, so which points survived is part of the answer.
         auto *go = new TGraph(); auto *gi = new TGraph();
         int no = 0, ni = 0;
         for (size_t i = 0; i < fPra.arc.size(); ++i) {
            if (fPra.inlier[i]) gi->SetPoint(ni++, fPra.arc[i], fPra.z[i]);
            else                go->SetPoint(no++, fPra.arc[i], fPra.z[i]);
         }
         if (no) { go->SetMarkerStyle(20); go->SetMarkerSize(0.35); go->SetMarkerColor(kGray + 1);
                   go->Draw("P same"); }
         if (ni) { gi->SetMarkerStyle(20); gi->SetMarkerSize(0.5); gi->SetMarkerColor(kAzure + 1);
                   gi->Draw("P same"); }
         // the fitted line, drawn across the frame
         double t0 = (alo - ma - fPra.px) / (std::fabs(fPra.dirX) > 1e-9 ? fPra.dirX : 1e-9);
         double t1 = (ahi + ma - fPra.px) / (std::fabs(fPra.dirX) > 1e-9 ? fPra.dirX : 1e-9);
         auto *ln = new TGraph();
         ln->SetPoint(0, fPra.px + t0 * fPra.dirX, fPra.py + t0 * fPra.dirY);
         ln->SetPoint(1, fPra.px + t1 * fPra.dirX, fPra.py + t1 * fPra.dirY);
         ln->SetLineColor(kRed + 1); ln->SetLineWidth(2); ln->Draw("L same");
         // first point of the array: where the unwrapping started
         auto *m0 = new TMarker(fPra.arc.front(), fPra.z.front(), 29);
         m0->SetMarkerColor(kGreen + 2); m0->SetMarkerSize(2.4); m0->Draw();
      }
      c->Update();
   }

   // ---- slots ----
   void OnRefit() { Refit(); }
   /// view-only: redraw with the current selection, WITHOUT refitting
   void OnRedraw() { Draw(CurrentSelection()); }
   /// Re-read the track so the reordering is redone from the untouched cluster array -- reversing
   /// in place twice would just flip back and forth.
   void OnReorder() { Load(); Refit(); }
   void OnReload() { fRunTag = fRunEntry->GetText(); fEntry = (Long64_t)fEntryNum->GetNumber();
                     fTid = (int)fTidNum->GetNumber(); fI0->SetNumber(0); Load(); Refit(); }
   void OnStep(int d)
   {
      if (fList.empty()) return;
      fIdx = (fIdx + d + (int)fList.size()) % fList.size();
      fRunTag = "run_" + std::get<0>(fList[fIdx]) + "_multifit";
      fEntry = std::get<1>(fList[fIdx]); fTid = std::get<2>(fList[fIdx]);
      fRunEntry->SetText(fRunTag); fEntryNum->SetNumber(fEntry); fTidNum->SetNumber(fTid);
      fI0->SetNumber(0);
      Load(); Refit();
   }
   void OnNext() { OnStep(+1); }
   void OnPrev() { OnStep(-1); }
   void OnAllModes()
   {
      int keep = fOrder->GetSelected();
      for (int m = 0; m < 3; ++m) { fOrder->Select(m, kFALSE); Refit(); }
      fOrder->Select(keep, kFALSE);
   }
   void OnSavePng() { fCanvas->GetCanvas()->SaveAs(Form("plots/gui_fit_%s_e%lld_t%d.png",
                                                        fRunTag.Data(), fEntry, fTid)); }
   /// Fit the LONGEST prefix (from cluster 0 -- the vertex, in mode 3) still under the chi2/ndf
   /// threshold. NOT the best chi2: chi2 falls monotonically with fewer points, so "best" picks
   /// the fewest clusters -- on run_0020/9789 that is 15 of 298 and a 13%-biased KE.
   /// The production FIND rule, reproduced exactly (AtGenfitter::SetFindMode(1)):
   ///   full length is kept IF it is already under c2Max -- it is the longest prefix AND passing,
   ///   so nothing shorter can beat it (a smaller chi2 below the limit is fewer points, not a
   ///   better measurement). If it FAILS, walk down and keep the LONGEST rung under c2Max; and
   ///   only if NOTHING gets under it does the outright MINIMUM chi2/ndf decide, so a hopeless
   ///   track still ends on its best available fit rather than on the last rung tried.
   void OnLongest()
   {
      if (fClusters.empty()) return;
      const double c2Max = fC2Max->GetNumber();
      const int n = fClusters.size();
      // down to 5 %. The 25 % floor was inherited from the original coarse ladder and had no
      // evidence behind it: on run_0016 e11488 a better fit sits at 14 % of the track, which a
      // floor of 25 could never find. The 8-cluster guard below still stops it going absurd.
      const int step = 5, floorPct = 5;
      int bestPct = -1;                 // longest passing
      int minPct = -1; double minC2 = 1e30;   // fallback
      for (int pct = 100; pct >= floorPct; pct -= step) {
         int m = (int)std::lround(n * pct / 100.0);
         if (m < 8) break;
         fI0->SetNumber(0); fI1->SetNumber(m - 1);
         Refit();
         if (fLastC2 <= 0) continue;
         if (fLastC2 < minC2) { minC2 = fLastC2; minPct = pct; }
         if (fLastC2 < c2Max && bestPct < 0) { bestPct = pct; break; }  // first = longest
      }
      const int keep = (bestPct > 0) ? bestPct : minPct;
      if (keep < 0) {
         fLog->AddLine(Form("  no prefix from 100 %% down to %d %% returned a fit at all", floorPct));
         fLog->ShowBottom();
         return;
      }
      // re-apply the winner: the loop leaves its LAST attempt selected, which would otherwise be
      // shown and fitted instead of the chosen one.
      const int m = (int)std::lround(n * keep / 100.0);
      fI0->SetNumber(0); fI1->SetNumber(m - 1);
      Refit();
      if (bestPct > 0)
         fLog->AddLine(Form("  LONGEST PASSING: %d %% = %d clusters, chi2/ndf %.3f < %.3f",
                            keep, m, fLastC2, c2Max));
      else
         fLog->AddLine(Form("  nothing reached chi2/ndf < %.3f -- FALLBACK to the MINIMUM: "
                            "%d %% = %d clusters, chi2/ndf %.3f", c2Max, keep, m, fLastC2));
      fLog->ShowBottom();
   }

   /// Rubber-band selection. Drag with the left button on ANY of the three projections; the
   /// clusters inside the box are kept or dropped according to the radio setting. The box is
   /// applied in the two coordinates of THAT pad only, so a cut made on the pad plane does not
   /// silently constrain z, and a cut in z does not constrain the transverse plane.
   void HandleEvent(Int_t event, Int_t px, Int_t py, TObject *)
   {
      if (fClusters.empty()) return;
      int pane = -1;
      for (int i = 0; i < 3; ++i) if (gPad == fPads[i]) pane = i;
      if (pane < 0) return;
      if (event == kButton1Down) {
         fDragX0 = gPad->AbsPixeltoX(px); fDragY0 = gPad->AbsPixeltoY(py);
         fDragging = true; fDragPane = pane; return;
      }
      if (event != kButton1Up || !fDragging || pane != fDragPane) return;
      fDragging = false;
      double x1 = gPad->AbsPixeltoX(px), y1 = gPad->AbsPixeltoY(py);
      double xlo = std::min(fDragX0, x1), xhi = std::max(fDragX0, x1);
      double ylo = std::min(fDragY0, y1), yhi = std::max(fDragY0, y1);
      // a click, not a drag: ignore rather than wipe the selection on a stray click
      if (std::fabs(xhi - xlo) < 1e-6 || std::fabs(yhi - ylo) < 1e-6) return;
      const int ax[3] = {2, 2, 0}, ay[3] = {0, 1, 1};
      auto co = [&](const ROOT::Math::XYZPoint &p, int a) { return a == 0 ? p.X() : (a == 1 ? p.Y() : p.Z()); };
      int n = 0;
      const bool keepMode = fKeepBox->IsOn();
      for (size_t i = 0; i < fClusters.size(); ++i) {
         auto p = fClusters[i].GetPosition();
         bool in = co(p, ax[pane]) >= xlo && co(p, ax[pane]) <= xhi &&
                   co(p, ay[pane]) >= ylo && co(p, ay[pane]) <= yhi;
         if (!in) { if (keepMode) { if (fKeep[i]) { fKeep[i] = 0; ++n; } } continue; }
         if (!keepMode && fKeep[i]) { fKeep[i] = 0; ++n; }
      }
      fLog->AddLine(Form("  box on pane %d [%.0f,%.0f]x[%.0f,%.0f] -> %s %d clusters, %d now selected",
                         pane, xlo, xhi, ylo, yhi, keepMode ? "kept only, dropped" : "dropped",
                         n, (int)std::count(fKeep.begin(), fKeep.end(), 1)));
      fLog->ShowBottom();
      Refit();
   }
   void OnResetSel() { fKeep.assign(fClusters.size(), 1); fI0->SetNumber(0);
                       fI1->SetNumber(fClusters.empty() ? 0 : fClusters.size() - 1);
                       fLog->AddLine("  selection reset: all clusters"); Refit(); }
   void OnInvertSel() { for (auto &k : fKeep) k = !k;
                        fLog->AddLine(Form("  selection inverted -> %d clusters",
                                           (int)std::count(fKeep.begin(), fKeep.end(), 1))); Refit(); }
   /// Drop clusters whose charge is below the box value -- the other selection people actually
   /// want on a stopping track, where the tail is a handful of tiny clusters.
   void OnChargeCut()
   {
      double qmin = fQMin->GetNumber(); int n = 0;
      for (size_t i = 0; i < fClusters.size(); ++i)
         if (fKeep[i] && fClusters[i].GetCharge() < qmin) { fKeep[i] = 0; ++n; }
      fLog->AddLine(Form("  charge < %.0f -> dropped %d, %d left", qmin, n,
                         (int)std::count(fKeep.begin(), fKeep.end(), 1)));
      Refit();
   }

   /// HARD exit. Never gApplication->Terminate() -- see the file header.
   void Quit() { gSystem->Exit(0); }

private:
   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1500, 780);
      main->SetWindowName("a2091 15C+d  --  live single-track genfit refit");
      main->Connect("CloseWindow()", "DpFitGui", this, "Quit()");
      main->DontCallClose();

      auto *top = new TGHorizontalFrame(main);
      auto add = [&](TGCompositeFrame *p, TGWindow *w, int l = 4) {
         p->AddFrame((TGFrame *)w, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, l, 2, 4, 4)); };
      add(top, new TGLabel(top, "run"));
      fRunEntry = new TGTextEntry(top, fRunTag.Data()); fRunEntry->Resize(150, 22); add(top, fRunEntry, 1);
      add(top, new TGLabel(top, "entry"));
      fEntryNum = new TGNumberEntry(top, fEntry, 8, -1, TGNumberFormat::kNESInteger); add(top, fEntryNum, 1);
      add(top, new TGLabel(top, "trk"));
      fTidNum = new TGNumberEntry(top, fTid, 3, -1, TGNumberFormat::kNESInteger); add(top, fTidNum, 1);
      auto *bLoad = new TGTextButton(top, "load"); bLoad->Connect("Clicked()", "DpFitGui", this, "OnReload()"); add(top, bLoad);
      auto *bPrev = new TGTextButton(top, "< prev"); bPrev->Connect("Clicked()", "DpFitGui", this, "OnPrev()"); add(top, bPrev, 12);
      auto *bNext = new TGTextButton(top, "next >"); bNext->Connect("Clicked()", "DpFitGui", this, "OnNext()"); add(top, bNext, 1);
      main->AddFrame(top, new TGLayoutHints(kLHintsExpandX));

      auto *par = new TGHorizontalFrame(main);
      add(par, new TGLabel(par, "ordering"));
      fOrder = new TGComboBox(par); fOrder->AddEntry("z-sort (production)", 0);
      fOrder->AddEntry("cluster order", 1); fOrder->AddEntry("arc walk", 2);
      fOrder->AddEntry("vertex-first", 3);
      fOrder->Select(0, kFALSE); fOrder->Resize(150, 22); add(par, fOrder, 1);
      add(par, new TGLabel(par, "measSigma"));
      fSigma = new TGNumberEntry(par, 4.0, 5, -1, TGNumberFormat::kNESRealOne); add(par, fSigma, 1);
      add(par, new TGLabel(par, "iter"));
      fMinIt = new TGNumberEntry(par, 2, 3, -1, TGNumberFormat::kNESInteger); add(par, fMinIt, 1);
      fMaxIt = new TGNumberEntry(par, 5, 3, -1, TGNumberFormat::kNESInteger); add(par, fMaxIt, 1);
      add(par, new TGLabel(par, "clusters"));
      fI0 = new TGNumberEntry(par, 0, 5, -1, TGNumberFormat::kNESInteger); add(par, fI0, 1);
      fI1 = new TGNumberEntry(par, 9999, 5, -1, TGNumberFormat::kNESInteger); add(par, fI1, 1);
      fMatFx = new TGCheckButton(par, "matFX"); fMatFx->SetOn(); add(par, fMatFx, 10);
      fCatima = new TGCheckButton(par, "CATIMA"); fCatima->SetOn(); add(par, fCatima, 2);
      fBackFix = new TGCheckButton(par, "bwSeedFix"); fBackFix->SetOn(); add(par, fBackFix, 2);
      fSpySeed = new TGCheckButton(par, "spyralSeed"); add(par, fSpySeed, 2);
      main->AddFrame(par, new TGLayoutHints(kLHintsExpandX));

      auto *sel = new TGHorizontalFrame(main);
      add(sel, new TGLabel(sel, "drag a box on any panel:"));
      fKeepBox = new TGCheckButton(sel, "box KEEPS (else drops)"); fKeepBox->SetOn(); add(sel, fKeepBox, 2);
      auto *bRes = new TGTextButton(sel, "reset selection"); bRes->Connect("Clicked()", "DpFitGui", this, "OnResetSel()"); add(sel, bRes, 10);
      auto *bInv = new TGTextButton(sel, "invert"); bInv->Connect("Clicked()", "DpFitGui", this, "OnInvertSel()"); add(sel, bInv, 2);
      add(sel, new TGLabel(sel, "drop charge <"), 12);
      fQMin = new TGNumberEntry(sel, 0, 7, -1, TGNumberFormat::kNESInteger); add(sel, fQMin, 1);
      auto *bQ = new TGTextButton(sel, "apply"); bQ->Connect("Clicked()", "DpFitGui", this, "OnChargeCut()"); add(sel, bQ, 2);
      fSelLbl = new TGLabel(sel, "selected: -"); add(sel, fSelLbl, 14);
      main->AddFrame(sel, new TGLayoutHints(kLHintsExpandX));

      auto *par2 = new TGHorizontalFrame(main);
      par = par2;
      auto *bFit = new TGTextButton(par, "  REFIT  "); bFit->Connect("Clicked()", "DpFitGui", this, "OnRefit()"); add(par, bFit, 12);
      add(par, new TGLabel(par, "longest chi2<"), 10);
      // 0.1, the production FIND threshold. NOT the chi2<5 production CUT: measSigma 4 mm against
      // ~0.6 mm residuals puts the median chi2/ndf at 0.09 (1.48 backward), so a limit of 5 would
      // accept essentially everything and the button would stop at 100 % on every track.
      fC2Max = new TGNumberEntry(par, 0.1, 5, -1, TGNumberFormat::kNESRealThree); add(par, fC2Max, 1);
      auto *bLong = new TGTextButton(par, "find"); bLong->Connect("Clicked()", "DpFitGui", this, "OnLongest()"); add(par, bLong, 1);
      fVertexFirst = new TGCheckButton(par, "production order (z-sort)"); fVertexFirst->SetOn();
      fVertexFirst->Connect("Clicked()", "DpFitGui", this, "OnReorder()"); add(par, fVertexFirst, 8);
      fDirPhys = new TGCheckButton(par, "dir from physics");
      fDirPhys->Connect("Clicked()", "DpFitGui", this, "OnReorder()"); add(par, fDirPhys, 2);
      // VIEW-ONLY toggles: they call Draw(), never Refit(), so switching a layer on cannot change
      // the fit you are looking at. "hits" shows the raw AtHit cloud under the clusters; "order"
      // joins the clusters in STORED ARRAY ORDER, which is the thing to judge -- a good ordering
      // draws a single clean stroke along the track, a broken one draws long jumps across it.
      fShowHits = new TGCheckButton(par, "hits");
      fShowHits->Connect("Clicked()", "DpFitGui", this, "OnRedraw()"); add(par, fShowHits, 8);
      fShowOrder = new TGCheckButton(par, "order");
      fShowOrder->Connect("Clicked()", "DpFitGui", this, "OnRedraw()"); add(par, fShowOrder, 2);
      auto *bAll = new TGTextButton(par, "all 3 orderings"); bAll->Connect("Clicked()", "DpFitGui", this, "OnAllModes()"); add(par, bAll, 2);
      auto *bPng = new TGTextButton(par, "PNG"); bPng->Connect("Clicked()", "DpFitGui", this, "OnSavePng()"); add(par, bPng, 2);
      auto *bQuit = new TGTextButton(par, "Quit"); bQuit->Connect("Clicked()", "DpFitGui", this, "Quit()"); add(par, bQuit, 2);
      main->AddFrame(par, new TGLayoutHints(kLHintsExpandX));

      fCanvas = new TRootEmbeddedCanvas("c", main, 1480, 430);
      fCanvas->GetCanvas()->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", "DpFitGui", this,
                                    "HandleEvent(Int_t,Int_t,Int_t,TObject*)");
      main->AddFrame(fCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 4));
      fLog = new TGTextView(main, 1480, 190);
      main->AddFrame(fLog, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 4));
      fLog->AddLine("Each REFIT appends one line, so a scan over parameters stays on screen and can be compared.");
      fLog->AddLine("NO CONTROL TRACK IS ESTABLISHED FOR C15d YET -- verify one fit against fit_batch.sh before trusting a number here.");

      main->MapSubwindows(); main->Resize(main->GetDefaultSize()); main->MapWindow();
   }

   TString fRunTag, fGfDir, fEloss;
   /// Which fit file the CLUSTERS are read from -- <run>_genfit_<fSpecies>.root. The refit itself
   /// is always a proton (pdg 2212 below), so reading the deuteron-gated file means refitting
   /// those clusters under the proton hypothesis. Point this at "p" once a proton pass exists.
   TString fSpecies{"d"};
   Long64_t fEntry; int fTid; int fIdx = 0;
   TFile *fFile = nullptr;
   AtTrack fTrack;
   PraTheta fPra;   //<! AtPRA's arc-vs-z construction, recomputed for the theta pane
   //<! THE FIT AS STORED IN THE PRODUCTION FILE -- i.e. what FIND actually chose. The GUI's own
   //<! Refit() is a LIVE full-length fit and never applies FIND, so without this the window shows
   //<! a fit the production never made, identically whichever gfDir it is pointed at.
   std::vector<std::array<double,3>> fStoredPts;
   double fStoredKE{-1}, fStoredTh{-1}, fStoredC2{-1}, fStoredVxy{-1};
   int fStoredNpts{0};
   std::vector<AtHitCluster> fClusters;
   std::vector<std::array<double,3>> fFitPts;
   std::vector<std::tuple<TString,Long64_t,int>> fList;
   TRootEmbeddedCanvas *fCanvas = nullptr;
   TGTextView *fLog = nullptr;
   TGTextEntry *fRunEntry = nullptr;
   TGNumberEntry *fEntryNum=nullptr, *fTidNum=nullptr, *fSigma=nullptr, *fMinIt=nullptr,
                 *fMaxIt=nullptr, *fI0=nullptr, *fI1=nullptr;
   TGComboBox *fOrder = nullptr;
   TGCheckButton *fKeepBox = nullptr, *fDirPhys = nullptr, *fVertexFirst = nullptr;
   TGCheckButton *fShowHits = nullptr, *fShowOrder = nullptr;
   TGNumberEntry *fC2Max = nullptr;
   double fLastC2 = -1;
   TGNumberEntry *fQMin = nullptr;
   TGLabel *fSelLbl = nullptr;
   std::vector<char> fKeep;
   std::vector<int> fUsedIdx;
   TVirtualPad *fPads[4] = {nullptr, nullptr, nullptr};
   double fDragX0 = 0, fDragY0 = 0; bool fDragging = false; int fDragPane = -1;
   TGCheckButton *fMatFx=nullptr, *fCatima=nullptr, *fBackFix=nullptr, *fSpySeed=nullptr;
};

void gui_fit_C15d(TString runTag = "run_0013", Long64_t entry = 0, int tid = 0,
                  TString gfDir = "/home/yassid/C15d_fit/", TString listFile = "",
                  // Which fit file the CLUSTERS are read from. The refit below is always a PROTON
                  // (pdg 2212), so reading a deuteron-gated file means refitting those clusters
                  // under the proton hypothesis -- deliberate, and the only C15d fits that exist
                  // today. Point this at "p" once a proton pass has been run.
                  TString species = "d")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gStyle->SetOptStat(0);
   new DpFitGui(runTag, entry, tid, gfDir, listFile, species);
   printf("\n  window is open; Quit or [X] closes it.\n");
   // NO gApplication->Run() here -- see the file header.
}
