/// @file adaptive_fit_dp.C
/// @brief Adaptive cluster selection for genfit, as an OPTION, with the three candidate rules
///        run side by side on the same tracks so the choice is measured and not argued.
///
/// RULE 0  production        every cluster, exactly as prod_dp_catima.sh (the control)
/// RULE 1  pre-filter        drop clusters on CLUSTER properties only, never on the fit:
///                             charge < qFrac * median(charge)  OR  isolated in space
///                           Unbiased by construction: no fit result enters the decision.
/// RULE 2  trimmed           pre-filter, then fit -> drop clusters more than nSigma * rms from
///                           the trajectory -> refit, up to nPass times. Standard robust fitting;
///                           the criterion is "inconsistent with the trajectory".
/// RULE 3  largest-passing   the LARGEST prefix fraction whose chi2/ndf stays under a threshold.
///
/// WHY NOT "scan the fractions and keep the best chi2": chi2/ndf FALLS with fewer points -- on
/// run_0020/9789 it is 0.00 at 5% of the clusters and 0.07 at 80%. Best-chi2 therefore selects
/// 15 clusters of 298 and returns KE 2.490 against 2.807 from 80%, a 13% bias with a much larger
/// true uncertainty. Rule 3 takes the largest subset that still passes, i.e. the MOST information
/// rather than the least, which is why it is the fraction rule implemented here.
///
///   root -b -q 'adaptive_fit_dp.C("run_0020_multifit",9789,0)'          // one track, all rules
///   root -b -q 'adaptive_fit_dp.C("","",0,"list.txt","out.csv")'        // batch over a list
#include <algorithm>
#include <fstream>
#include <numeric>
#include <vector>

struct DpFitOut { double ke, th, c2n, ndf, rms; int nUsed, nPts; bool ok; };

static std::unique_ptr<EventFit::AtGenfitter> dpMakeFitter(TString elossPath, bool useClusterOrder)
{
   auto f = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 1.00782503207, 1,
                                                    elossPath.Data(), kFALSE, 2, 5);
   f->SetCatimaMaterial(kTRUE, kTRUE);
   f->SetCatimaELoss(kTRUE, kFALSE);
   f->SetELossHybrid(kTRUE, 6.61e-5);
   f->SetZPadPlane(1000.0);
   f->SetMeasSigma(4.0);
   f->SetSeedFromSpyral(kFALSE);
   f->SetMatEffectsFallback(kFALSE);
   f->SetThetaWindow(10.0, 170.0);
   f->SetBackwardSeedFix(kTRUE);
   f->SetBackExtrapToAxis(kTRUE);
   f->SetUseClusterOrder(useClusterOrder);
   f->Init();
   return f;
}

static DpFitOut dpFit(const AtTrack &proto, const std::vector<AtHitCluster> &use, TString elossPath,
                      std::vector<double> *residOut = nullptr)
{
   DpFitOut o{0,0,0,0,0,(int)use.size(),0,false};
   if (use.size() < 4) return o;
   auto fitter = dpMakeFitter(elossPath, kFALSE);
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
   o.ke = k.kineticEnergy; o.th = k.theta * TMath::RadToDeg();
   o.c2n = ndf > 0 ? c2 / ndf : 1e9; o.ndf = ndf; o.nPts = sp.size(); o.ok = true;
   double s2 = 0;
   if (residOut) residOut->assign(use.size(), 0.0);
   for (size_t i = 0; i < use.size(); ++i) {
      auto q = use[i].GetPosition(); double bd = 1e30;
      for (auto &m : sp) { double dz = q.Z() - (1000.0 - m.Z());
         double d = (q.X()-m.X())*(q.X()-m.X()) + (q.Y()-m.Y())*(q.Y()-m.Y()) + dz*dz;
         if (d < bd) bd = d; }
      s2 += bd;
      if (residOut) (*residOut)[i] = std::sqrt(bd);
   }
   o.rms = std::sqrt(s2 / use.size());
   return o;
}

/// RULE 1. Cluster properties only. `qFrac` of the track's own MEDIAN charge (a median, not a
/// mean: the Bragg peak is a factor-4 outlier and would drag a mean up and take real track with
/// it). `gapMM` rejects clusters whose nearest neighbour is far away -- the junk tail jumps 30 mm
/// where the track steps 1.7 mm, so it is isolated in a way real track never is.
static std::vector<AtHitCluster> dpPreFilter(const std::vector<AtHitCluster> &all, double qFrac,
                                             double gapMM, int *nDropQ = nullptr, int *nDropG = nullptr)
{
   std::vector<double> q;
   for (auto &c : all) q.push_back(c.GetCharge());
   std::vector<double> qs = q;
   std::nth_element(qs.begin(), qs.begin() + qs.size()/2, qs.end());
   const double qmed = qs[qs.size()/2];
   const double qcut = qFrac * qmed;
   int dq = 0, dg = 0;
   std::vector<AtHitCluster> out;
   for (size_t i = 0; i < all.size(); ++i) {
      if (q[i] < qcut) { ++dq; continue; }
      auto p = all[i].GetPosition();
      double nn = 1e30;
      for (size_t j = 0; j < all.size(); ++j) {
         if (j == i || q[j] < qcut) continue;
         auto r = all[j].GetPosition();
         double d = (p.X()-r.X())*(p.X()-r.X()) + (p.Y()-r.Y())*(p.Y()-r.Y()) + (p.Z()-r.Z())*(p.Z()-r.Z());
         if (d < nn) nn = d;
      }
      if (std::sqrt(nn) > gapMM) { ++dg; continue; }
      out.push_back(all[i]);
   }
   if (nDropQ) *nDropQ = dq;
   if (nDropG) *nDropG = dg;
   return out;
}

void adaptive_fit_dp(TString runTag = "run_0020_multifit", Long64_t entry = 9789, int tid = 0,
                     TString listFile = "", TString csvOut = "",
                     double qFrac = 0.25, double gapMM = 12.0, double nSigma = 3.0, int nPass = 3,
                     double c2Max = 5.0,
                     TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                     TString geoName = "ATTPC_D300torr_v2_geomanager.root",
                     TString eLossTable = "proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf = TFile::Open(dir + "/geometry/" + geoName); gf->Get("FAIRGeom"); }
   TString elossPath = dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable;

   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   if (listFile.Length()) {
      std::ifstream in(listFile.Data()); std::string r; Long64_t e; int t;
      while (in >> r >> e >> t) jobs.push_back({TString("run_") + r.c_str() + "_multifit", e, t});
   } else jobs.push_back({runTag, entry, tid});

   std::ofstream csv;
   if (csvOut.Length()) { csv.open(csvOut.Data());
      csv << "run,entry,tid,nclus,ke0,th0,c20,ke1,th1,c21,n1,ke2,th2,c22,n2,ke3,th3,c23,n3\n"; }
   if (!csvOut.Length())
      printf("\n  %-22s %6s | %8s %7s %9s | %8s %7s %9s %5s | %8s %7s %9s %5s | %8s %7s %9s %5s\n",
             "track","nclus","KE_prod","th","chi2","KE_pre","th","chi2","n","KE_trim","th","chi2","n",
             "KE_frac","th","chi2","n");

   TString openRun; TFile *f = nullptr; TTree *tree = nullptr; TClonesArray *te = nullptr;
   long done = 0;
   for (auto &j : jobs) {
      TString rt = std::get<0>(j); Long64_t en = std::get<1>(j); int ti = std::get<2>(j);
      if (rt != openRun) {
         if (f) f->Close();
         f = TFile::Open(gfDir + rt + "_genfitter_p.root");
         if (!f || f->IsZombie()) { f = nullptr; continue; }
         tree = (TTree *)f->Get("cbmsim"); te = nullptr;
         tree->SetBranchAddress("AtTrackingEvent", &te);
         openRun = rt;
      }
      tree->GetEntry(en);
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;
      auto trks = ev->GetTrackArray();
      AtTrack *src = nullptr;
      for (auto &tr : trks) if (tr.GetTrackID() == ti) src = &tr;
      if (!src) continue;
      std::vector<AtHitCluster> all = *src->GetHitClusterArray();

      DpFitOut o0 = dpFit(*src, all, elossPath);                             // rule 0
      int dq = 0, dg = 0;
      auto pre = dpPreFilter(all, qFrac, gapMM, &dq, &dg);
      DpFitOut o1 = dpFit(*src, pre, elossPath);                             // rule 1
      // rule 2: trim on residuals, starting from the pre-filtered set
      std::vector<AtHitCluster> cur = pre;
      DpFitOut o2 = o1;
      for (int pass = 0; pass < nPass; ++pass) {
         std::vector<double> res;
         DpFitOut o = dpFit(*src, cur, elossPath, &res);
         if (!o.ok) break;
         o2 = o;
         double cut = nSigma * o.rms;
         std::vector<AtHitCluster> keep;
         for (size_t i = 0; i < cur.size(); ++i) if (res[i] <= cut) keep.push_back(cur[i]);
         if (keep.size() == cur.size() || keep.size() < 8) break;
         cur = keep;
      }
      // rule 3: largest prefix fraction still under c2Max
      DpFitOut o3{}; o3.ok = false;
      for (int fr = 95; fr >= 20; fr -= 5) {
         int m = std::max(6, (int)std::lround(all.size() * fr / 100.0));
         if (m > (int)all.size()) continue;
         std::vector<AtHitCluster> use(all.begin(), all.begin() + m);
         DpFitOut o = dpFit(*src, use, elossPath);
         if (o.ok && o.c2n < c2Max) { o3 = o; break; }
      }
      if (csvOut.Length()) {
         csv << rt << "," << en << "," << ti << "," << all.size()
             << "," << o0.ke << "," << o0.th << "," << o0.c2n
             << "," << o1.ke << "," << o1.th << "," << o1.c2n << "," << o1.nUsed
             << "," << o2.ke << "," << o2.th << "," << o2.c2n << "," << o2.nUsed
             << "," << o3.ke << "," << o3.th << "," << o3.c2n << "," << o3.nUsed << "\n";
         if (++done % 20 == 0) printf("  %ld / %zu\n", done, jobs.size());
      } else {
         auto row = [](const DpFitOut &o) {
            return o.ok ? TString::Format("%8.3f %7.2f %9.2f %5d", o.ke, o.th, o.c2n, o.nUsed)
                        : TString("       -       -         -     -"); };
         printf("  %-22s %6zu | %s | %s | %s | %s\n", (rt + Form(" e%lld", en)).Data(), all.size(),
                TString::Format("%8.3f %7.2f %9.2f", o0.ke, o0.th, o0.c2n).Data(),
                row(o1).Data(), row(o2).Data(), row(o3).Data());
         printf("      pre-filter dropped %d on charge, %d isolated (of %zu)\n", dq, dg, all.size());
      }
   }
   if (f) f->Close();
   if (csvOut.Length()) { csv.close(); printf("  wrote %s (%ld tracks)\n", csvOut.Data(), done); }
}
