/// @file dmax_arbiter_dp.C
/// @brief Fit-free radius arbiter for the 3D-vs-Spyral Brho disagreement: d_max = the largest distance
///        a track reaches from the beam axis is the circle's DIAMETER for a track launched from the
///        axis (0.93 x 2R for a stopping spiral, measured on (d,t)). Compares R_3D and R_Spyral
///        (R = Brho sin(theta) / B) against d_max/2 on the tail list and on well-measured tracks.
///   root -l -b -q 'dmax_arbiter_dp.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016_guard.root")'
#include "pid_dd_referee.C"
void dmax_arbiter_dp(TString pidFile, TString reco = "/mnt/f/a1975/reco_d2_tau/run_0016_multifit_reco.root",
                     double Eb = 184.25)
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kFatal;
   TFile f(pidFile);
   auto *w = (TTree *)f.Get("w");
   TString S = "va&&v3&&sp3<6&&nh>=30&&sq3h>0";
   auto med = [&](const char *v) {
      w->Draw(v, S, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      return x[x.size() / 2];
   };
   double mS = med("sqa"), mH = med("sq3h");
   int va, v3, nh, tid;
   Long64_t ev;
   double br, b3, sqa, sq3h, pol, th3, sp3, R3;
   w->SetBranchAddress("ev", &ev); w->SetBranchAddress("tid", &tid); w->SetBranchAddress("va", &va);
   w->SetBranchAddress("v3", &v3); w->SetBranchAddress("nh", &nh); w->SetBranchAddress("br", &br);
   w->SetBranchAddress("b3", &b3); w->SetBranchAddress("sqa", &sqa); w->SetBranchAddress("sq3h", &sq3h);
   w->SetBranchAddress("pol", &pol); w->SetBranchAddress("th3", &th3); w->SetBranchAddress("sp3", &sp3);
   w->SetBranchAddress("R3", &R3);
   struct T { Long64_t ev; int tid; bool tail; double R3, RS; };
   std::vector<T> L;
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      if (!va || !v3 || sp3 >= 6 || nh < 30 || sq3h <= 0) continue;
      double thS = (180 - pol) * M_PI / 180;
      if (!(std::abs(exDD(Eb, thS, br)) < 1.5 && std::abs(exDD(Eb, th3 * M_PI / 180, b3)) < 1.5)) continue;
      double lS = std::log(sqa / mS * br), lH = std::log(sq3h / mH * b3);
      if (lS < -1.429 - 0.45) continue;
      bool tail = lH < -1.470 - 0.45;
      L.push_back({ev, tid, tail, R3, br * std::abs(std::sin(thS)) / 2.85 * 1000});
   }
   TFile fr(reco);
   auto *t = (TTree *)fr.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   std::vector<double> r3T, rST, r3G, rSG;
   int closer3 = 0, closerS = 0;
   for (auto &x : L) {
      t->GetEntry(x.ev);
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr.GetTrackID() != x.tid) continue;
         // d_max: robust (98th percentile) distance from the axis, so one stray hit cannot set it
         std::vector<double> rho;
         for (auto &h : tr.GetHitArray()) rho.push_back(std::hypot(h->GetPosition().X(), h->GetPosition().Y()));
         std::sort(rho.begin(), rho.end());
         double half = 0.5 * rho[rho.size() * 98 / 100] / 0.93;
         if (x.tail) {
            r3T.push_back(x.R3 / half); rST.push_back(x.RS / half);
            (std::abs(std::log(x.R3 / half)) < std::abs(std::log(x.RS / half)) ? closer3 : closerS)++;
         } else { r3G.push_back(x.R3 / half); rSG.push_back(x.RS / half); }
      }
   }
   auto q = [](const char *l, std::vector<double> v) {
      std::sort(v.begin(), v.end());
      int n = v.size();
      if (n > 3) printf("  %-34s n=%4d p16 %.2f med %.2f p84 %.2f\n", l, n, v[n * 16 / 100], v[n / 2], v[n * 84 / 100]);
   };
   printf("radius / (d_max / 2 / 0.93)   (1.00 = consistent with the track's own extent)\n");
   q("TAIL: R_3D", r3T); q("TAIL: R_Spyral", rST); q("well measured: R_3D", r3G); q("well measured: R_Spyral", rSG);
   printf("TAIL tracks where 3D is closer to d_max: %d, Spyral closer: %d\n", closer3, closerS);
}
