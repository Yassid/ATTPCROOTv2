/// @file gate_band_dt.C
/// @brief Is the loss in the (d,t) data gate real, or is the sim's triton band just off the gate? For truth-matched sim tritons at gain 100000 and gated data tracks, per Brho:
///        gate edges in sqrt(dE/dx), sim and data medians and spreads, the share of sim tracks below and above the gate, and the share inside after scaling sim sqrt(dE/dx)
///        to the data median (band alignment only; the band width stays the sim's).
///   root -l -b -q 'gate_band_dt.C()'
namespace gbd {
double pc(std::vector<double> v, double f) { if (v.empty()) return 0; std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, size_t(f * v.size()))]; }
}
void gate_band_dt(TString reco = "/mnt/f/a1975/sim_rebuild/reco_g100000/t_ex2_reco.root", TString reco2 = "/mnt/f/a1975/sim_rebuild/reco_g100000/t_gs_reco.root")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   auto gate = AtTools::AtParticleID::LoadJSON((TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1975/D2_UKF/pid/triton_d2_dv1104.json").Data());
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);
   const int NB = 9;
   const double bl[NB + 1] = {0.15, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.85, 1.05};
   std::vector<double> S[NB], D[NB];
   for (const char *run : {"0016", "0026", "0058", "0088"}) {
      TFile f(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", run));
      auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         if (t->GetEntry(e) <= 0 || !pa->GetEntries()) continue;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            auto r = spy.Estimate(tr);
            if (!r.valid || !gate.IsInside(r.sqrtdEdx, r.brho)) continue;
            for (int k = 0; k < NB; ++k) if (r.brho >= bl[k] && r.brho < bl[k + 1]) D[k].push_back(r.sqrtdEdx);
         }
      }
   }
   std::vector<double> Sb[NB];
   for (TString fn : {reco, reco2}) {
      TFile f(fn);
      auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      for (Long64_t e = 0; e < std::min(8000LL, t->GetEntries()); ++e) {
         if (t->GetEntry(e) <= 0 || !pa->GetEntries()) continue;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            int n = 0, nt = 0;
            for (auto &h : tr.GetHitArray()) { auto &ms = h->GetMCSimPointArray(); if (ms.empty()) continue; ++n; if (ms[0].A == 3 && ms[0].Z == 1) ++nt; }
            if (!n || nt <= n / 2) continue;
            auto r = spy.Estimate(tr);
            if (!r.valid) continue;
            for (int k = 0; k < NB; ++k) if (r.brho >= bl[k] && r.brho < bl[k + 1]) { S[k].push_back(r.sqrtdEdx); Sb[k].push_back(r.brho); }
         }
      }
   }
   printf("Brho bin    | gate lo..hi  | data N med | sim N med p16..p84 | med ratio s/d | in gate below above | in gate after alignment\n");
   for (int k = 0; k < NB; ++k) {
      const double bc = 0.5 * (bl[k] + bl[k + 1]);
      double lo = 1e9, hi = -1;
      for (double x = 0; x < 60; x += 0.02) if (gate.IsInside(x, bc)) { lo = std::min(lo, x); hi = std::max(hi, x); }
      if (S[k].size() < 20) continue;
      const double ms = gbd::pc(S[k], .5), md = gbd::pc(D[k], .5);
      long in = 0, be = 0, ab = 0, inA = 0;
      for (size_t i = 0; i < S[k].size(); ++i) {
         const double x = S[k][i];
         if (gate.IsInside(x, Sb[k][i])) ++in; else if (x < lo) ++be; else if (x > hi) ++ab;
         if (md > 0 && gate.IsInside(x * md / ms, Sb[k][i])) ++inA;
      }
      const double N = S[k].size();
      printf("%.2f-%.2f | %5.2f..%5.2f | %5zu %5.2f | %5zu %5.2f %5.2f..%5.2f | %5.2f | %5.1f %5.1f %5.1f | %5.1f\n", bl[k], bl[k + 1], lo, hi, D[k].size(), md, S[k].size(), ms, gbd::pc(S[k], .16), gbd::pc(S[k], .84), md > 0 ? ms / md : 0., 100 * in / N, 100 * be / N, 100 * ab / N, 100 * inA / N);
   }
}
