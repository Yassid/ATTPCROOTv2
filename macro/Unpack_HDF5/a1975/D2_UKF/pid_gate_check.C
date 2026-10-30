/// @file pid_gate_check.C
/// @brief GATE-LEVEL gain of a PID estimator, on KINEMATICALLY certain species (no estimator grades
///        itself): deuterons = on the 16C(d,d) elastic locus under both geometries (theta ~78-86 deg);
///        protons = BACKWARD (theta_lab > 100 deg under both geometries) -- elastic deuterons cannot go
///        beyond 90 deg and (d,t) tritons stay forward, so backward light particles are (d,p) protons.
///        For each estimator, gates of half-width w in ln(sqrt(dE/dx)*Brho) around its own peak:
///        efficiency on the true species and contamination by the other one, scanned in w; then the
///        efficiency at FIXED contamination (1 % and 3 %).
///        vars = "name:sqrtdEdx:brho;..."
///   root -l -b -q 'pid_gate_check.C("/mnt/f/a1975/pid3d_final/pv_all.root","Spyral:sqa:br;3D:s0:b0")'
#include "pid_dd_referee.C"
void pid_gate_check(TString in, TString vars, double Eb = 184.25)
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   int va, v3, nh;
   double br, b3, pol, th3, sp3, sq3h;
   w->SetBranchAddress("va", &va); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("nh", &nh);
   w->SetBranchAddress("br", &br); w->SetBranchAddress("b3", &b3); w->SetBranchAddress("pol", &pol);
   w->SetBranchAddress("th3", &th3); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("sq3h", &sq3h);
   std::vector<Long64_t> D, Pr;
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      if (!va || !v3 || sp3 >= 6 || nh < 30 || sq3h <= 0) continue;
      double thS = 180 - pol;
      if (std::abs(exDD(Eb, thS * M_PI / 180, br)) < 1.5 && std::abs(exDD(Eb, th3 * M_PI / 180, b3)) < 1.5) D.push_back(i);
      else if (thS > 100 && th3 > 100) Pr.push_back(i);
   }
   printf("kinematic samples: deuterons (elastic locus) %zu, protons (backward > 100 deg) %zu\n", D.size(), Pr.size());
   TString base = "va&&v3&&sp3<6&&nh>=30&&sq3h>0";
   for (auto *o : *vars.Tokenize(";")) {
      auto *fl = ((TObjString *)o)->GetString().Tokenize(":");
      TString nm = ((TObjString *)fl->At(0))->GetString();
      TTreeFormula fq("q", ((TObjString *)fl->At(1))->GetString(), w), fb("b", ((TObjString *)fl->At(2))->GetString(), w),
         fs("s", base, w);
      std::vector<double> all;
      for (Long64_t i = 0; i < w->GetEntries(); ++i) {
         w->GetEntry(i);
         if (fs.EvalInstance() && fq.EvalInstance() > 0) all.push_back(fq.EvalInstance());
      }
      std::sort(all.begin(), all.end());
      double medq = all[all.size() / 2];
      auto lnm = [&](Long64_t i) {
         w->GetEntry(i);
         double q = fq.EvalInstance(), b = fb.EvalInstance();
         return (q > 0 && b > 0) ? std::log(q / medq * b) : -99.0;
      };
      std::vector<double> ld, lp;
      for (auto i : D) ld.push_back(lnm(i));
      for (auto i : Pr) lp.push_back(lnm(i));
      auto medOf = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; };
      double cd = medOf(ld), cp = medOf(lp);
      auto frac = [](const std::vector<double> &v, double c, double hw) {
         long k = 0;
         for (double x : v) k += std::abs(x - c) < hw;
         return double(k) / v.size();
      };
      printf("\n%s: deuteron peak %.3f, proton peak %.3f (separation %.3f)\n", nm.Data(), cd, cp, cd - cp);
      printf("  gate half-width | d-gate: eff(d)  cont(p) | p-gate: eff(p)  cont(d)\n");
      // scan; remember best efficiency at fixed contamination
      double effD1 = 0, effD3 = 0, effP1 = 0, effP3 = 0;
      for (double hw = 0.05; hw <= 0.60001; hw += 0.025) {
         double eD = frac(ld, cd, hw), cP = frac(lp, cd, hw), eP = frac(lp, cp, hw), cD = frac(ld, cp, hw);
         if (cP <= 0.01) effD1 = std::max(effD1, eD);
         if (cP <= 0.03) effD3 = std::max(effD3, eD);
         if (cD <= 0.01) effP1 = std::max(effP1, eP);
         if (cD <= 0.03) effP3 = std::max(effP3, eP);
         if (std::fmod(hw + 1e-9, 0.1) < 0.026)
            printf("      %.3f       |   %5.1f %%  %5.2f %%   |   %5.1f %%  %5.2f %%\n", hw, 100 * eD, 100 * cP, 100 * eP, 100 * cD);
      }
      printf("  BEST at fixed contamination: d-gate eff %.1f %% (cont<=1 %%), %.1f %% (<=3 %%) | p-gate eff %.1f %% (<=1 %%), %.1f %% (<=3 %%)\n",
             100 * effD1, 100 * effD3, 100 * effP1, 100 * effP3);
   }
}
