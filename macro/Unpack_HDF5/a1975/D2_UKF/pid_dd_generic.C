/// @file pid_dd_generic.C
/// @brief Deuteron kinematic referee for ANY list of PID variants on a FIXED sample: tracks on the
///        16C(d,d) elastic locus under Spyral's AND the reference 3D geometry (br/pol, b3/th3). Reports
///        each variant's ln-mass half-width, far tail and proton-like fraction on that same sample.
///        vars = "name:sqrtdEdx:brho;..." (branch expressions of tree "w").
///   root -l -b -q 'pid_dd_generic.C("/mnt/f/a1975/pid3d_var/pv_all.root","Spyral:sqa:br;v0:s0:b0")'
#include "pid_dd_referee.C"
void pid_dd_generic(TString in, TString vars, double Eb = 184.25)
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   std::vector<std::array<TString, 3>> V;
   for (auto *o : *vars.Tokenize(";")) {
      auto *fl = ((TObjString *)o)->GetString().Tokenize(":");
      V.push_back({((TObjString *)fl->At(0))->GetString(), ((TObjString *)fl->At(1))->GetString(),
                   ((TObjString *)fl->At(2))->GetString()});
   }
   TString base = "va&&v3&&sp3<6&&nh>=30&&sq3h>0";
   // fixed referee sample
   int va, v3, nh;
   double br, b3, pol, th3, sp3, sq3h;
   w->SetBranchAddress("va", &va); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("nh", &nh);
   w->SetBranchAddress("br", &br); w->SetBranchAddress("b3", &b3); w->SetBranchAddress("pol", &pol);
   w->SetBranchAddress("th3", &th3); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("sq3h", &sq3h);
   std::vector<Long64_t> sample;
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      if (!va || !v3 || sp3 >= 6 || nh < 30 || sq3h <= 0) continue;
      if (std::abs(exDD(Eb, (180 - pol) * M_PI / 180, br)) < 1.5 && std::abs(exDD(Eb, th3 * M_PI / 180, b3)) < 1.5)
         sample.push_back(i);
   }
   printf("deuteron referee sample (fixed): %zu\n", sample.size());
   for (auto &v : V) {
      TTreeFormula fq("q", v[1], w), fb("b", v[2], w);
      // normalisation: median of this variant's sqrt(dE/dx) over the base selection
      std::vector<double> all;
      TTreeFormula fsel("s", base, w);
      for (Long64_t i = 0; i < w->GetEntries(); ++i) {
         w->GetEntry(i);
         if (fsel.EvalInstance() && fq.EvalInstance() > 0) all.push_back(fq.EvalInstance());
      }
      std::sort(all.begin(), all.end());
      double medq = all[all.size() / 2];
      std::vector<double> x;
      for (auto i : sample) {
         w->GetEntry(i);
         double q = fq.EvalInstance(), bb = fb.EvalInstance();
         if (q > 0 && bb > 0) x.push_back(std::log(q / medq * bb));
      }
      std::sort(x.begin(), x.end());
      int n = x.size();
      double m = x[n / 2], hw = 0.5 * (x[n * 84 / 100] - x[n * 16 / 100]);
      int far = 0, prot = 0;
      for (double y : x) { far += std::abs(y - m) > 0.3; prot += y < m - 0.45; }
      printf("  %-24s n=%6d  half-width %.4f  far(>0.3) %5.2f %%  proton-like %5.2f %%\n", v[0].Data(), n, hw,
             100.0 * far / n, 100.0 * prot / n);
   }
}
