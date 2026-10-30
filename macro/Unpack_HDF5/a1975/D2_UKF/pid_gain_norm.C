/// @file pid_gain_norm.C
/// @brief PER-RUN GAIN NORMALISATION of the PID dE/dx, self-calibrated on the data.
///
/// The deuteron band of a1975 D2 drifts through the campaign: its ln(sqrt(dE/dx)*Brho) mode goes from
/// 1.79 (run 0016) to 1.67 (run 0097) -- the gain falls by ~20 %. Pooling runs smears every band and a
/// single gate is off-centre for some runs (a run-dependent efficiency). The fix: per run, measure the
/// band position on a FIXED-kinematics slice (theta_lab 78-86 deg, Brho 0.15-0.35 T m: dominated by
/// elastic deuterons, so momentum does not enter) and scale sqrt(dE/dx) so it lands at a common
/// reference. Done for each estimator with its own factors (gain affects both).
/// Writes a merged tree with `run`, the factors, and corrected sqrt(dE/dx) columns `<name>_g`.
///
///   root -l -b -q 'pid_gain_norm.C("/mnt/f/a1975/pid3d_final/pv_","0016,0020,...","out.root","sqa:br;s0:b0")'
void pid_gain_norm(TString prefix, TString runsCSV, TString out, TString cols = "sqa:br;s0:b0",
                   TString calSel = "va&&v3&&sp3<6&&nh>=30&&th3>78&&th3<86")
{
   std::vector<std::pair<TString, TString>> C;
   for (auto *o : *cols.Tokenize(";")) {
      auto *f = ((TObjString *)o)->GetString().Tokenize(":");
      C.push_back({((TObjString *)f->At(0))->GetString(), ((TObjString *)f->At(1))->GetString()});
   }
   // pass 1: per-run band position (mode of a fine histogram, robust) for each estimator
   std::vector<TString> runs;
   for (auto *o : *runsCSV.Tokenize(",")) runs.push_back(((TObjString *)o)->GetString());
   std::map<TString, std::vector<double>> mode; // per column, per run
   for (auto &r : runs) {
      TFile f(prefix + r + ".root");
      auto *w = (TTree *)f.Get("w");
      for (auto &[sq, br] : C) {
         TH1D h("h", "", 400, -2, 6);
         w->Draw(Form("log(%s*%s)>>h", sq.Data(), br.Data()),
                 calSel + Form("&&%s>0.15&&%s<0.35&&%s>0", br.Data(), br.Data(), sq.Data()), "goff");
         // mode refined by the mean of the +-0.05 window around the peak bin
         int ib = h.GetMaximumBin();
         double s = 0, sw = 0;
         for (int b = ib - 6; b <= ib + 6; ++b) { s += h.GetBinCenter(b) * h.GetBinContent(b); sw += h.GetBinContent(b); }
         mode[sq].push_back(sw > 0 ? s / sw : 0);
      }
   }
   std::map<TString, double> ref;
   for (auto &[sq, br] : C) ref[sq] = mode[sq][0]; // reference = the first run
   printf("per-run band position (ln m) and applied sqrt(dE/dx) factor:\n");
   for (size_t i = 0; i < runs.size(); ++i) {
      printf("  run %s", runs[i].Data());
      for (auto &[sq, br] : C) printf("  %s %.3f (x%.3f)", sq.Data(), mode[sq][i], std::exp(ref[sq] - mode[sq][i]));
      printf("\n");
   }
   // pass 2: merged tree with corrected columns
   TFile fo(out, "RECREATE");
   TTree *mt = nullptr;
   int runNo;
   std::map<TString, double> g, corr;
   for (size_t i = 0; i < runs.size(); ++i) {
      TFile f(prefix + runs[i] + ".root");
      auto *w = (TTree *)f.Get("w");
      // bind the raw columns FIRST, so the output tree is linked to these addresses (binding them
      // after CloneTree/CopyAddresses left the output pointing at stale buffers: raw columns = garbage)
      std::map<TString, double> raw;
      for (auto &[sq, br] : C) w->SetBranchAddress(sq, &raw[sq]);
      fo.cd();
      if (!mt) {
         mt = w->CloneTree(0);
         mt->Branch("run", &runNo);
         for (auto &[sq, br] : C) {
            mt->Branch(sq + "_g", &corr[sq]);
            mt->Branch("g_" + sq, &g[sq]);
         }
      } else
         w->CopyAddresses(mt);
      runNo = runs[i].Atoi();
      for (Long64_t j = 0; j < w->GetEntries(); ++j) {
         w->GetEntry(j);
         for (auto &[sq, br] : C) {
            g[sq] = std::exp(ref[sq] - mode[sq][i]);
            corr[sq] = raw[sq] > 0 ? raw[sq] * g[sq] : raw[sq];
         }
         mt->Fill();
      }
      w->ResetBranchAddresses();
   }
   fo.cd();
   mt->Write();
   fo.Close();
   printf("wrote %s\n", out.Data());
}
