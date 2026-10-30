/// @file pid_angle_table.C
/// @brief Proton ln-mass (ln sqrt(dE/dx)*Brho, Z=1) per theta_lab bin for several PID estimators on the
///        SAME tracks, protons selected by the UNION of each estimator's window so no estimator is
///        graded on its own selection (grading Spyral inside Spyral's window hid its near-90 failure).
///        A good estimator is flat in theta and narrow in every bin.
///   root -l -b -q 'pid_angle_table.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root")'
void pid_angle_table(TString in, TString sel = "va&&vw&&v3&&sp3<6&&nh>=30&&sq3h>0", TString vars = "")
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   std::vector<std::array<std::string, 3>> V0;
   if (vars.Length()) {
      for (auto *o : *vars.Tokenize(";")) {
         auto *f = ((TObjString *)o)->GetString().Tokenize(":");
         V0.push_back({((TObjString *)f->At(0))->GetString().Data(), ((TObjString *)f->At(1))->GetString().Data(),
                       ((TObjString *)f->At(2))->GetString().Data()});
      }
   }
   std::vector<std::array<const char *, 3>> V = {
      {"Spyral amp", "sqa", "br"}, {"3D helix path", "sq3h", "b3"}};
   if (!V0.empty()) {
      V.clear();
      for (auto &v : V0) V.push_back({v[0].c_str(), v[1].c_str(), v[2].c_str()});
   }
   std::vector<TString> alias;
   TString uni = "(0";
   for (size_t k = 0; k < V.size(); ++k) {
      w->Draw(V[k][1], sel, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      TString a = Form("lm%zu", k);
      w->SetAlias(a, Form("log(%s/%g*%s)", V[k][1], x[x.size() / 2], V[k][2]));
      // window: around this estimator's own proton mode
      TH1D h(Form("hm%zu", k), "", 120, -3, 0);
      w->Draw(a + Form(">>hm%zu", k), sel, "goff");
      // proton = the peak ~ln2 below the deuteron (the global mode)
      double md = h.GetBinCenter(h.GetMaximumBin());
      h.GetXaxis()->SetRangeUser(md - 0.95, md - 0.45);
      double mp = h.GetBinCenter(h.GetMaximumBin());
      uni += Form("||(%s>%g&&%s<%g)", a.Data(), mp - 0.3, a.Data(), mp + 0.25);
      alias.push_back(a);
   }
   uni += ")";
   TString P = sel + "&&" + uni;
   printf("protons (union of windows): %lld\n  theta    n ", w->GetEntries(P));
   for (auto &v : V) printf(" | %-20.20s", v[0]);
   printf("\n");
   double tot[8] = {0};
   long ntot = 0;
   for (int a = 20; a < 160; a += 20) {
      TString c = P + Form("&&th3>=%d&&th3<%d", a, a + 20);
      long n = w->GetEntries(c);
      if (n < 8) continue;
      printf("  %3d-%3d %4ld", a, a + 20, n);
      for (size_t k = 0; k < V.size(); ++k) {
         w->Draw(alias[k], c, "goff");
         std::vector<double> x(w->GetV1(), w->GetV1() + n);
         std::sort(x.begin(), x.end());
         double hw = 0.5 * (x[n * 84 / 100] - x[n * 16 / 100]);
         printf(" | %+.3f [%.3f]      ", x[n / 2], hw);
         tot[k] += hw * hw * n;
      }
      ntot += n;
      printf("\n");
   }
   printf("  rms of half-widths (n-weighted):");
   for (size_t k = 0; k < V.size(); ++k) printf("  %s %.3f", V[k][0], std::sqrt(tot[k] / ntot));
   printf("\n");
}
