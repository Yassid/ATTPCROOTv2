/// @file pid_proton_fixed.C
/// @brief Proton ln-mass width per theta_lab bin on a FIXED, method-symmetric proton sample: tracks that
///        BOTH reference estimators (Spyral and the 3D helix-path PID) put inside a generous proton
///        window (+-0.35 around each one's proton mode). Unlike the union selection, the sample does
///        not change when variants are added, so variants can be tuned against each other.
///        vars = "name:sqrtdEdx:brho;..." (any branch expressions of the tree "w").
///   root -l -b -q 'pid_proton_fixed.C("/mnt/f/a1975/pid3d_multi/pw_all.root","3D:sq3h:b3;Spyral:sqa:br")'
void pid_proton_fixed(TString in, TString vars, TString sel = "va&&v3&&sp3<6&&nh>=30&&sq3h>0")
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   auto lnm = [&](const char *sq, const char *br, TString &alias, const char *nm) {
      w->Draw(sq, sel, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      alias = Form("lm_%s", nm);
      w->SetAlias(alias, Form("log((%s)/%g*(%s))", sq, x[x.size() / 2], br));
      TH1D h(Form("h_%s", nm), "", 120, -3, 0);
      w->Draw(alias + Form(">>h_%s", nm), sel, "goff");
      double md = h.GetBinCenter(h.GetMaximumBin());
      h.GetXaxis()->SetRangeUser(md - 0.95, md - 0.45);
      return h.GetBinCenter(h.GetMaximumBin()); // proton mode
   };
   TString aS, a3;
   double pS = lnm("sqa", "br", aS, "refS"), p3 = lnm("sq3h", "b3", a3, "ref3");
   TString P = sel + Form("&&abs(%s-(%g))<0.35&&abs(%s-(%g))<0.35", aS.Data(), pS, a3.Data(), p3);
   std::vector<std::pair<TString, TString>> V;
   int iv = 0;
   for (auto *o : *vars.Tokenize(";")) {
      auto *fl = ((TObjString *)o)->GetString().Tokenize(":");
      TString nm = ((TObjString *)fl->At(0))->GetString(), al;
      TString sq = ((TObjString *)fl->At(1))->GetString(), br = ((TObjString *)fl->At(2))->GetString();
      lnm(sq, br, al, Form("v%d", iv++)); // alias names must be plain identifiers
      V.push_back({nm, al});
   }
   printf("FIXED proton sample (both references call it a proton): %lld\n  theta     n ", w->GetEntries(P));
   for (auto &v : V) printf("| %-16.16s ", v.first.Data());
   printf("\n");
   std::vector<double> acc(V.size(), 0);
   long ntot = 0;
   for (int a = 20; a < 160; a += 20) {
      TString c = P + Form("&&th3>=%d&&th3<%d", a, a + 20);
      long n = w->GetEntries(c);
      if (n < 20) continue;
      printf("  %3d-%3d %5ld ", a, a + 20, n);
      for (size_t k = 0; k < V.size(); ++k) {
         w->Draw(V[k].second, c, "goff");
         std::vector<double> x(w->GetV1(), w->GetV1() + n);
         std::sort(x.begin(), x.end());
         double hw = 0.5 * (x[n * 84 / 100] - x[n * 16 / 100]);
         acc[k] += hw * hw * n;
         printf("| %.4f           ", hw);
      }
      ntot += n;
      printf("\n");
   }
   printf("  rms (n-weighted)  ");
   for (size_t k = 0; k < V.size(); ++k) printf("| %.4f           ", std::sqrt(acc[k] / ntot));
   printf("\n");
}
