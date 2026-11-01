/// @file pid_dt_slices_dp.C
/// @brief Deuteron/triton separation per Brho slice for several PID variants on the same tracks. In each slice the
///        ln-mass ln(sqrt(dE/dx)/median*Brho) - deuteron mode is histogrammed (log y); one row per slice, one
///        column per variant. Printed: triton-peak position, valley/triton-peak ratio (lower = better separated).
///        vars = "name:sqrtdEdx:brho;..." (tree "w", e.g. the merged pid_padgain_test_dp.C output).
///   root -l -b -q 'pid_dt_slices_dp.C("/mnt/f/a1975/pid3d_padgain/pv_gainnorm_full.root","Spyral:sqp_g:br;3D:s1_g:b1","plots/pid_dt_slices.png")'
void pid_dt_slices_dp(TString in, TString vars, TString png, TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30")
{
   gStyle->SetOptStat(0);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   std::vector<std::array<TString, 3>> V;
   for (auto *o : *vars.Tokenize(";")) {
      auto *fl = ((TObjString *)o)->GetString().Tokenize(":");
      V.push_back({((TObjString *)fl->At(0))->GetString(), ((TObjString *)fl->At(1))->GetString(),
                   ((TObjString *)fl->At(2))->GetString()});
   }
   // slices in the SPYRAL Brho (common track selection for every variant)
   const double sl[] = {0.2, 0.3, 0.4, 0.55, 0.7, 1.0};
   const int ns = 5;
   auto *c = new TCanvas("c", "c", 600 * V.size(), 320 * ns);
   c->Divide(V.size(), ns);
   for (size_t v = 0; v < V.size(); ++v) {
      auto &[nm, sq, br] = V[v];
      TString cut = sel + Form("&&%s>0&&%s>0", sq.Data(), br.Data());
      w->Draw(sq, cut, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      double med = x[x.size() / 2];
      TH1D hm("hm", "", 800, -6, 4);
      w->Draw(Form("log(%s/%g*%s)>>hm", sq.Data(), med, br.Data()), cut, "goff");
      hm.Smooth(2);
      double dmode = hm.GetBinCenter(hm.GetMaximumBin());
      printf("%s\n", nm.Data());
      for (int s = 0; s < ns; ++s) {
         c->cd(1 + v + s * V.size());
         gPad->SetLogy();
         gPad->SetGrid();
         auto *h = new TH1D(Form("h%zu_%d", v, s), Form("%s, B#rho_{Spyral} %.2f-%.2f;ln-mass - d mode;tracks", nm.Data(), sl[s], sl[s + 1]), 90, -1.0, 0.8);
         w->Draw(Form("log(%s/%g*%s)-%g>>%s", sq.Data(), med, br.Data(), dmode, h->GetName()),
                 cut + Form("&&br>%g&&br<%g", sl[s], sl[s + 1]), "goff");
         h->SetLineColor(v ? kRed + 1 : kBlue + 1);
         h->SetLineWidth(2);
         h->SetMinimum(0.8);
         h->Draw("hist");
         // triton peak: max in [0.25,0.65]; valley: min of 3-bin running mean in [0.1, tpos]
         auto mean3 = [&](int b) { return (h->GetBinContent(b - 1) + h->GetBinContent(b) + h->GetBinContent(b + 1)) / 3; };
         int b0 = h->FindBin(0.25), b1 = h->FindBin(0.65), bt = b0;
         for (int b = b0; b <= b1; ++b) if (mean3(b) > mean3(bt)) bt = b;
         int bv = h->FindBin(0.08);
         for (int b = h->FindBin(0.08); b < bt; ++b) if (mean3(b) < mean3(bv)) bv = b;
         double nd = h->Integral(h->FindBin(-0.15), h->FindBin(0.15)), nt = h->Integral(h->FindBin(0.25), h->FindBin(0.6));
         printf("  Brho %.2f-%.2f  n %6.0f  d(+-0.15) %6.0f  t-window(0.25-0.6) %5.0f  t-peak at %+.2f  valley at %+.2f  valley/peak %.2f\n",
                sl[s], sl[s + 1], h->GetEntries(), nd, nt, h->GetBinCenter(bt), h->GetBinCenter(bv),
                mean3(bt) > 0 ? mean3(bv) / mean3(bt) : -1);
      }
   }
   c->SaveAs(png);
}
