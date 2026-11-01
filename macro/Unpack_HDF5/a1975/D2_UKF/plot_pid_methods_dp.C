/// @file plot_pid_methods_dp.C
/// @brief Side-by-side comparison of the a1975 D2 PID methods on the SAME tracks: Spyral and the 3D
///        estimator (AtHelix3D), each raw, with the per-run gain correction, and with run + pad-gain map.
///        Page 1: the six PID planes (sqrt(dE/dx) vs Brho, framed on the data). Page 2: the ln-mass
///        variable ln(sqrt(dE/dx)/median * Brho), shifted so the deuteron mode sits at 0, Spyral raw vs
///        Spyral best vs 3D best (log y, to show the tails).
///        Input: merged tree from pid_padgain_test_dp.C + pid_gain_norm.C (columns sqa/sqa_g/sqp_g, s0/s0_g/s1_g).
///   root -l -b -q 'plot_pid_methods_dp.C("/mnt/f/a1975/pid3d_padgain/pv_gainnorm_full.root","plots/pid_methods")'
void plot_pid_methods_dp(TString in, TString pngBase, TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   struct M { const char *name, *sq, *br; };
   const M m[6] = {{"Spyral, raw", "sqa", "br"},          {"Spyral, run gain", "sqa_g", "br"},
                   {"Spyral, run + pad gain", "sqp_g", "br"}, {"3D, raw", "s0", "b0"},
                   {"3D, run gain", "s0_g", "b0"},          {"3D, run + pad gain", "s1_g", "b1"}};
   auto quant = [&](const char *expr, double q) {
      w->Draw(expr, sel + Form("&&%s>0", expr), "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      return x[size_t(q * (x.size() - 1))];
   };
   // page 1: planes, x framed per estimator family on its raw 99.5th percentile
   auto *c = new TCanvas("c", "c", 1800, 1100);
   c->Divide(3, 2);
   for (int i = 0; i < 6; ++i) {
      c->cd(i + 1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.12);
      double xmax = 1.15 * quant(m[i].sq, 0.995);
      auto *h = new TH2D(Form("h%d", i), Form("%s;#sqrt{dE/dx} [arb.];B#rho [T m]", m[i].name), 300, 0, xmax, 300, 0,
                         1.0);
      w->Draw(Form("%s:%s>>h%d", m[i].br, m[i].sq, i), sel + Form("&&%s>0", m[i].sq), "goff");
      h->Draw("colz");
   }
   c->SaveAs(pngBase + "_planes.png");
   // page 2: ln-mass, deuteron mode at 0
   auto *c2 = new TCanvas("c2", "c2", 1000, 700);
   gPad->SetLogy();
   gPad->SetGrid();
   const int pick[3] = {0, 2, 5};
   const int col[3] = {kGray + 2, kBlue + 1, kRed + 1};
   auto *lg = new TLegend(0.12, 0.70, 0.45, 0.88);
   double ymax = 0;
   std::vector<TH1D *> hs;
   for (int k = 0; k < 3; ++k) {
      const M &mm = m[pick[k]];
      double med = quant(mm.sq, 0.5);
      TString cut = sel + Form("&&%s>0&&%s>0", mm.sq, mm.br);
      // deuteron = the dominant peak, located on a wide fine histogram, then subtracted in the expression
      TH1D hm(Form("lm%d", k), "", 800, -6, 4);
      w->Draw(Form("log(%s/%g*%s)>>lm%d", mm.sq, med, mm.br, k), cut, "goff");
      hm.Smooth(2);
      double dmode = hm.GetBinCenter(hm.GetMaximumBin());
      auto *hs_ = new TH1D(Form("ls%d", k), ";ln(#sqrt{dE/dx}/median #times B#rho) - deuteron mode;tracks", 110, -1.4, 0.8);
      w->Draw(Form("log(%s/%g*%s)-%g>>ls%d", mm.sq, med, mm.br, dmode, k), cut, "goff");
      hs_->SetLineColor(col[k]);
      hs_->SetLineWidth(k ? 2 : 3);
      ymax = std::max(ymax, hs_->GetMaximum());
      hs.push_back(hs_);
      lg->AddEntry(hs_, mm.name, "l");
   }
   for (int k = 0; k < 3; ++k) {
      hs[k]->SetMaximum(3 * ymax);
      hs[k]->SetMinimum(1);
      hs[k]->GetXaxis()->SetRangeUser(-1.4, 0.8);
      hs[k]->Draw(k ? "hist same" : "hist");
   }
   lg->Draw();
   TLatex t;
   t.SetTextSize(0.035);
   t.DrawLatexNDC(0.60, 0.84, "deuterons");
   t.DrawLatexNDC(0.25, 0.55, "protons");
   c2->SaveAs(pngBase + "_lnmass.png");
}
