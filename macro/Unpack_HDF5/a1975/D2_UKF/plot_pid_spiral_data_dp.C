/// @file plot_pid_spiral_data_dp.C
/// @brief Data (runs 0016/0026/0058/0088): the spiral (shrinking-radius) Brho against the circle (bA) and Spyral Brho,
///        all with the same 3D dE/dx s1_g (2026-11-01). Top: PID planes (Brho on y). Bottom: ln-mass per slice of each
///        method's own Brho, each shifted so its deuteron mode is 0, log y; lines at p, t, 4He.
///        Inputs: pid_sbpid_dp.C merged tree + spiral_fit_dp.C data tree (same row order, friend "s").
///   root -l -b -q 'plot_pid_spiral_data_dp.C("plots/pid_spiral_data.png")'
void plot_pid_spiral_data_dp(TString png, TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root",
                             TString sp = "/mnt/f/a1975/pid3d_spiral/sp_all.root", TString spB = "s.b20")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("s=w", sp);
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&br>0&&" + spB + ">0";
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   const int nm = 3;
   TString B[nm] = {"bA", "br", spB};
   const char *name[nm] = {"circle + arc angle (adopted)", "Spyral B#rho", "spiral (shrinking radius)"};
   const int col[nm] = {kBlue + 2, kRed + 1, kGreen + 2};
   auto *c = new TCanvas("c", "c", 2000, 1300);
   auto *top = new TPad("top", "", 0, 0.45, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.45); bot->Draw();
   top->Divide(nm, 1);
   for (int m = 0; m < nm; ++m) {
      top->cd(m + 1); gPad->SetLogz(); gPad->SetRightMargin(0.12);
      auto *h2 = new TH2D(Form("p%d", m), Form("%s;#sqrt{dE/dx} / median (3D);B#rho [T m]", name[m]), 300, 0, 2.6, 300, 0, 1.0);
      w->Draw(Form("%s:s1_g/%g>>p%d", B[m].Data(), md, m), sel, "goff");
      h2->Draw("colz");
   }
   const double bl[] = {0.15, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80};
   bot->Divide(6, 1);
   for (int b = 0; b < 6; ++b) {
      bot->cd(b + 1); gPad->SetLogy(); gPad->SetGridx();
      TH1D *h[nm];
      double ymax = 0;
      for (int m = 0; m < nm; ++m) {
         TString L = Form("log(s1_g/%g*%s)", md, B[m].Data());
         TString cut = sel + Form("&&%s>%g&&%s<%g", B[m].Data(), bl[b], B[m].Data(), bl[b + 1]);
         TH1D tmp(Form("t%d%d", b, m), "", 600, -4, 2);
         w->Draw(L + Form(">>t%d%d", b, m), cut, "goff");
         tmp.Smooth(2);
         double d = tmp.GetBinCenter(tmp.GetMaximumBin());
         h[m] = new TH1D(Form("h%d%d", b, m), Form("B#rho %.2f-%.2f T m;ln mass - ln m_{d};tracks", bl[b], bl[b + 1]), 160, -1.2, 1.2);
         w->Draw(Form("%s-(%g)>>h%d%d", L.Data(), d, b, m), cut, "goff");
         h[m]->SetLineColor(col[m]); h[m]->SetLineWidth(2);
         ymax = std::max(ymax, h[m]->GetMaximum());
      }
      h[0]->SetMaximum(3 * ymax); h[0]->SetMinimum(0.8);
      h[0]->Draw("hist");
      for (int m = 1; m < nm; ++m) h[m]->Draw("hist same");
      for (double v : {std::log(938.272 / 1875.613), std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)}) {
         auto *l = new TLine(v, 0.8, v, 3 * ymax); l->SetLineStyle(3); l->Draw();
      }
      if (b == 0) {
         auto *lg = new TLegend(0.12, 0.75, 0.75, 0.89);
         for (int m = 0; m < nm; ++m) lg->AddEntry(h[m], name[m], "l");
         lg->Draw();
      }
   }
   c->SaveAs(png);
}
