/// @file plot_pid_best_dp.C
/// @brief The adopted species-blind a1975 D2 PID (2026-11-01): 3D sqrt(dE/dx) with run + pad gain (s1_g) against the
///        gated Brho -- hit-refit circle + arc-length angle (bA) below bSwitch, Spyral (br) above. Left: the PID plane, Brho on y
///        (log z, framed on the data). Right: the ln-mass variable ln(sqrt(dE/dx)/median * Brho) per Brho range, deuteron
///        mode at 0, lines at ln(m/m_d) for p, t and 4He (Z = 2 lands at ln(m_alpha/m_d) in this variable).
///   root -l -b -q 'plot_pid_best_dp.C("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root","plots/pid_best.png")'
void plot_pid_best_dp(TString in, TString png, double bSwitch = 0.45,
                      TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&br>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   TString B = Form("(bA<%g?bA:br)", bSwitch);
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double med = x[x.size() / 2];
   const double yhi = x[x.size() * 995 / 1000] * 1.1;
   TString L = Form("log(s1_g/%g*%s)", med, B.Data());
   TH1D hd("hd", "", 600, -3, 3);
   w->Draw(L + ">>hd", sel + Form("&&%s>0.15&&%s<0.35", B.Data(), B.Data()), "goff");
   hd.Smooth(2);
   const double dmode = hd.GetBinCenter(hd.GetMaximumBin());
   auto *c = new TCanvas("c", "c", 1800, 900);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz(); gPad->SetRightMargin(0.12);
   auto *h2 = new TH2D("h2", Form("a1975 D2, runs 0016/0026/0058/0088: B#rho (hit-refit+arc below %.2f, Spyral above) "
                                  "vs 3D #sqrt{dE/dx} (run+pad gain);#sqrt{dE/dx} [arb.];B#rho [T m]", bSwitch),
                       400, 0, yhi, 400, 0, 1.2);
   w->Draw(Form("%s:s1_g>>h2", B.Data()), sel, "goff");
   h2->Draw("colz");
   auto *ls = new TLine(0, bSwitch, yhi, bSwitch); ls->SetLineStyle(2); ls->SetLineColor(kRed); ls->Draw();
   c->cd(2);
   auto *pr = new TPad("pr", "", 0, 0, 1, 1); pr->Draw(); pr->cd(); pr->Divide(2, 2);
   const double bl[] = {0.1, 0.3, 0.45, 0.6, 0.8};
   const double lm[] = {std::log(938.272 / 1875.613), 0, std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)};
   const char *ln[] = {"p", "d", "t", "^{4}He"};
   for (int b = 0; b < 4; ++b) {
      pr->cd(b + 1);
      gPad->SetLogy(); gPad->SetGridx();
      auto *h = new TH1D(Form("h_%d", b), Form("B#rho %.2f-%.2f T m (%s);ln mass - ln m_{d};tracks", bl[b], bl[b + 1],
                                              bl[b] < bSwitch ? "hit-refit+arc" : "Spyral B#rho"), 240, -1.2, 1.2);
      w->Draw(Form("%s-(%g)>>h_%d", L.Data(), dmode, b), sel + Form("&&%s>%g&&%s<%g", B.Data(), bl[b], B.Data(), bl[b + 1]), "goff");
      h->SetLineWidth(2); h->SetLineColor(kBlue + 2); h->SetMinimum(0.8); h->SetMaximum(3 * h->GetMaximum());
      h->Draw("hist");
      for (int k = 0; k < 4; ++k) {
         auto *l = new TLine(lm[k], 0.8, lm[k], h->GetMaximum()); l->SetLineStyle(2); l->SetLineColor(kRed + 1); l->Draw();
         auto *t = new TLatex(lm[k] + 0.02, h->GetMaximum() * 0.4, ln[k]); t->SetTextSize(0.05); t->SetTextColor(kRed + 1); t->Draw();
      }
   }
   c->SaveAs(png);
}
