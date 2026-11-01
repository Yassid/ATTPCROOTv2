/// @file plot_pid_best_dp.C
/// @brief The adopted species-blind a1975 D2 PID (revised 2026-11-01 after the d + t simulations): 3D sqrt(dE/dx) with
///        run + pad gain (s1_g) against the gated Brho, the region chosen on the circle Brho bA:
///          bA <  b1 : circle + arc-length angle (bA)          -- best low-Brho d/t in the d + t sims
///          b1..b2   : spiral (shrinking-radius) fit, s.b21     -- best at 0.35-0.60 (sim and data)
///          bA >= b2 : Spyral (br)
///        Left: the PID plane, Brho on y, region boundaries dashed. Right: ln-mass ln(sqrt(dE/dx)/median * Brho) per
///        Brho range, each slice shifted so ITS deuteron mode is 0 (the three estimators have different scales), lines at
///        ln(m/m_d) for p, t and 4He (Z = 2 lands at ln(m_alpha/m_d) in this variable).
///        Inputs: pid_sbpid_dp.C merged tree + spiral_fit_dp.C data tree (same rows, friend "s").
///   root -l -b -q 'plot_pid_best_dp.C("plots/pid_best.png")'
void plot_pid_best_dp(TString png, double b1 = 0.35, double b2 = 0.60, TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root",
                      TString sp = "/mnt/f/a1975/pid3d_spiral/sp_all.root", TString spB = "s.b21",
                      TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&br>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("s=w", sp);
   // selector sum, not nested ?: (TTreeFormula returns nothing for the nested ternary)
   TString B = Form("((bA<%g)*bA+(bA>=%g&&bA<%g&&%s>0)*%s+(bA>=%g&&bA<%g&&%s<=0)*bA+(bA>=%g)*br)", b1, b1, b2, spB.Data(),
                    spB.Data(), b1, b2, spB.Data(), b2);
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double med = x[x.size() / 2];
   const double xhi = x[x.size() * 995 / 1000] * 1.1;
   TString L = Form("log(s1_g/%g*%s)", med, B.Data());
   auto *c = new TCanvas("c", "c", 1800, 900);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz(); gPad->SetRightMargin(0.12);
   auto *h2 = new TH2D("h2", Form("a1975 D2, runs 0016/0026/0058/0088: circle < %.2f, spiral %.2f-%.2f, Spyral above;"
                                  "#sqrt{dE/dx} (3D, run+pad gain) [arb.];B#rho [T m]", b1, b1, b2),
                       400, 0, xhi, 400, 0, 1.2);
   w->Draw(Form("%s:s1_g>>h2", B.Data()), sel, "goff");
   h2->Draw("colz");
   for (double b : {b1, b2}) { auto *l = new TLine(0, b, xhi, b); l->SetLineStyle(2); l->SetLineColor(kRed); l->Draw(); }
   c->cd(2);
   auto *pr = new TPad("pr", "", 0, 0, 1, 1); pr->Draw(); pr->cd(); pr->Divide(3, 2);
   const double bl[] = {0.10, 0.25, 0.35, 0.45, 0.60, 0.80, 1.00};
   const char *est[] = {"circle", "circle", "spiral", "spiral", "Spyral", "Spyral"};
   const double lm[] = {std::log(938.272 / 1875.613), 0, std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)};
   const char *ln[] = {"p", "d", "t", "^{4}He"};
   for (int b = 0; b < 6; ++b) {
      pr->cd(b + 1);
      gPad->SetLogy(); gPad->SetGridx();
      TString cut = sel + Form("&&%s>%g&&%s<%g", B.Data(), bl[b], B.Data(), bl[b + 1]);
      TH1D tmp(Form("t%d", b), "", 600, -4, 2);
      w->Draw(L + Form(">>t%d", b), cut, "goff");
      tmp.Smooth(2);
      const double d = tmp.GetBinCenter(tmp.GetMaximumBin());
      auto *h = new TH1D(Form("h_%d", b), Form("B#rho %.2f-%.2f T m (%s);ln mass - ln m_{d};tracks", bl[b], bl[b + 1], est[b]),
                         240, -1.2, 1.2);
      w->Draw(Form("%s-(%g)>>h_%d", L.Data(), d, b), cut, "goff");
      h->SetLineWidth(2); h->SetLineColor(kBlue + 2); h->SetMinimum(0.8); h->SetMaximum(3 * h->GetMaximum());
      h->Draw("hist");
      for (int k = 0; k < 4; ++k) {
         auto *l = new TLine(lm[k], 0.8, lm[k], h->GetMaximum()); l->SetLineStyle(2); l->SetLineColor(kRed + 1); l->Draw();
         auto *t = new TLatex(lm[k] + 0.02, h->GetMaximum() * 0.4, ln[k]); t->SetTextSize(0.06); t->SetTextColor(kRed + 1); t->Draw();
      }
   }
   c->SaveAs(png);
}
