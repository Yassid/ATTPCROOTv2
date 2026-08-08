/// @file fresco_cmp_C14.C
/// @brief Compare the measured 14C(p,p') angular distributions against FRESCO DWBA shapes.
///
/// Neither fitter's elastic distribution can be trusted a priori -- UKF and GENFIT disagree
/// qualitatively between 60 and 140 deg -- so an independent prediction is the only arbiter.
/// FRESCO runs in normal kinematics (p on 14C) at the lab proton energy giving the same E_cm as
/// 14C at 161 MeV on a proton: E_lab(p) = 11.581 MeV. Koning-Delaroche 2003 proton global OMP,
/// computed in ~/a1954_C14_fresco/kd_params.py (NOT copied from the a2091 15C input, whose real
/// depth uses the neutron asymmetry term).
///
/// SHAPE ONLY: FRESCO is scaled by one factor per curve, fitted over a forward window where the
/// acceptance is well measured and the statistics are large. No luminosity is claimed.
///
/// Elastic data comes from the no-chi2 caches (the chi2 cut was shown to distort 20-40 deg);
/// inelastic from the chi2<5 caches (there the cut improves fitter agreement).
///
///   root -b -q 'fresco_cmp_C14.C()'

static TGraph *loadFresco(TString path, double scale = 1.0)
{
   auto *g = new TGraph();
   std::ifstream in(path.Data());
   if (!in) {
      printf("\033[1;31mcannot read %s\033[0m\n", path.Data());
      return nullptr;
   }
   double th, xs;
   int n = 0;
   while (in >> th >> xs)
      g->SetPoint(n++, th, xs * scale);
   return g;
}

/// least-squares factor that puts the HISTOGRAM onto the graph over [lo, hi], i.e. the number
/// the data must be multiplied by so that data*s ~ FRESCO. (The inner fit solves s' with
/// s'*FRESCO ~ data; the returned value is its reciprocal.)
static double fitScale(TH1D *h, TGraph *g, double lo, double hi)
{
   double sn = 0, sd = 0;
   for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double c = h->GetBinCenter(b), y = h->GetBinContent(b), e = h->GetBinError(b);
      if (c < lo || c > hi || y <= 0 || e <= 0)
         continue;
      double f = g->Eval(c);
      if (f <= 0)
         continue;
      sn += y * f / (e * e);
      sd += f * f / (e * e);
   }
   double sPrime = sd > 0 ? sn / sd : 1.0;
   return sPrime > 0 ? 1.0 / sPrime : 1.0;
}

void fresco_cmp_C14(TString frDir = "/home/yassid/a1954_C14_fresco/outputs", Double_t normLo = 20.0,
                    Double_t normHi = 50.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   auto grab = [&](const char *f, const char *obj, const char *as) -> TH1D * {
      TFile *fi = TFile::Open(here + "/plots/" + f);
      if (!fi || fi->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", f);
         return nullptr;
      }
      auto *h = (TH1D *)fi->Get(obj);
      if (!h) {
         printf("\033[1;31mno %s in %s\033[0m\n", obj, f);
         return nullptr;
      }
      auto *c = (TH1D *)h->Clone(as);
      c->SetDirectory(nullptr);
      fi->Close();
      return c;
   };
   // acceptance-corrected dsigma/dOmega (the sin(theta) factor is already divided out there)
   TH1D *elU = grab("angdist_gs_ukf_nc.root", "dsigma_dOmega", "elU");
   TH1D *elG = grab("angdist_gs_gf_nc.root", "dsigma_dOmega", "elG");
   TH1D *inU = grab("angdist_ex1_ukf.root", "dsigma_dOmega", "inU");
   TH1D *inG = grab("angdist_ex1_gf.root", "dsigma_dOmega", "inG");
   if (!elU || !elG || !inU || !inG)
      return;

   TGraph *frEl = loadFresco(frDir + "/p14C_el_161_dsdo.dat");
   TGraph *frEl155 = loadFresco(frDir + "/p14C_el_155_dsdo.dat");
   TGraph *frL1 = loadFresco(frDir + "/p14C_inel_161_6094_L1_dsdo_ex2.dat");
   TGraph *frL3 = loadFresco(frDir + "/p14C_inel_161_6728_L3_dsdo_ex2.dat");
   TGraph *frL2 = loadFresco(frDir + "/p14C_inel_161_7012_L2_dsdo_ex2.dat");
   if (!frEl || !frL1)
      return;

   const double sU = fitScale(elU, frEl, normLo, normHi);
   const double sG = fitScale(elG, frEl, normLo, normHi);
   printf("\nelastic normalisation over %.0f-%.0f deg:  UKF x %.4g   GENFIT x %.4g\n", normLo, normHi, sU, sG);

   printf("\n===== elastic: measured shape vs FRESCO (both scaled to FRESCO) =====\n");
   printf("  theta_cm   FRESCO 161  FRESCO 155 |   UKF/norm   ratio  |  GENFIT/norm   ratio\n");
   for (int b = 1; b <= elU->GetNbinsX(); ++b) {
      double c = elU->GetBinCenter(b);
      if (c < 15 || c > 155)
         continue;
      double f = frEl->Eval(c), f5 = frEl155->Eval(c);
      double u = elU->GetBinContent(b) * sU, g = elG->GetBinContent(b) * sG;
      printf("  %3.0f-%3.0f  %11.4g %11.4g | %10.4g  %6.2f  | %11.4g  %6.2f\n", elU->GetBinLowEdge(b),
             elU->GetBinLowEdge(b) + elU->GetBinWidth(b), f, f5, u, f > 0 ? u / f : 0, g, f > 0 ? g / f : 0);
   }

   TCanvas *c1 = new TCanvas("c1", "vs FRESCO", 1500, 950);
   c1->Divide(2, 2);
   auto styleH = [](TH1D *h, int col, int mk) {
      h->SetMarkerStyle(mk);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
      h->SetMarkerSize(1.2);
   };
   auto styleG = [](TGraph *g, int col, int ls) {
      g->SetLineColor(col);
      g->SetLineWidth(3);
      g->SetLineStyle(ls);
   };
   styleH(elU, kAzure + 2, 20);
   styleH(elG, kRed + 1, 21);
   styleH(inU, kAzure + 2, 20);
   styleH(inG, kRed + 1, 21);
   styleG(frEl, kBlack, 1);
   styleG(frEl155, kGray + 2, 2);
   styleG(frL1, kBlack, 1);
   styleG(frL3, kGreen + 3, 7);
   styleG(frL2, kViolet + 1, 5);

   // ---- panel 1: elastic, log
   c1->cd(1);
   gPad->SetLogy();
   auto *fr = new TH1D("fr", "elastic: data vs FRESCO (KD03, E_{lab}(p)=11.58 MeV);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr, shape]",
                       1, 15, 155);
   fr->SetMinimum(0.3);
   fr->SetMaximum(3e4);
   fr->Draw();
   frEl->Draw("L same");
   frEl155->Draw("L same");
   auto scaled = [](TH1D *h, double s, const char *nm) {
      auto *c = (TH1D *)h->Clone(nm);
      c->Scale(s);
      return c;
   };
   TH1D *elUs = scaled(elU, sU, "elUs"), *elGs = scaled(elG, sG, "elGs");
   elUs->Draw("E1 same");
   elGs->Draw("E1 same");
   auto *lg = new TLegend(0.48, 0.62, 0.89, 0.88);
   lg->AddEntry(frEl, "FRESCO elastic, 161 MeV beam", "l");
   lg->AddEntry(frEl155, "FRESCO elastic, 155 MeV beam", "l");
   lg->AddEntry(elUs, "UKF (acc-corrected, no #chi^{2})", "lp");
   lg->AddEntry(elGs, "GENFIT (acc-corrected, no #chi^{2})", "lp");
   lg->SetTextSize(0.032);
   lg->Draw();

   // ---- panel 2: data / FRESCO
   c1->cd(2);
   gPad->SetLogy();
   auto ratio = [&](TH1D *h, double s, const char *nm) {
      auto *r = (TH1D *)h->Clone(nm);
      r->Reset();
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double c = h->GetBinCenter(b), f = frEl->Eval(c);
         if (c < 15 || c > 155 || f <= 0 || h->GetBinContent(b) <= 0)
            continue;
         r->SetBinContent(b, h->GetBinContent(b) * s / f);
         r->SetBinError(b, h->GetBinError(b) * s / f);
      }
      return r;
   };
   TH1D *rU = ratio(elU, sU, "rU"), *rG = ratio(elG, sG, "rG");
   styleH(rU, kAzure + 2, 20);
   styleH(rG, kRed + 1, 21);
   rU->SetTitle("elastic: data / FRESCO;#theta_{cm} [deg];data / FRESCO");
   rU->GetXaxis()->SetRangeUser(15, 155);
   rU->SetMinimum(0.02);
   rU->SetMaximum(50);
   rU->Draw("E1");
   rG->Draw("E1 same");
   auto *one = new TLine(15, 1, 155, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();
   auto *lg2 = new TLegend(0.15, 0.72, 0.50, 0.88);
   lg2->AddEntry(rU, "UKF", "lp");
   lg2->AddEntry(rG, "GENFIT", "lp");
   lg2->Draw();

   // ---- panel 3: inelastic vs the three multipoles
   c1->cd(3);
   gPad->SetLogy();
   const double iU = fitScale(inU, frL1, 40, 100), iG = fitScale(inG, frL1, 40, 100);
   const double i3 = fitScale(inU, frL3, 40, 100), i2 = fitScale(inU, frL2, 40, 100);
   printf("\ninelastic normalisation over 40-100 deg: UKF/L1 x %.4g  GENFIT/L1 x %.4g\n", iU, iG);
   auto *fr3 = new TH1D("fr3", "E_{x} 5.5-7.0 group vs single-multipole DWBA;#theta_{cm} [deg];d#sigma/d#Omega [shape]",
                        1, 15, 155);
   fr3->SetMinimum(2e-3);
   fr3->SetMaximum(0.5);
   fr3->Draw();
   frL1->Draw("L same");
   auto *gL3 = (TGraph *)frL3->Clone("gL3");
   auto *gL2 = (TGraph *)frL2->Clone("gL2");
   for (int i = 0; i < gL3->GetN(); ++i) {
      gL3->SetPointY(i, gL3->GetPointY(i) * i3 / iU);
      gL2->SetPointY(i, gL2->GetPointY(i) * i2 / iU);
   }
   gL3->Draw("L same");
   gL2->Draw("L same");
   TH1D *inUs = scaled(inU, iU, "inUs"), *inGs = scaled(inG, iU, "inGs");
   inUs->Draw("E1 same");
   inGs->Draw("E1 same");
   auto *lg3 = new TLegend(0.45, 0.15, 0.89, 0.40);
   lg3->AddEntry(frL1, "FRESCO 6.094 (1^{-}, L=1)", "l");
   lg3->AddEntry(gL3, "FRESCO 6.728 (3^{-}, L=3)", "l");
   lg3->AddEntry(gL2, "FRESCO 7.012 (2^{+}, L=2)", "l");
   lg3->AddEntry(inUs, "UKF data", "lp");
   lg3->AddEntry(inGs, "GENFIT data", "lp");
   lg3->SetTextSize(0.030);
   lg3->Draw();

   // ---- panel 4: elastic linear, zoomed on the structure region
   c1->cd(4);
   auto *fr4 = new TH1D("fr4", "elastic, 55-155 deg (linear);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr, shape]", 1, 55,
                        155);
   double m = 0;
   for (int b = 1; b <= elUs->GetNbinsX(); ++b)
      if (elUs->GetBinCenter(b) > 55)
         m = std::max(m, elUs->GetBinContent(b));
   fr4->SetMinimum(0);
   fr4->SetMaximum(std::max(m, 25.0) * 1.5);
   fr4->Draw();
   frEl->Draw("L same");
   frEl155->Draw("L same");
   elUs->Draw("E1 same");
   elGs->Draw("E1 same");
   lg->Draw();

   TString png = here + "/plots/fresco_cmp_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
