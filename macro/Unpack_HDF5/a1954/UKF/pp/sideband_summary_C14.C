/// @file sideband_summary_C14.C
/// @brief UKF and GENFIT elastic distributions from the sideband method, together, against FRESCO.
///
/// The deliverable version of Fig. 8 (lower) of Ayyad et al., EPJ A 59:294 (2023): both fitters,
/// acceptance-corrected, extracted by locus-tracking + sideband subtraction, each normalised to
/// the DWBA over the same forward window. Run elastic_sideband_C14.C for "ukf" and "gf" first.
///
///   root -b -q 'sideband_summary_C14.C()'

void sideband_summary_C14(Double_t normLo = 20.0, Double_t normHi = 50.0,
                          TString frFile = "/home/yassid/a1954_C14_fresco/outputs/p14C_el_161_dsdo.dat")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   auto grab = [&](const char *f, const char *obj, const char *as) -> TH1D * {
      TFile *fi = TFile::Open(here + "/plots/" + f);
      if (!fi || fi->IsZombie()) {
         printf("\033[1;31mmissing %s -- run elastic_sideband_C14.C first\033[0m\n", f);
         return nullptr;
      }
      auto *h = (TH1D *)fi->Get(obj);
      if (!h)
         return nullptr;
      auto *c = (TH1D *)h->Clone(as);
      c->SetDirectory(nullptr);
      fi->Close();
      return c;
   };
   TH1D *sU = grab("elastic_sideband_gs_ukf.root", "dsigma_dOmega", "sU");
   TH1D *sG = grab("elastic_sideband_gs_gf.root", "dsigma_dOmega", "sG");
   TH1D *wU = grab("elastic_sideband_gs_ukf.root", "dsigma_dOmega_window", "wU");
   if (!sU || !sG || !wU)
      return;

   auto *fr = new TGraph();
   {
      std::ifstream in(frFile.Data());
      double th, xs;
      int n = 0;
      while (in >> th >> xs)
         fr->SetPoint(n++, th, xs);
   }
   auto norm = [&](TH1D *h) {
      double sn = 0, sd = 0;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double c = h->GetBinCenter(b), y = h->GetBinContent(b), e = h->GetBinError(b);
         if (c < normLo || c > normHi || y <= 0 || e <= 0)
            continue;
         double f = fr->Eval(c);
         sn += y * f / (e * e);
         sd += f * f / (e * e);
      }
      return sn > 0 ? sd / sn : 1.0;
   };
   double kU = norm(sU), kG = norm(sG), kW = norm(wU);
   sU->Scale(kU);
   sG->Scale(kG);
   wU->Scale(kW);

   printf("\n  theta_cm     FRESCO |   UKF sb   ratio |  GENFIT sb  ratio | UKF fixed-window  ratio\n");
   double cU = 0, cG = 0, cW = 0;
   int nU = 0, nG = 0, nW = 0;
   for (int b = 1; b <= sU->GetNbinsX(); ++b) {
      double c = sU->GetBinCenter(b);
      if (c < 18 || c > 148)
         continue;
      double f = fr->Eval(c), u = sU->GetBinContent(b), g = sG->GetBinContent(b), w = wU->GetBinContent(b);
      if (u <= 0 && g <= 0 && w <= 0)
         continue;
      printf("  %3.0f-%3.0f %10.4g | %8.4g %6.2f | %9.4g %6.2f | %14.4g %6.2f\n", sU->GetBinLowEdge(b),
             sU->GetBinLowEdge(b) + sU->GetBinWidth(b), f, u, f > 0 ? u / f : 0, g, f > 0 ? g / f : 0, w,
             f > 0 ? w / f : 0);
      if (f > 0 && u > 0) { cU += std::pow(std::log(u / f), 2); ++nU; }
      if (f > 0 && g > 0) { cG += std::pow(std::log(g / f), 2); ++nG; }
      if (f > 0 && w > 0) { cW += std::pow(std::log(w / f), 2); ++nW; }
   }
   printf("\nrms of ln(data/FRESCO):  UKF sideband %.3f (%d bins) | GENFIT sideband %.3f (%d bins) | "
          "UKF fixed window %.3f (%d bins)\n",
          nU ? std::sqrt(cU / nU) : 0, nU, nG ? std::sqrt(cG / nG) : 0, nG, nW ? std::sqrt(cW / nW) : 0, nW);

   TCanvas *c1 = new TCanvas("c1", "summary", 1300, 560);
   c1->Divide(2, 1);
   auto style = [](TH1D *h, int col, int m) {
      h->SetMarkerStyle(m);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
      h->SetMarkerSize(1.2);
   };
   style(sU, kAzure + 2, 20);
   style(sG, kRed + 1, 21);
   style(wU, kGray + 2, 25);

   c1->cd(1);
   gPad->SetLogy();
   auto *frm = new TH1D("frm", "^{14}C(p,p) elastic;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr, shape]", 1, 15, 150);
   frm->SetMinimum(0.5);
   frm->SetMaximum(3e3);
   frm->Draw();
   fr->SetLineWidth(3);
   fr->SetLineColor(kBlack);
   fr->Draw("L same");
   wU->Draw("E1 same");
   sU->Draw("E1 same");
   sG->Draw("E1 same");
   auto *lg = new TLegend(0.42, 0.66, 0.89, 0.88);
   lg->AddEntry(fr, "FRESCO DWBA (KD03)", "l");
   lg->AddEntry(sU, "UKF, sideband-subtracted", "lp");
   lg->AddEntry(sG, "GENFIT, sideband-subtracted", "lp");
   lg->AddEntry(wU, "UKF, fixed |E_{x}|<0.6 (old)", "lp");
   lg->SetTextSize(0.033);
   lg->Draw();

   c1->cd(2);
   gPad->SetLogy();
   auto ratio = [&](TH1D *h, const char *nm) {
      auto *r = (TH1D *)h->Clone(nm);
      r->Reset();
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double c = h->GetBinCenter(b), f = fr->Eval(c);
         if (c < 15 || c > 150 || f <= 0 || h->GetBinContent(b) <= 0)
            continue;
         r->SetBinContent(b, h->GetBinContent(b) / f);
         r->SetBinError(b, h->GetBinError(b) / f);
      }
      return r;
   };
   TH1D *rU = ratio(sU, "rU"), *rG = ratio(sG, "rG"), *rW = ratio(wU, "rW");
   style(rU, kAzure + 2, 20);
   style(rG, kRed + 1, 21);
   style(rW, kGray + 2, 25);
   rU->SetTitle("data / DWBA;#theta_{cm} [deg];ratio");
   rU->GetXaxis()->SetRangeUser(15, 150);
   rU->SetMinimum(0.05);
   rU->SetMaximum(20);
   rU->Draw("E1");
   rG->Draw("E1 same");
   rW->Draw("E1 same");
   auto *one = new TLine(15, 1, 150, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();
   lg->Draw();

   TString png = here + "/plots/sideband_summary_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
