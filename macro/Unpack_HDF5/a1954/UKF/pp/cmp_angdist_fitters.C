/// @file cmp_angdist_fitters.C
/// @brief Overlay the UKF and GENFIT acceptance-corrected 14C(p,p') angular distributions.
///
/// Both fitters measure the same physics, so their acceptance-CORRECTED shapes must agree even
/// though their raw yields do not. Where they disagree, one of the two acceptances is wrong.
/// Run apply_acceptance_C14.C with outTag "_ukf" and "_gf" first; this reads the four
/// plots/angdist_{gs,ex1}_{ukf,gf}.root it writes.
///
/// Panels: (1) the two acceptances, (2) the two raw yields, (3) the corrected dsigma/dOmega
/// shapes area-normalised over the common range, (4) the GENFIT/UKF ratio of the corrected
/// shapes -- flat at 1 means the correction closed the gap.
///
/// `variant` selects which pair to compare: "" reads angdist_<level>_{ukf,gf}.root, "_nc" reads
/// the no-chi2-cut pair angdist_<level>_{ukf,gf}_nc.root.
///
///   root -b -q 'cmp_angdist_fitters.C("gs")'
///   root -b -q 'cmp_angdist_fitters.C("gs",20,150,"_nc")'

void cmp_angdist_fitters(TString level = "gs", Double_t cmMin = 20.0, Double_t cmMax = 150.0, TString variant = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fu = TFile::Open(here + "/plots/angdist_" + level + "_ukf" + variant + ".root");
   TFile *fg = TFile::Open(here + "/plots/angdist_" + level + "_gf" + variant + ".root");
   if (!fu || fu->IsZombie() || !fg || fg->IsZombie()) {
      printf("\033[1;31mrun apply_acceptance_C14.C with outTag _ukf and _gf first\033[0m\n");
      return;
   }
   auto get = [](TFile *f, const char *n, const char *as) {
      auto *h = (TH1D *)f->Get(n);
      return h ? (TH1D *)h->Clone(as) : nullptr;
   };
   TH1D *accU = get(fu, "acceptance", "accU"), *accG = get(fg, "acceptance", "accG");
   TH1D *rawU = get(fu, "raw", "rawU"), *rawG = get(fg, "raw", "rawG");
   TH1D *dsU = get(fu, "dsigma_dOmega", "dsU"), *dsG = get(fg, "dsigma_dOmega", "dsG");
   if (!accU || !accG || !dsU || !dsG) {
      printf("\033[1;31mmissing histograms in one of the files\033[0m\n");
      return;
   }

   // area-normalise the corrected shapes over the bins BOTH fitters populate, so the comparison
   // is of shape only -- neither carries a luminosity
   double sU = 0, sG = 0;
   for (int b = 1; b <= dsU->GetNbinsX(); ++b) {
      double c = dsU->GetBinCenter(b);
      if (c < cmMin || c > cmMax || dsU->GetBinContent(b) <= 0 || dsG->GetBinContent(b) <= 0)
         continue;
      sU += dsU->GetBinContent(b);
      sG += dsG->GetBinContent(b);
   }
   auto *nU = (TH1D *)dsU->Clone("nU");
   auto *nG = (TH1D *)dsG->Clone("nG");
   if (sU > 0) nU->Scale(1.0 / sU);
   if (sG > 0) nG->Scale(1.0 / sG);

   auto *ratio = (TH1D *)nG->Clone("ratio");
   ratio->Reset();
   printf("\n===== level %s%s : GENFIT vs UKF, acceptance-corrected =====\n", level.Data(), variant.Data());
   printf("  theta_cm    accU    accG   accG/accU     rawU    rawG   rawG/rawU    shapeG/shapeU\n");
   for (int b = 1; b <= nU->GetNbinsX(); ++b) {
      double c = nU->GetBinCenter(b);
      if (c < cmMin || c > cmMax)
         continue;
      double u = nU->GetBinContent(b), g = nG->GetBinContent(b);
      double au = accU->GetBinContent(b), ag = accG->GetBinContent(b);
      double ru = rawU->GetBinContent(b), rg = rawG->GetBinContent(b);
      if (u > 0 && g > 0) {
         double r = g / u;
         double eu = nU->GetBinError(b) / u, eg = nG->GetBinError(b) / g;
         ratio->SetBinContent(b, r);
         ratio->SetBinError(b, r * std::sqrt(eu * eu + eg * eg));
      }
      printf("  %3.0f-%3.0f  %6.3f  %6.3f  %8.3f  %7.0f %7.0f  %8.3f  %13s\n", nU->GetBinLowEdge(b),
             nU->GetBinLowEdge(b) + nU->GetBinWidth(b), au, ag, au > 0 ? ag / au : 0, ru, rg,
             ru > 0 ? rg / ru : 0,
             (u > 0 && g > 0) ? TString::Format("%.3f", g / u).Data() : "-");
   }

   TCanvas *c = new TCanvas("c", "fitter comparison", 1500, 950);
   c->Divide(2, 2);
   auto style = [](TH1D *h, int col, int mk) {
      h->SetMarkerStyle(mk);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
   };
   style(accU, kAzure + 2, 20);
   style(accG, kRed + 1, 21);
   style(rawU, kAzure + 2, 20);
   style(rawG, kRed + 1, 21);
   style(nU, kAzure + 2, 20);
   style(nG, kRed + 1, 21);
   style(ratio, kBlack, 20);

   c->cd(1);
   accU->SetTitle(TString::Format("simulated acceptance (%s);#theta_{cm} [deg];acceptance", level.Data()));
   accU->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   accU->GetYaxis()->SetRangeUser(0, 1.15);
   accU->Draw("E1");
   accG->Draw("E1 same");
   auto *l1 = new TLegend(0.55, 0.18, 0.88, 0.35);
   l1->AddEntry(accU, "UKF", "lp");
   l1->AddEntry(accG, "GENFIT", "lp");
   l1->Draw();

   c->cd(2);
   gPad->SetLogy();
   rawU->SetTitle("measured yield;#theta_{cm} [deg];counts");
   rawU->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   rawU->Draw("E1");
   rawG->Draw("E1 same");
   l1->Draw();

   c->cd(3);
   gPad->SetLogy();
   nU->SetTitle("corrected d#sigma/d#Omega, area-normalised;#theta_{cm} [deg];shape [a.u.]");
   nU->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   nU->Draw("E1");
   nG->Draw("E1 same");
   l1->Draw();

   c->cd(4);
   ratio->SetTitle("GENFIT / UKF of the corrected shape;#theta_{cm} [deg];ratio");
   ratio->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   ratio->GetYaxis()->SetRangeUser(0, 3);
   ratio->Draw("E1");
   auto *one = new TLine(cmMin - 5, 1, cmMax + 5, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();

   TString png = here + "/plots/cmp_angdist_" + level + variant + ".png";
   c->SaveAs(png);
   printf("wrote %s\n\n", png.Data());
}
