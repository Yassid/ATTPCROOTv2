/// @file dt_summary.C
/// @brief Summary view of the a1975 16C(d,t)15C channel from the ex_dt_a1975.C cache.
///
/// Three things this makes visible, in order of how much they cost the spectrum:
///
/// 1. ACCEPTANCE, not calibration, is what wrecks the integrated spectrum. The (d,t) g.s.
///    locus rises very steeply forward (theta 27 deg -> KE_t ~62 MeV), so the 2-20 MeV KE
///    box the old analysis used only touches the locus above theta ~45 deg. Everything the
///    same box admits below 45 deg is far off-locus and reconstructs at Ex 8-20 MeV -- that
///    is the broad "background" bump, and it outnumbers the real peak ~4:1.
/// 2. The BEAM GATE. The D2 runs do have an ion chamber (unpackIC_d2.C); the beam is a
///    cocktail with a second component at IC ~2060 in comparable numbers to the 16C at
///    ~1150. Nothing in the D2 chain has ever gated on it.
/// 3. What is left after both: the 15C peaks and their width.
///
///   root -b -q 'dt_summary.C("/path/dt_kin_full.root")'
///
/// Cache branches used: ke, theta, ex, chi2ndf, vertexz, ic (ic<0 = no IC file for that run).

void dt_summary(TString cache = "triton_kin_dt.root", TString plotOut = "plots/dt_summary.png",
                double icLo = 900, double icHi = 1300, double thLo = 45, double thHi = 56, double keLo = 2,
                double keHi = 20, double chi2Cut = 5)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree pk in %s\n", cache.Data()); return; }
   printf("cache: %lld candidates\n", t->GetEntries());

   TString onLocus = Form("theta>%g&&theta<%g&&ke>%g&&ke<%g&&chi2ndf<%g", thLo, thHi, keLo, keHi, chi2Cut);
   TString offLocus = Form("theta>25&&theta<%g&&ke>%g&&ke<%g&&chi2ndf<%g", thLo, keLo, keHi, chi2Cut);
   TString icIn = Form("ic>%g&&ic<%g", icLo, icHi);
   TString icOut = Form("ic>=0&&(ic<=%g||ic>=%g)", icLo, icHi);

   auto *cv = new TCanvas("cv", "dt_summary", 1700, 1050);
   cv->Divide(3, 2);

   cv->cd(1); gPad->SetLogz();
   auto *hk = new TH2F("hk", "KE vs #theta_{lab}, all candidates;#theta_{lab} [deg];KE_{t} [MeV]", 90, 0, 90, 120, 0, 90);
   t->Draw("ke:theta>>hk", Form("chi2ndf<%g", chi2Cut), "colz");
   auto *bx = new TBox(thLo, keLo, thHi, keHi);
   bx->SetFillStyle(0); bx->SetLineColor(kGreen + 1); bx->SetLineWidth(3); bx->Draw();

   cv->cd(2);
   auto *hIC = new TH1F("hIC", "IC of events with a triton candidate;IC max ADC;", 200, 0, 2500);
   t->Draw("ic>>hIC", "ic>=0", "hist");
   auto *l1 = new TLine(icLo, 0, icLo, hIC->GetMaximum()); l1->SetLineColor(kRed); l1->Draw();
   auto *l2 = new TLine(icHi, 0, icHi, hIC->GetMaximum()); l2->SetLineColor(kRed); l2->Draw();

   cv->cd(3);
   auto *hOff = new TH1F("hOff", "OFF-locus 25<#theta<45, same KE box;E_{x}(^{15}C) [MeV];", 90, -6, 25);
   t->Draw("ex>>hOff", offLocus, "hist");

   cv->cd(4);
   auto *hOn = new TH1F("hOn", "ON-locus (green box);E_{x}(^{15}C) [MeV];tritons", 90, -6, 20);
   t->Draw("ex>>hOn", onLocus, "hist");

   cv->cd(5);
   auto *hOnIn = new TH1F("hOnIn", "on-locus: IC in gate (black) / outside (red);E_{x} [MeV];", 90, -6, 20);
   t->Draw("ex>>hOnIn", onLocus + "&&" + icIn, "hist");
   auto *hOnOut = new TH1F("hOnOut", "", 90, -6, 20);
   hOnOut->SetLineColor(kRed);
   t->Draw("ex>>hOnOut", onLocus + "&&" + icOut, "hist same");

   cv->cd(6);
   auto *hFin = (TH1F *)hOnIn->Clone("hFin");
   hFin->SetTitle("FINAL: on-locus + beam gate;E_{x}(^{15}C) [MeV];tritons");
   hFin->Draw("hist");
   for (double lv : {0.0, 0.740, 3.103, 4.220, 4.657}) {
      auto *l = new TLine(lv, 0, lv, hFin->GetMaximum());
      l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw();
   }
   int mb = hFin->GetMaximumBin();
   double x0 = hFin->GetBinCenter(mb);
   TF1 g("g", "gaus", x0 - 1.5, x0 + 1.5);
   hFin->Fit(&g, "RQ0");
   printf("\non-locus  : %.0f   (IC in gate %.0f, outside %.0f, no IC %.0f)\n", hOn->GetEntries(),
          hOnIn->GetEntries(), hOnOut->GetEntries(),
          hOn->GetEntries() - hOnIn->GetEntries() - hOnOut->GetEntries());
   printf("off-locus : %.0f  <- the Ex 8-20 background\n", hOff->GetEntries());
   printf("final peak: centroid %.2f MeV, sigma %.2f, FWHM %.2f\n", g.GetParameter(1), fabs(g.GetParameter(2)),
          2.355 * fabs(g.GetParameter(2)));
   g.SetLineColor(kRed);
   g.DrawCopy("same");
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
