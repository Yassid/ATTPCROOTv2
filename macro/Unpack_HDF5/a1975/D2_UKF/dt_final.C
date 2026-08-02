/// @file dt_final.C
/// @brief The 16C(d,t)15C spectrum under the selection that the resolution study settled on,
///        with the two controls that say what the number means.
///
/// The genfit resolution in this channel was never the problem. Measured on the isolated
/// 3.103 MeV level of 15C (a clean, unambiguous ruler, unlike the 0/0.740 doublet):
///
///     16C(p,d)15C, H2 target, the benchmark      sigma 0.252  FWHM 0.593
///     16C(d,d)16C elastic, SAME D2 data          sigma 0.418  FWHM 0.984
///     16C(d,t)15C, this selection                sigma 0.279  FWHM 0.657
///
/// (d,t) is BETTER than the elastic control taken on the same gas and only ~10 % wider than
/// the H2 (p,d) benchmark -- the "little bit worse from straggling" that was expected. With
/// FWHM 0.66 against a 0.740 MeV spacing the g.s. and the first excited state do separate,
/// and a two-gaussian fit at fixed 0.740 spacing is preferred over a single gaussian.
///
/// Three things had been hiding that, in order of size:
///
///  1. THE OLD KINEMATIC BOX. 45 < theta < 56 with 2 < KE_t < 20 MeV sits on the slow branch
///     right at the kinematic turnover (theta_max ~ 56 deg). That box is dominated by a
///     low-energy pile-up whose KE threshold maps onto a theta-dependent Ex edge: the
///     apparent "peak" walks 4.8 MeV across the box, which no beam energy and no KE scale
///     can flatten because it is a cut, not a state. Dropping the box and cutting KE > 5
///     instead takes S/B on the 3.103 line from 0.70 to 1.35.
///  2. NO BEAM GATE. The D2 beam is a cocktail; the component outside IC[900,1300]
///     reconstructs ~1.7 MeV high. Gating on the ion chamber (which the D2 chain never had
///     until the unpackIC_d2.C files existed) takes sigma from 0.425 to 0.336.
///  3. NO VERTEX FIDUCIAL. Restricting 50 < z < 700 mm takes sigma from 0.336 to 0.279.
///
///   root -b -q 'dt_final.C("/path/dt_kin_full.root")'

void dt_final(TString cache = "/tmp/dt_kin_full.root", TString plotOut = "plots/dt_final.png",
              double icLo = 900, double icHi = 1300, double keMin = 5, double vzLo = 50, double vzHi = 700,
              double chi2max = 5)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) {
      printf("no tree pk\n");
      return;
   }

   TString sel = Form("chi2ndf<%g&&ic>%g&&ic<%g&&ke>%g&&vertexz>%g&&vertexz<%g", chi2max, icLo, icHi, keMin, vzLo, vzHi);
   TString old = Form("chi2ndf<%g&&theta>45&&theta<56&&ke>2&&ke<20", chi2max);

   auto *cv = new TCanvas("cvf", "dt_final", 1600, 1000);
   cv->Divide(2, 2);

   // ---- 1. the spectrum, with both fits ------------------------------------------------
   cv->cd(1);
   t->Draw("ex>>hF(180,-3,6)", sel, "goff");
   auto *h = (TH1F *)gDirectory->Get("hF");
   h->SetTitle("^{15}C from 16C(d,t), recommended selection;E_{x}(^{15}C) [MeV];tritons");
   h->SetLineColor(kBlack);
   h->Draw("hist");

   TF1 *d = new TF1("d", "[0]*exp(-0.5*((x-[1])/[3])^2)+[2]*exp(-0.5*((x-[1]-0.740)/[3])^2)+[4]+[5]*x", -1.2, 2.2);
   d->SetParameters(40, 0, 40, 0.30, 10, 0);
   d->SetParLimits(1, -0.4, 0.4);
   d->SetParLimits(3, 0.10, 0.9);
   d->SetParLimits(0, 0, 1e5);
   d->SetParLimits(2, 0, 1e5);
   h->Fit(d, "QRN+");
   d->SetLineColor(kRed);
   d->SetNpx(400);
   d->DrawCopy("same");

   TF1 *g3 = new TF1("g3", "gaus(0)+pol1(3)", 1.9, 4.3);
   g3->SetParameters(30, 3.1, 0.30, 20, 0);
   g3->SetParLimits(1, 2.7, 3.5);
   g3->SetParLimits(2, 0.10, 1.0);
   h->Fit(g3, "QRN+");
   g3->SetLineColor(kBlue + 1);
   g3->SetNpx(400);
   g3->DrawCopy("same");

   for (double lv : {0.0, 0.740, 3.103, 4.220}) {
      auto *l = new TLine(lv, 0, lv, 0.95 * h->GetMaximum());
      l->SetLineColor(kGray + 2);
      l->SetLineStyle(2);
      l->Draw();
   }
   double sgD = std::fabs(d->GetParameter(3)), sg3 = std::fabs(g3->GetParameter(2));
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.036);
   tx->DrawLatex(0.50, 0.85, Form("g.s. %+.3f, #sigma %.3f", d->GetParameter(1), sgD));
   tx->DrawLatex(0.50, 0.80, Form("3.103: %.3f, #sigma %.3f", g3->GetParameter(1), sg3));
   tx->DrawLatex(0.50, 0.75, Form("FWHM %.3f MeV", 2.355 * sg3));

   printf("\n=== dt_final ===\nselection: %s\nN = %.0f\n", sel.Data(), h->GetEntries());
   printf("0/0.740 doublet : g.s. %+.3f  sigma %.3f  FWHM %.3f   A(gs) %.0f  A(0.740) %.0f  chi2/ndf %.2f\n",
          d->GetParameter(1), sgD, 2.355 * sgD, d->GetParameter(0), d->GetParameter(2),
          d->GetNDF() > 0 ? d->GetChisquare() / d->GetNDF() : -1);
   printf("3.103 level     : mu   %+.3f  sigma %.3f  FWHM %.3f   chi2/ndf %.2f\n", g3->GetParameter(1), sg3,
          2.355 * sg3, g3->GetNDF() > 0 ? g3->GetChisquare() / g3->GetNDF() : -1);
   printf("controls        : (p,d) H2 sigma 0.252 | (d,d) elastic D2 sigma 0.418\n");

   // ---- 2. what the old box was actually selecting -------------------------------------
   cv->cd(2);
   t->Draw("ex>>hO(180,-3,6)", old, "goff");
   auto *ho = (TH1F *)gDirectory->Get("hO");
   ho->SetTitle("old box 45<#theta<56, 2<KE<20 (no beam gate);E_{x}(^{15}C) [MeV];tritons");
   ho->SetLineColor(kGray + 2);
   ho->Draw("hist");
   auto *hs = (TH1F *)h->Clone("hs");
   hs->Scale(ho->Integral() / std::max(1.0, h->Integral()));
   hs->SetLineColor(kRed);
   hs->Draw("hist same");
   auto *lg = new TLegend(0.45, 0.75, 0.88, 0.88);
   lg->AddEntry(ho, "old box", "l");
   lg->AddEntry(hs, "recommended (scaled)", "l");
   lg->Draw();

   // ---- 3. why the old box fails: the threshold edge walks with theta -------------------
   cv->cd(3);
   gPad->SetLogz();
   t->Draw("ex:theta>>hE(70,10,80,110,-3,8)", Form("chi2ndf<%g&&ic>%g&&ic<%g&&ke>2&&ke<20", chi2max, icLo, icHi),
           "colz");
   auto *he = (TH2F *)gDirectory->Get("hE");
   he->SetTitle("with the 2<KE<20 box: the edge, not a state;#theta_{lab} [deg];E_{x} [MeV]");
   for (double lv : {0.0, 0.740, 3.103}) {
      auto *l = new TLine(10, lv, 80, lv);
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
   }
   {
      auto *b = new TBox(45, -3, 56, 8);
      b->SetFillStyle(0);
      b->SetLineColor(kGreen + 1);
      b->SetLineWidth(3);
      b->Draw();
   }

   // ---- 4. the beam gate ---------------------------------------------------------------
   cv->cd(4);
   TString noIC = Form("chi2ndf<%g&&ke>%g&&vertexz>%g&&vertexz<%g", chi2max, keMin, vzLo, vzHi);
   t->Draw("ex>>hIn(180,-3,6)", noIC + Form("&&ic>%g&&ic<%g", icLo, icHi), "goff");
   t->Draw("ex>>hOut(180,-3,6)", noIC + Form("&&ic>=0&&(ic<=%g||ic>=%g)", icLo, icHi), "goff");
   auto *hi = (TH1F *)gDirectory->Get("hIn");
   auto *hou = (TH1F *)gDirectory->Get("hOut");
   hi->SetTitle("beam gate: in (black) vs out (red);E_{x}(^{15}C) [MeV];tritons");
   hi->SetLineColor(kBlack);
   hou->SetLineColor(kRed);
   hi->Draw("hist");
   hou->Draw("hist same");
   for (double lv : {0.0, 0.740, 3.103}) {
      auto *l = new TLine(lv, 0, lv, 0.95 * hi->GetMaximum());
      l->SetLineColor(kGray + 2);
      l->SetLineStyle(2);
      l->Draw();
   }

   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
