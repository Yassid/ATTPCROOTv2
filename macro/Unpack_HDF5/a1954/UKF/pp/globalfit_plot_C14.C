/// @file globalfit_plot_C14.C
/// @brief Figure of the explorer's global fit (same model and fit as globalfit_check_C14.C), one column per cache.
///
/// Each column: Ex spectrum with the total model and every component; pulls (n - m)/sqrt(m) below, or deviance
/// residuals when there is no background (the model vanishes between peaks and Pearson diverges).
/// withGS adds the ground state with its OWN free centroid and free sigma (the g.s. centroid walks with angle, and
/// the level width law extrapolated to E = 0 gives 0.057 MeV, below the measured width). gsGauss = 2 adds a second,
/// broad g.s. gaussian (free amplitude, centroid, sigma) for the resolution spread. gsCmLo/gsCmHi take the events
/// below gsExCut from their own theta_cm window (e.g. 0-60, where the g.s. is narrow) and those above from cmLo-cmHi.
/// Curves are snapshotted into TGraphs right after each fit (a TF1 on a C++ function re-reads the
/// namespace globals lazily at draw time).
///
///   root -l -b -q globalfit_plot_C14.C                          # s013 vs tc, reference fit -> plots/globalfit_C14_s013_vs_tc.png
///   root -l -b -q 'globalfit_plot_C14.C("tc", true)'            # ADOPTED g.s. fit -> plots/globalfit_C14_tc_gs2_gscm0-60_nobg.png
/// Earlier variants (args: exLo, exHi, useBg 0/1, gsGauss, gsCmLo, gsCmHi; gsCmLo -1 = g.s. from cmLo-cmHi):
///   'globalfit_plot_C14.C("tc", true, -3.0, 10.4, 1, 1, -1, -1)'   # one g.s. gaussian + linear bg -> _gs.png
///   'globalfit_plot_C14.C("tc", true, -3.0, 10.4, 0, 1, -1, -1)'   # one g.s. gaussian, no bg     -> _gs_nobg.png
///   'globalfit_plot_C14.C("tc", true, -3.0, 10.4, 0, 2, -1, -1)'   # two gaussians, g.s. 20-140  -> _gs2_nobg.png
namespace gp {
int NL = 5;
bool gs = false, gs2 = false;
double E[6] = {6.091, 6.728, 7.012, 7.341, 8.317}, SG[6];
double shift = -0.011, muN = 9.178, sgN = 0.296;
TH1D *hps = nullptr;
// parameters: 0..NL-1 level amplitudes (0 = g.s. core when gs), NL 14N, NL+1 PS, NL+2/NL+3 background,
// NL+4 g.s. sigma, NL+5 g.s. mu, NL+6 broad g.s. amplitude, NL+7 broad mu, NL+8 broad sigma
double gauss(double x, double mu, double sg) { return std::exp(-0.5 * std::pow((x - mu) / sg, 2)); }
double model(double *x, double *p)
{
   double s = p[NL + 2] + p[NL + 3] * x[0];
   for (int i = 0; i < NL; ++i)
      s += (gs && i == 0) ? p[0] * gauss(x[0], p[NL + 5], p[NL + 4]) : p[i] * gauss(x[0], E[i] + shift, SG[i]);
   if (gs2) s += p[NL + 6] * gauss(x[0], p[NL + 7], p[NL + 8]);
   s += p[NL] * gauss(x[0], muN, sgN);
   if (hps) s += p[NL + 1] * hps->Interpolate(x[0]);
   return s;
}
} // namespace gp

void globalfit_plot_C14(TString which = "s013,tc", Bool_t withGS = false, Double_t exLoIn = -999, Double_t exHi = 10.4,
                        Int_t useBgIn = -1, Int_t gsGauss = 2, Double_t gsCmLo = 0, Double_t gsCmHi = 60,
                        Double_t gsExCut = 3.0, Double_t cmLo = 20, Double_t cmHi = 140, Double_t binW = 0.05,
                        Double_t vzLo = 10, Double_t vzHi = 490, Double_t chi2Cut = 5.0, Double_t sig0 = 0.132,
                        Double_t dSig = 0.0123)
{
   // defaults: without the g.s., the reference fit (5.0-10.4 MeV, linear background); with it, the ADOPTED g.s.
   // configuration (2026-09-13): -3.0-10.4 MeV, no background, two g.s. gaussians, g.s. events from theta_cm 0-60
   const double exLo = exLoIn < -900 ? (withGS ? -3.0 : 5.0) : exLoIn;
   const bool useBg = useBgIn < 0 ? !withGS : useBgIn != 0;
   const bool splitCm = withGS && gsCmLo >= 0 && gsCmHi > gsCmLo;
   gp::gs = withGS;
   gp::gs2 = withGS && gsGauss == 2;
   gStyle->SetOptStat(0);
   gStyle->SetTitleFontSize(gp::gs2 ? 0.038 : 0.05);
   gStyle->SetTitleAlign(13);
   gStyle->SetTitleX(0.12);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fp = TFile::Open(here + "/plots/phasespace_1n_C14.root");
   gp::hps = (TH1D *)fp->Get("hPS");
   gp::hps->SetDirectory(nullptr);
   gp::hps->Scale(1.0 / gp::hps->Integral());

   const double lit[5] = {6.091, 6.728, 7.012, 7.341, 8.317};
   gp::NL = withGS ? 6 : 5;
   for (int i = 0; i < gp::NL; ++i) gp::E[i] = withGS ? (i == 0 ? 0.0 : lit[i - 1]) : lit[i];
   for (int i = 0; i < gp::NL; ++i) gp::SG[i] = sig0 + dSig * (gp::E[i] - 6.094);

   // reference palette, light mode: inelastic levels slots 1-5, 14N slot 6, phase space slot 7, g.s. slot 8
   const char *lvHex[5] = {"#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"};
   const char *hexN = "#008300", *hexPS = "#4a3aa7", *hexGS = "#e34948";
   std::map<TString, TString> label = {{"s013", "s013 (0.130 deg/MeV, page default)"}, {"tc", "tc (0.056 deg/MeV)"}};
   std::unique_ptr<TObjArray> list(which.Tokenize(","));
   const int ncol = list->GetEntries();
   const int NL = gp::NL, NPAR = NL + 4 + (withGS ? 2 : 0) + (gp::gs2 ? 3 : 0);
   const int nb = (int)std::lround((exHi - exLo) / binW), NP = 1500;
   const double sq = std::sqrt(2 * TMath::Pi());
   const int nFit = NPAR - (useBg ? 0 : 2); // free parameters, for ndf

   auto *c = new TCanvas("cgp", "", ncol == 1 ? 1100 : 1800, 1000);
   c->SetFillColor(TColor::GetColor("#fcfcfb"));
   for (int k = 0; k < ncol; ++k) {
      TString tag = ((TObjString *)list->At(k))->GetString();
      TFile *fd = TFile::Open(here + "/plots/proton_kin_cat5_" + tag + ".root");
      auto *t = (TNtuple *)fd->Get("pk");
      auto *h = new TH1D(Form("h%d", k), "", nb, exLo, exHi);
      TString common = Form("chi2ndf<%g && vertexz>=%g && vertexz<=%g", chi2Cut, vzLo, vzHi);
      TString sel = splitCm ? Form("%s && ((ex<%g && thcm>=%g && thcm<=%g) || (ex>=%g && thcm>=%g && thcm<=%g))",
                                   common.Data(), gsExCut, gsCmLo, gsCmHi, gsExCut, cmLo, cmHi)
                            : Form("%s && thcm>=%g && thcm<=%g", common.Data(), cmLo, cmHi);
      t->Draw(Form("ex>>h%d", k), sel, "goff");

      auto *F = new TF1(Form("F%d", k), gp::model, exLo, exHi, NPAR);
      for (int i = 0; i < NL; ++i) {
         F->SetParameter(i, std::max(1.0, h->GetBinContent(h->FindBin(gp::E[i] + gp::shift)) * 0.5));
         F->SetParLimits(i, 0, 1e5);
      }
      F->SetParameter(NL, 20);                       F->SetParLimits(NL, 0, 1e5);
      F->SetParameter(NL + 1, h->Integral() * 0.05); F->SetParLimits(NL + 1, 0, 1e7);
      F->SetParameter(NL + 2, 2);                    F->SetParameter(NL + 3, 0);
      if (!useBg) { F->FixParameter(NL + 2, 0); F->FixParameter(NL + 3, 0); }
      if (withGS) {
         F->SetParameter(NL + 4, 0.15); F->SetParLimits(NL + 4, 0.03, 1.0);
         F->SetParameter(NL + 5, 0.0);  F->SetParLimits(NL + 5, -1.0, 1.0);
      }
      if (gp::gs2) {
         F->SetParameter(NL + 4, 0.14);                    F->SetParLimits(NL + 4, 0.03, 0.5);
         F->SetParameter(NL + 6, 0.05 * h->GetMaximum());  F->SetParLimits(NL + 6, 0, 1e5);
         F->SetParameter(NL + 7, -0.5);                    F->SetParLimits(NL + 7, -2.0, 1.0);
         F->SetParameter(NL + 8, 0.7);                     F->SetParLimits(NL + 8, 0.2, 10.0);
      }
      TFitResultPtr r = h->Fit(F, "RQNSL");
      for (int pass = 0; pass < 3 && ((int)r->Status() != 0 || r->CovMatrixStatus() < 2); ++pass) r = h->Fit(F, "RQNSL");

      // snapshot: total, each level, broad g.s., 14N, PS, background
      double p[16];
      for (int j = 0; j < NPAR; ++j) p[j] = F->GetParameter(j);
      auto sigOf = [&](int i) { return (withGS && i == 0) ? p[NL + 4] : gp::SG[i]; };
      auto muOf = [&](int i) { return (withGS && i == 0) ? p[NL + 5] : gp::E[i] + gp::shift; };
      auto *gTot = new TGraph(NP), *gBg = new TGraph(NP), *gN = new TGraph(NP), *gPS = new TGraph(NP), *gGS2 = new TGraph(NP);
      std::vector<TGraph *> gL(NL);
      for (auto &g : gL) g = new TGraph(NP);
      for (int ip = 0; ip < NP; ++ip) {
         double x = exLo + (exHi - exLo) * ip / (NP - 1), bg = p[NL + 2] + p[NL + 3] * x;
         gTot->SetPoint(ip, x, gp::model(&x, p));
         gBg->SetPoint(ip, x, bg);
         for (int i = 0; i < NL; ++i) gL[i]->SetPoint(ip, x, p[i] * gp::gauss(x, muOf(i), sigOf(i)));
         gGS2->SetPoint(ip, x, gp::gs2 ? p[NL + 6] * gp::gauss(x, p[NL + 7], p[NL + 8]) : 0.0);
         gN->SetPoint(ip, x, p[NL] * gp::gauss(x, gp::muN, gp::sgN));
         gPS->SetPoint(ip, x, p[NL + 1] * gp::hps->Interpolate(x));
      }
      // Pearson (what the reference prints) and Baker-Cousins chi2_lambda (the likelihood statistic, finite as m -> 0)
      double pear = 0, bc = 0; int nfree = 0, nTiny = 0;
      auto *hPull = new TH1D(Form("pull%d", k), "", nb, exLo, exHi);
      for (int b = 1; b <= nb; ++b) {
         double x = h->GetBinCenter(b), m = std::max(gp::model(&x, p), 1e-12), n = h->GetBinContent(b);
         double dev = 2 * (m - n + (n > 0 ? n * std::log(n / m) : 0));
         pear += (n - m) * (n - m) / m; bc += dev; ++nfree;
         if (m < 0.1 && n > 0) ++nTiny;
         hPull->SetBinContent(b, useBg ? (n - m) / std::sqrt(m) : (n >= m ? 1 : -1) * std::sqrt(std::max(dev, 0.0)));
      }
      double psInt = 0;
      for (int b = 1; b <= nb; ++b) { double x = h->GetBinCenter(b); psInt += p[NL + 1] * gp::hps->Interpolate(x); }
      // areas A*sigma*sqrt(2pi)/binW; where sigma floats too, propagate the covariance
      TMatrixDSym V = r->GetCovarianceMatrix();
      auto areaAS = [&](int iA, int iS) { return p[iA] * p[iS] * sq / binW; };
      auto errAS = [&](int iA, int iS) {
         double A = p[iA], s = p[iS];
         return std::sqrt(s * s * V(iA, iA) + A * A * V(iS, iS) + 2 * A * s * V(iA, iS)) * sq / binW;
      };
      auto area = [&](int i) { return (withGS && i == 0) ? areaAS(0, NL + 4) : p[i] * gp::SG[i] * sq / binW; };
      auto areaErr = [&](int i) { return (withGS && i == 0) ? errAS(0, NL + 4) : F->GetParError(i) * gp::SG[i] * sq / binW; };
      double gsTot = 0, gsTotErr = 0;
      if (withGS) {
         std::vector<std::pair<int, double>> grad = {{0, p[NL + 4]}, {NL + 4, p[0]}};
         if (gp::gs2) { grad.push_back({NL + 6, p[NL + 8]}); grad.push_back({NL + 8, p[NL + 6]}); }
         double var = 0;
         for (auto &a : grad) for (auto &b : grad) var += a.second * b.second * V(a.first, b.first);
         gsTot = area(0) + (gp::gs2 ? areaAS(NL + 6, NL + 8) : 0);
         gsTotErr = std::sqrt(var) * sq / binW;
      }

      // pads
      c->cd();
      double x0 = double(k) / ncol, x1 = double(k + 1) / ncol;
      auto *pt = new TPad(Form("pt%d", k), "", x0, 0.30, x1, 1.0);
      auto *pb = new TPad(Form("pb%d", k), "", x0, 0.0, x1, 0.30);
      for (auto *pd : {pt, pb}) { pd->SetFillColor(TColor::GetColor("#fcfcfb")); pd->SetLeftMargin(0.13); pd->SetRightMargin(0.03); pd->Draw(); }
      pt->SetBottomMargin(0.02); pt->SetTopMargin(0.09); pb->SetTopMargin(0.03); pb->SetBottomMargin(0.30);

      pt->cd();
      pt->SetGridy(); gStyle->SetGridColor(TColor::GetColor("#e1e0d9"));
      if (withGS) pt->SetLogy();
      TString cmNote = splitCm ? Form("#theta_{cm} %g-%g#circ (g.s. %g-%g#circ)", cmLo, cmHi, gsCmLo, gsCmHi)
                               : Form("#theta_{cm} %g-%g#circ", cmLo, cmHi);
      TString gsNote = "";
      if (gp::gs2)
         gsNote = Form("   g.s. core #mu %.3f #sigma %.3f, broad #mu %.3f #sigma %.3f%s", p[NL + 5], p[NL + 4], p[NL + 7],
                       p[NL + 8], useBg ? "" : ", no bg");
      else if (withGS)
         gsNote = Form("   g.s. #mu %.3f  #sigma %.3f MeV%s", p[NL + 5], p[NL + 4], useBg ? "" : "   no background");
      h->SetTitle(Form("%s   %s%s;;counts / %.0f keV", label[tag].Data(), cmNote.Data(), gsNote.Data(), binW * 1000));
      h->SetMarkerStyle(20); h->SetMarkerSize(withGS ? 0.5 : 0.8); h->SetLineColor(TColor::GetColor("#0b0b0b"));
      if (withGS) { h->SetMinimum(0.5); h->SetMaximum(h->GetMaximum() * 20); }
      else { h->SetMaximum(h->GetMaximum() * 1.25); h->SetMinimum(0); }
      h->GetXaxis()->SetLabelSize(0);
      h->GetYaxis()->SetTitleSize(0.05); h->GetYaxis()->SetLabelSize(0.045); h->GetYaxis()->SetTitleOffset(1.0);
      h->Draw("E1");
      // components below 0.3 counts are not drawn on the log scale
      auto clip = [&](TGraph *g) { if (withGS) for (int ip = 0; ip < g->GetN(); ++ip) if (g->GetY()[ip] < 0.3) g->SetPoint(ip, g->GetX()[ip], 0.3); };
      if (useBg) { gBg->SetLineColor(TColor::GetColor("#898781")); gBg->SetLineStyle(2); gBg->SetLineWidth(2); clip(gBg); gBg->Draw("L"); }
      gPS->SetLineColor(TColor::GetColor(hexPS)); gPS->SetLineWidth(2); gPS->SetLineStyle(7); clip(gPS); gPS->Draw("L");
      gN->SetLineColor(TColor::GetColor(hexN)); gN->SetLineWidth(2); clip(gN); gN->Draw("L");
      if (gp::gs2) { gGS2->SetLineColor(TColor::GetColor(hexGS)); gGS2->SetLineWidth(2); gGS2->SetLineStyle(2); clip(gGS2); gGS2->Draw("L"); }
      for (int i = 0; i < NL; ++i) {
         gL[i]->SetLineColor(TColor::GetColor((withGS && i == 0) ? hexGS : lvHex[i - (withGS ? 1 : 0)]));
         gL[i]->SetLineWidth(2); clip(gL[i]); gL[i]->Draw("L");
      }
      gTot->SetLineColor(TColor::GetColor("#0b0b0b")); gTot->SetLineWidth(withGS ? 2 : 3); clip(gTot); gTot->Draw("L");
      h->Draw("E1 same");

      auto *leg = withGS ? new TLegend(0.30, 0.55, 0.97, 0.89) : new TLegend(0.56, 0.40, 0.96, 0.89);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(withGS ? 0.030 : 0.036);
      if (withGS) leg->SetNColumns(2);
      leg->AddEntry(h, Form("data  N = %.0f", h->Integral()), "pe");
      leg->AddEntry(gTot, useBg ? Form("total   #chi^{2}/ndf %.0f/%d = %.2f", pear, nfree - nFit, pear / (nfree - nFit))
                                : Form("total   #chi^{2}_{#lambda}/ndf %.0f/%d = %.2f", bc, nfree - nFit, bc / (nfree - nFit)), "l");
      for (int i = 0; i < NL; ++i) {
         if (withGS && i == 0) {
            if (gp::gs2) {
               leg->AddEntry(gL[i], Form("g.s. core   %.0f #pm %.0f", area(i), areaErr(i)), "l");
               leg->AddEntry(gGS2, Form("g.s. broad   %.0f #pm %.0f", areaAS(NL + 6, NL + 8), errAS(NL + 6, NL + 8)), "l");
            } else
               leg->AddEntry(gL[i], Form("g.s.   %.0f #pm %.0f", area(i), areaErr(i)), "l");
         } else
            leg->AddEntry(gL[i], Form("%.3f MeV   %.0f #pm %.0f", gp::E[i], area(i), areaErr(i)), "l");
      }
      leg->AddEntry(gN, Form("^{14}N blend   %.0f #pm %.0f", p[NL] * gp::sgN * sq / binW, F->GetParError(NL) * gp::sgN * sq / binW), "l");
      leg->AddEntry(gPS, Form("1n phase space   %.0f", psInt), "l");
      if (useBg) leg->AddEntry(gBg, "linear background", "l");
      if (gp::gs2) leg->AddEntry((TObject *)nullptr, Form("g.s. total   %.0f #pm %.0f", gsTot, gsTotErr), "");
      leg->Draw();

      pb->cd();
      pb->SetGridy();
      hPull->SetFillColor(TColor::GetColor("#898781")); hPull->SetLineColor(TColor::GetColor("#898781"));
      hPull->SetTitle(useBg ? ";E_{x} (MeV);(n#minusm)/#sqrt{m}" : ";E_{x} (MeV);deviance residual");
      double pr = withGS ? std::max(6.0, 1.1 * std::max(hPull->GetMaximum(), -hPull->GetMinimum())) : 6.0;
      hPull->SetMinimum(-pr); hPull->SetMaximum(pr);
      for (auto *ax : {hPull->GetXaxis(), hPull->GetYaxis()}) { ax->SetLabelSize(0.10); ax->SetTitleSize(0.11); }
      hPull->GetYaxis()->SetTitleOffset(0.55); hPull->GetYaxis()->SetTitleSize(0.09); hPull->GetYaxis()->SetNdivisions(505);
      hPull->Draw("hist");

      printf("%s  chi2/ndf %.1f/%d = %.2f  status %d cov %d\n", label[tag].Data(), pear, nfree - nFit, pear / (nfree - nFit),
             (int)r->Status(), r->CovMatrixStatus());
      printf("  BakerCousins %.1f/%d = %.2f ; bins with m<0.1 but n>0: %d\n", bc, nfree - nFit, bc / (nfree - nFit), nTiny);
      for (int i = 0; i < NL; ++i) printf("  LEVEL %.3f area %.1f err %.1f sigma %.4f\n", gp::E[i], area(i), areaErr(i), sigOf(i));
      if (withGS)
         printf("  GS mu %.4f +- %.4f  sigma %.4f +- %.4f\n", p[NL + 5], F->GetParError(NL + 5), p[NL + 4], F->GetParError(NL + 4));
      if (gp::gs2)
         printf("  GS broad area %.1f +- %.1f  mu %.4f +- %.4f  sigma %.4f +- %.4f ; g.s. total %.1f +- %.1f\n",
                areaAS(NL + 6, NL + 8), errAS(NL + 6, NL + 8), p[NL + 7], F->GetParError(NL + 7), p[NL + 8],
                F->GetParError(NL + 8), gsTot, gsTotErr);
      printf("  N14 %.1f  PS %.1f  BG %.4f + %.4f x\n", p[NL] * gp::sgN * sq / binW, psInt, p[NL + 2], p[NL + 3]);
      if (withGS) {
         double pg = 0, pi = 0, bg_ = 0, bi = 0, bm = 0; int ng = 0, ni = 0, nm = 0;
         for (int b = 1; b <= nb; ++b) {
            double x = h->GetBinCenter(b), m = std::max(gp::model(&x, p), 1e-12), n = h->GetBinContent(b);
            double dev = 2 * (m - n + (n > 0 ? n * std::log(n / m) : 0));
            if (x < 2.5) { pg += (n - m) * (n - m) / m; bg_ += dev; ++ng; }
            else if (x >= 5.0) { pi += (n - m) * (n - m) / m; bi += dev; ++ni; }
            else { bm += dev; ++nm; }
         }
         printf("  Ex<2.5: Pearson %.1f  BC %.1f / %d bins ; 2.5-5.0: BC %.1f / %d bins (counts %.0f) ; Ex>=5.0: Pearson %.1f  BC %.1f / %d bins\n",
                pg, bg_, ng, bm, nm, h->Integral(h->FindBin(2.501), h->FindBin(4.999)), pi, bi, ni);
      }
   }
   TString out = here + "/plots/globalfit_C14_" + which.ReplaceAll(",", "_vs_");
   if (withGS) out += gp::gs2 ? "_gs2" : "_gs";
   if (splitCm) out += Form("_gscm%g-%g", gsCmLo, gsCmHi);
   if (!useBg) out += "_nobg";
   c->SaveAs(out + ".png");
}
