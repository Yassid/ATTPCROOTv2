/// @file panels_pd.C
/// @brief All four 16C(p,d)15C states, one panel each, against the report's CRC calculations.
///
/// The curves are the COUPLED REACTION CHANNEL calculations from the theory report
/// ("Results of 16C(p,p') and (d,d') Analysis", 2025-09-05), digitised from the vector content of
/// its figures into pp/crc/ -- see pp/crc/README.txt. They REPLACE the FRESCO DWBA I built, which
/// used unexamined global optical potentials.
///
/// THEY ARE PLOTTED ABSOLUTE, WITH NO FITTED SCALE. The report fitted its <16C|15C+n> spectroscopic
/// amplitudes to the measured (p,d) angular distributions of the same 11.5 MeV/nucleon data set, so
/// the curves already carry their normalisation. Rescaling them onto our points would throw away
/// exactly the information worth having: the ratio of our absolute cross-section to theirs is an
/// independent check on our luminosity (168.3 mb^-1, itself carried over from the (p,p) elastic).
/// That ratio is printed per bin and summarised per state.
///
/// The old FRESCO comparison is still available in pp/dwba_compare_pd.C.
///
/// There is no curve for the ~4.9 structure. 15C is unbound above Sn = 1.218 MeV, so the
/// ~3.4 and ~4.9 structures are resonances: a transfer DWBA to a bound single-particle state does
/// not describe them, and FRESCO would need a continuum/resonance treatment (a bin in the n+14C
/// continuum, or an R-matrix form factor). Their panels therefore show DATA ONLY, and the panel
/// says so rather than leaving an empty frame to be misread as a failed calculation.
///
/// The two unbound cross-sections are Gaussian areas on top of the fitted continuum. With the
/// phase-space scale now pinned by the 1.4-2.9 MeV gap they are far better determined than they
/// were, but they remain lineshape-dependent in a way the bound states are not.
///
///   root -b -q 'pp/panels_pd.C()'

static TGraph *pnReadFresco(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   if (!in) return g;
   double a, x;
   while (in >> a >> x)
      if (x > 0) g->SetPoint(g->GetN(), a, x);
   return g;
}

void panels_pd(TString dataFile = "plots/fit_pd_ps_gap2.root",
               TString frescoDir = "/home/yassid/fair_install/16C_dp/fresco/outputs",
               TString stem = "pd15C_11MeV", Double_t fitLo = 20, Double_t fitHi = 80,
               TString png = "plots/panels_pd.png")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const char *tag[4] = {"gs", "ex1", "ex2", "ex3"};
   const char *ttl[4] = {"^{15}C g.s.  1/2^{+}", "^{15}C 0.740  5/2^{+}",
                         "^{15}C 3.103  1/2^{-}  (fitted at 3.38)", "^{15}C ~4.9  (no calculation)"};
   const char *crc[4] = {"gs", "e074", "e310", ""};
   int col[4] = {kBlack, kRed + 1, kGreen + 3, kOrange + 7};
   int mrk[4] = {20, 24, 21, 25};

   TFile *fd = TFile::Open(here + "/" + dataFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mcannot open %s\033[0m\n", dataFile.Data()); return; }
   TH1D *D[4], *S[4];
   for (int g = 0; g < 4; ++g) {
      D[g] = (TH1D *)fd->Get(Form("dsdo_%s", tag[g]));
      S[g] = (TH1D *)fd->Get(Form("dsdo_%s_syst", tag[g]));
      if (D[g]) { D[g] = (TH1D *)D[g]->Clone(Form("Dp%d", g)); D[g]->SetDirectory(nullptr); }
      if (S[g]) { S[g] = (TH1D *)S[g]->Clone(Form("Sp%d", g)); S[g]->SetDirectory(nullptr); }
   }
   fd->Close();

   // CRC curves, "Full" and "Reduced", for the first three states
   TGraph *T[4] = {nullptr, nullptr, nullptr, nullptr};   // Full
   TGraph *R[4] = {nullptr, nullptr, nullptr, nullptr};   // Reduced
   for (int g = 0; g < 4; ++g) {
      if (!crc[g][0]) continue;
      T[g] = pnReadFresco(here + "/crc/" + crc[g] + "_full.dat");
      R[g] = pnReadFresco(here + "/crc/" + crc[g] + "_reduced.dat");
      if (T[g]->GetN() < 10) { printf("\033[1;31m  missing crc/%s_full.dat\033[0m\n", crc[g]); T[g] = nullptr; }
   }

   TCanvas *c = new TCanvas("cpn", "pd panels", 1100, 850);
   c->Divide(2, 2, 0.001, 0.001);
   printf("\n  CRC curves plotted ABSOLUTE -- no scale fitted.\n");
   printf("  state      data/CRC     rms ln(data/CRC)  over %.0f-%.0f deg\n", fitLo, fitHi);
   for (int g = 0; g < 4; ++g) {
      c->cd(g + 1);
      gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
      if (!D[g]) { printf("  %-10s  no data histogram\n", tag[g]); continue; }

      double mx = 0, mn = 1e9;
      for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
         double v = D[g]->GetBinContent(b);
         if (v > 0) { mx = std::max(mx, v); mn = std::min(mn, v); }
      }
      auto *fr = new TH1F(Form("fr%d", g), Form("%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", ttl[g]),
                          100, 10, 85);
      fr->SetMinimum(std::max(0.05, mn * 0.35)); fr->SetMaximum(mx * 4.0);
      fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleSize(0.05);
      fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045);
      fr->Draw();

      if (S[g]) {
         auto *bnd = (TH1D *)D[g]->Clone(Form("bn%d", g));
         for (int b = 1; b <= bnd->GetNbinsX(); ++b) bnd->SetBinError(b, S[g]->GetBinError(b));
         bnd->SetFillColorAlpha(col[g], 0.28); bnd->SetLineColor(0); bnd->SetMarkerSize(0);
         bnd->Draw("E2 same");
      }

      double ratio = 0, rms = 0;
      if (T[g]) {
         // NO SCALE FITTED. Sum both over the overlap window and take the plain ratio.
         double sd = 0, st = 0, s2 = 0; int m = 0;
         for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
            double a = D[g]->GetBinCenter(b), d = D[g]->GetBinContent(b);
            if (d <= 0 || a < fitLo || a > fitHi) continue;
            double t = T[g]->Eval(a);
            if (t <= 0) continue;
            sd += d; st += t;
            double l = std::log(d / t); s2 += l * l; ++m;
         }
         ratio = st > 0 ? sd / st : 0;
         rms = m ? std::sqrt(s2 / m) : 0;
         for (int k = 0; k < 2; ++k) {
            TGraph *src = k ? R[g] : T[g];
            if (!src || src->GetN() < 10) continue;
            auto *tg = new TGraph();
            for (int i = 0; i < src->GetN(); ++i) {
               double x, y; src->GetPoint(i, x, y);
               if (x >= 10 && x <= 85) tg->SetPoint(tg->GetN(), x, y);
            }
            tg->SetLineColor(col[g]); tg->SetLineWidth(k ? 2 : 3); tg->SetLineStyle(k ? 2 : 1);
            tg->Draw("L same");
         }
         printf("  %-10s  %8.3f      %8.3f\n", tag[g], ratio, rms);
      } else {
         printf("  %-10s  no CRC curve in the report\n", tag[g]);
      }

      D[g]->SetMarkerStyle(mrk[g]); D[g]->SetMarkerColor(col[g]); D[g]->SetLineColor(col[g]);
      D[g]->SetMarkerSize(1.2); D[g]->SetLineWidth(2);
      D[g]->Draw("E1 same");

      auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.042);
      if (T[g]) {
         tx->DrawLatex(0.50, 0.86, Form("CRC: data/calc = %.2f", ratio));
         tx->SetTextSize(0.034); tx->SetTextColor(kGray + 2);
         tx->DrawLatex(0.50, 0.81, "full / reduced, absolute");
      } else {
         tx->SetTextColor(kGray + 2);
         tx->DrawLatex(0.40, 0.86, "no CRC curve");
         tx->DrawLatex(0.40, 0.81, "in the report");
      }
   }
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(here + "/" + png);
   printf("\n  wrote %s\n\n", png.Data());
}
