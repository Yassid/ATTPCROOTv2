/// @file panels_pd.C
/// @brief All four 16C(p,d)15C states, one panel each, with the DWBA where it exists.
///
/// The calculation covers ONLY the two bound states. 15C is unbound above Sn = 1.218 MeV, so the
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
   const char *ttl[4] = {"^{15}C g.s.  1/2^{+}   (#font[12]{l} = 0)", "^{15}C 0.740  5/2^{+}   (#font[12]{l} = 2)",
                         "^{15}C ~3.4  (unbound)", "^{15}C ~4.9  (unbound)"};
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

   // DWBA only for the two bound states. NOTE the file naming: run_fresco.sh assumes a (d,p)
   // reaction and labels by outgoing particle, so for (p,d) "_dsdo" is the g.s. TRANSFER and
   // "_state1" is the elastic.
   TGraph *T[4] = {pnReadFresco(frescoDir + "/" + stem + "_dsdo.dat"),
                   pnReadFresco(frescoDir + "/" + stem + "_state2.dat"), nullptr, nullptr};

   TCanvas *c = new TCanvas("cpn", "pd panels", 1100, 850);
   c->Divide(2, 2, 0.001, 0.001);
   printf("\n  state        scale (~S)   rms ln(data/theory) %.0f-%.0f\n", fitLo, fitHi);
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

      double sc = 1, rms = 0;
      if (T[g] && T[g]->GetN() > 10) {
         double num = 0, den = 0;
         for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
            double a = D[g]->GetBinCenter(b), d = D[g]->GetBinContent(b);
            if (d <= 0 || a < fitLo || a > fitHi) continue;
            double e = D[g]->GetBinError(b), sy = S[g] ? S[g]->GetBinError(b) : 0;
            double w = 1.0 / (e * e + sy * sy), t = T[g]->Eval(a);
            if (t <= 0) continue;
            num += w * d * t; den += w * t * t;
         }
         if (den > 0) sc = num / den;
         double s2 = 0; int m = 0;
         for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
            double a = D[g]->GetBinCenter(b), d = D[g]->GetBinContent(b);
            if (d <= 0 || a < fitLo || a > fitHi) continue;
            double t = sc * T[g]->Eval(a);
            if (t <= 0) continue;
            double l = std::log(d / t); s2 += l * l; ++m;
         }
         rms = m ? std::sqrt(s2 / m) : 0;
         auto *tg = new TGraph();
         for (int i = 0; i < T[g]->GetN(); ++i) {
            double x, y; T[g]->GetPoint(i, x, y);
            if (x >= 10 && x <= 85) tg->SetPoint(tg->GetN(), x, sc * y);
         }
         tg->SetLineColor(col[g]); tg->SetLineWidth(3); tg->Draw("L same");
         printf("  %-10s  %10.4g   %8.3f\n", tag[g], sc, rms);
      } else {
         printf("  %-10s  no DWBA (unbound - needs a continuum/resonance treatment)\n", tag[g]);
      }

      D[g]->SetMarkerStyle(mrk[g]); D[g]->SetMarkerColor(col[g]); D[g]->SetLineColor(col[g]);
      D[g]->SetMarkerSize(1.2); D[g]->SetLineWidth(2);
      D[g]->Draw("E1 same");

      auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.042);
      if (T[g] && T[g]->GetN() > 10)
         tx->DrawLatex(0.55, 0.85, Form("DWBA  S #approx %.3f", sc));
      else {
         tx->SetTextColor(kGray + 2);
         tx->DrawLatex(0.36, 0.86, "no DWBA: unbound,");
         tx->DrawLatex(0.36, 0.81, "needs a continuum");
         tx->DrawLatex(0.36, 0.76, "treatment");
      }
   }
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(here + "/" + png);
   printf("\n  wrote %s\n\n", png.Data());
}
