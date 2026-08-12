/// @file dwba_compare_pd.C
/// @brief 16C(p,d)15C measured angular distributions against a FRESCO transfer DWBA.
///
/// The DWBA is NOT one of the files in the Downloads DWBA folder -- those are elastic only
/// (16C_pp_*, 16C_dd_*). A transfer calculation had to be built: 16C_dp/fresco/inputs/
/// pd15C_11MeV.nin, adapted from the (d,p) template in the same directory.
///     entrance   p + 16C, Koning-Delaroche at Elab(p) = 11.636 MeV
///                (= 16C beam 185 MeV in inverse kinematics)
///     exit       d + 15C, Daehnick "B" at Ed = 10.1 MeV
///     Q          -2.0257 MeV  =  B(d) - Sn(16C)  =  2.2246 - 4.2503
///     overlaps   <16C(0+)|15C(1/2+) + n[2s1/2]>  l=0  BE 4.2503
///                <16C(0+)|15C(5/2+) + n[1d5/2]>  l=2  BE 4.9903
///
/// THE TWO STATES CARRY DIFFERENT L BY ANGULAR MOMENTUM ALONE. 16C is 0+, so coupling to 15C(1/2+)
/// forces j=1/2 on the picked-up neutron and hence l=0; coupling to 15C(5/2+) forces j=5/2 and
/// l=2. So the two measured shapes SHOULD differ, and that is the prediction being tested here --
/// it is not something put in by hand.
///
/// SPECTROSCOPIC FACTORS ARE NOT IN THE CALCULATION. The cfp amplitudes are unity, so each curve is
/// fitted to the data with ONE free scale, and that scale IS the spectroscopic factor (times any
/// residual normalisation). It is reported, but read it as indicative: it inherits the ~4 %
/// luminosity softness and every optical-model choice above.
///
///   root -b -q 'pp/dwba_compare_pd.C()'

/// run_fresco.sh already writes two-column "angle mb/sr" files.
static TGraph *readFresco(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   if (!in) { printf("\033[1;31mcannot open %s\033[0m\n", f.Data()); return g; }
   double a, x;
   while (in >> a >> x)
      if (x > 0) g->SetPoint(g->GetN(), a, x);
   return g;
}

void dwba_compare_pd(TString dataFile = "plots/fit_pd_ps_gap2.root",
                     TString frescoDir = "/home/yassid/fair_install/16C_dp/fresco/outputs",
                     TString stem = "pd15C_11MeV", Double_t fitLo = 20, Double_t fitHi = 80,
                     TString png = "plots/dwba_compare_pd.png")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + dataFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mcannot open %s\033[0m\n", dataFile.Data()); return; }
   TH1D *D[2] = {(TH1D *)fd->Get("dsdo_gs"), (TH1D *)fd->Get("dsdo_ex1")};
   TH1D *S[2] = {(TH1D *)fd->Get("dsdo_gs_syst"), (TH1D *)fd->Get("dsdo_ex1_syst")};
   for (int g = 0; g < 2; ++g) {
      if (!D[g]) { printf("\033[1;31mmissing dsdo histogram %d\033[0m\n", g); return; }
      D[g] = (TH1D *)D[g]->Clone(Form("D%d", g)); D[g]->SetDirectory(nullptr);
      if (S[g]) { S[g] = (TH1D *)S[g]->Clone(Form("S%d", g)); S[g]->SetDirectory(nullptr); }
   }
   fd->Close();

   // FILE NAMES ARE MISLEADING. run_fresco.sh labels channels by the outgoing particle assuming a
   // (d,p) reaction: an outgoing proton is "transfer", an outgoing deuteron is "elastic". For (p,d)
   // that is exactly backwards, so _state1.dat holds the p+16C ELASTIC (6.7e7 mb/sr at 1 deg, the
   // Rutherford divergence) and _dsdo.dat holds the g.s. TRANSFER. Taken from FRESCO's own channel
   // headers in the .out, not from the file names.
   TGraph *T[2] = {readFresco(frescoDir + "/" + stem + "_dsdo.dat"),     // g.s. transfer, l=0
                   readFresco(frescoDir + "/" + stem + "_state2.dat")};  // 0.740 transfer, l=2
   const char *nm[2] = {"g.s. (1/2^{+}), #font[12]{l} = 0", "0.740 (5/2^{+}), #font[12]{l} = 2"};
   const char *sn[2] = {"g.s.  l=0", "0.740 l=2"};
   for (int g = 0; g < 2; ++g)
      if (T[g]->GetN() < 10) { printf("\033[1;31mFRESCO curve %d has %d points\033[0m\n", g, T[g]->GetN()); return; }

   printf("\n  FRESCO: %d and %d angles\n", T[0]->GetN(), T[1]->GetN());
   double scale[2] = {1, 1}, rms[2] = {0, 0};
   for (int g = 0; g < 2; ++g) {
      // error-weighted least-squares scale over the fit window; stat and syst added in quadrature
      double num = 0, den = 0;
      for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
         double a = D[g]->GetBinCenter(b), d = D[g]->GetBinContent(b);
         if (d <= 0 || a < fitLo || a > fitHi) continue;
         double e = D[g]->GetBinError(b), sy = S[g] ? S[g]->GetBinError(b) : 0;
         double w = 1.0 / (e * e + sy * sy);
         double t = T[g]->Eval(a);
         if (t <= 0) continue;
         num += w * d * t; den += w * t * t;
      }
      if (den > 0) scale[g] = num / den;
      double s2 = 0; int m = 0;
      for (int b = 1; b <= D[g]->GetNbinsX(); ++b) {
         double a = D[g]->GetBinCenter(b), d = D[g]->GetBinContent(b);
         if (d <= 0 || a < fitLo || a > fitHi) continue;
         double t = scale[g] * T[g]->Eval(a);
         if (t <= 0) continue;
         double l = std::log(d / t); s2 += l * l; ++m;
      }
      rms[g] = m ? std::sqrt(s2 / m) : 0;
      printf("  %-10s  scale (~ S-factor) = %.4g    rms ln(data/theory) over %.0f-%.0f = %.3f\n", sn[g], scale[g],
             fitLo, fitHi, rms[g]);
   }

   printf("\n  theta_cm |   g.s. data / DWBA   |  0.740 data / DWBA\n");
   for (int b = 1; b <= D[0]->GetNbinsX(); ++b) {
      double a = D[0]->GetBinCenter(b);
      if (D[0]->GetBinContent(b) <= 0 && D[1]->GetBinContent(b) <= 0) continue;
      printf("  %5.1f    ", a);
      for (int g = 0; g < 2; ++g) {
         double d = D[g]->GetBinContent(b), t = scale[g] * T[g]->Eval(a);
         if (d > 0 && t > 0) printf("  %6.2f / %-6.2f = %4.2f  ", d, t, d / t);
         else printf("        --            ");
      }
      printf("\n");
   }

   TCanvas *c = new TCanvas("cpd", "pd DWBA", 1000, 760);
   gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   int col[2] = {kBlack, kRed + 1}, mk[2] = {20, 24};
   auto *fr = new TH1F("fr", "^{16}C(p,d)^{15}C  vs  FRESCO transfer DWBA;#theta_{cm} [deg];"
                             "d#sigma/d#Omega [mb/sr]", 100, 10, 90);
   fr->SetMinimum(0.3); fr->SetMaximum(30); fr->Draw();
   auto *lg = new TLegend(0.55, 0.70, 0.89, 0.88);
   lg->SetFillStyle(0);
   for (int g = 0; g < 2; ++g) {
      if (S[g]) { // systematic band behind the points
         auto *bnd = (TH1D *)D[g]->Clone(Form("b%d", g));
         for (int b = 1; b <= bnd->GetNbinsX(); ++b) bnd->SetBinError(b, S[g]->GetBinError(b));
         bnd->SetFillColorAlpha(col[g], 0.25); bnd->SetLineColor(0); bnd->SetMarkerSize(0);
         bnd->Draw("E2 same");
      }
      auto *tg = new TGraph();
      for (int i = 0; i < T[g]->GetN(); ++i) {
         double x, y; T[g]->GetPoint(i, x, y);
         if (x >= 10 && x <= 90) tg->SetPoint(tg->GetN(), x, scale[g] * y);
      }
      tg->SetLineColor(col[g]); tg->SetLineWidth(3); tg->SetLineStyle(g ? 2 : 1); tg->Draw("L same");
      D[g]->SetMarkerStyle(mk[g]); D[g]->SetMarkerColor(col[g]); D[g]->SetLineColor(col[g]);
      D[g]->SetMarkerSize(1.3); D[g]->SetLineWidth(2);
      D[g]->Draw("E1 same");
      lg->AddEntry(D[g], Form("%s   S ~ %.2f, rms %.2f", nm[g], scale[g], rms[g]), "lp");
   }
   lg->Draw();
   auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.026);
   tx->DrawLatex(0.14, 0.17, "FRESCO: KD p+^{16}C at 11.64 MeV, Daehnick d+^{15}C, Q = -2.03 MeV");
   tx->DrawLatex(0.14, 0.135, "bands = width-treatment systematic; one free scale per state");
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(here + "/" + png);
   printf("\n  wrote %s\n\n", png.Data());
}
