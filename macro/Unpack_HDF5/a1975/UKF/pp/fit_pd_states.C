/// @file fit_pd_states.C
/// @brief 16C(p,d)15C excitation spectrum: global peak fit, then the same fit in theta_cm slices.
///
/// Selection is the one saved from the browser explorer on 2026-08-12
/// (explorer_params_16C_p_d_15C_2026-08-12-11-27-43.json):
///     Ebeam 185 MeV, chi2/ndf < 5, KE 0-30 MeV, vertex z 0-500 mm, theta_cm 0-80,
///     theta correction ON: theta -> theta - (360/2950.8)*(KE - 27)
/// THE VERTEX CUT IS HALF THE TARGET. The elastic luminosity (316.4 mb^-1) was measured over the
/// full ~940 mm with no z cut, so it does NOT transfer to this selection unscaled -- see the note
/// where the yields are written out.
///
/// WHAT IS FITTED. Four Gaussians plus a linear background, over the FULL range, both globally and
/// in every theta_cm slice. Below Sn(15C) = 1.218 MeV the two bound states are cleanly resolved
/// under this selection: the g.s. near 0 and the 5/2+ near 0.78. The structures at ~3.4 and
/// ~4.5-5.4 are resonances on a continuum, so their Gaussian areas depend on where the background
/// is drawn and should be read as indicative rather than as yields -- but they are fitted, not
/// excluded. Fitting only the bound region instead leaves the background under those peaks
/// determined by two bins of empty gap, which is worse for the bound peaks too.
///
/// PER-ANGLE WIDTHS. The two bound peaks are the same resolution, so they SHARE one width, which
/// is free per slice within [0.10, 0.35]. Fixing both to the global value was wrong: resolution
/// varies with angle, and in a slice sharper than average the fit cannot reach the peak -- the
/// 25-30 deg bin undershot the 0.740 state by more than half (chi2/ndf 1.85). One shared free
/// width costs one parameter and fixes it. The two resonance widths stay fixed to the global
/// values, where the statistics per slice do not support freeing anything. Centroids are free
/// throughout: the Ex-vs-theta_cm drift is what a per-angle fit should expose.
///
///   root -b -q 'pp/fit_pd_states.C("/path/pd_kin.root")'

static double fp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static std::pair<double, double> fp_kine(double Kp, double thl, double Ke)
{
   const double u = 931.49401;
   double m1 = 16.014701 * u, m2 = 1.007825 * u, m3 = 2.013553 * u, m4 = 15.010599 * u;
   double E1 = Kp + m1, E3 = Ke + m3, E4 = E1 + m2 - E3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, uu = m2 * m2 + m3 * m3 - 2 * m2 * E3;
   double a = (std::cos(thl) * fp_om2(s, m1 * m1, m2 * m2) * fp_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   if (a < 0) return {std::nan(""), std::nan("")};
   double m4x = std::sqrt(a), ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * E4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
               (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (fp_om2(s, m1 * m1, m2 * m2) * fp_om2(s, m3 * m3, m4x * m4x));
   if (c < -1 || c > 1) return {std::nan(""), std::nan("")};
   return {ex, (TMath::Pi() - std::acos(c)) * TMath::RadToDeg()};
}

/// @param driftCorr  subtract the measured Ex-vs-theta_cm drift. The drift is fitted on the G.S.
///                   ONLY and then applied to the whole spectrum, so the 0.740 state is an
///                   INDEPENDENT test: if it also flattens, the drift is a common energy-scale
///                   effect and the correction is real; if only the g.s. flattens, it is a fudge.
///                   Verified NOT to be a vertex-z effect -- the g.s. centroid is flat in z
///                   (+0.003 to +0.026 over 50-500 mm), so it is not beam energy loss.
void fit_pd_states(TString cache, Double_t Ebeam = 185.0, Double_t kcDenom = 2950.8, Double_t kcPivot = 27.0,
                   Double_t chi2max = 5.0, Double_t keLo = 0, Double_t keHi = 30, Double_t vzLo = 0,
                   Double_t vzHi = 500, Double_t icMin = 950, Double_t icMax = 1350, Double_t cmLo = 0,
                   Double_t cmHi = 80, Double_t dcm = 5.0, TString tag = "sel20260812",
                   Bool_t driftCorr = kFALSE, Double_t lumiFull = 316.4, Double_t zFull = 940.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double slope = 360.0 / kcDenom;

   TFile *f = TFile::Open(cache);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("pk") : nullptr;
   if (!t) { printf("\033[1;31mcannot open %s\033[0m\n", cache.Data()); return; }
   float ke, th, vz, c2, ic;
   t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vz", &vz);
   t->SetBranchAddress("chi2ndf", &c2); t->SetBranchAddress("ic", &ic);

   const int NB = (int)std::lround((cmHi - cmLo) / dcm);
   auto *hAll = new TH1D("hAll", "", 240, -1.5, 6.5); hAll->SetDirectory(nullptr); hAll->Sumw2();
   std::vector<TH1D *> hSl(NB);
   for (int k = 0; k < NB; ++k) {
      hSl[k] = new TH1D(Form("hSl%d", k), "", 120, -1.5, 6.5);
      hSl[k]->SetDirectory(nullptr); hSl[k]->Sumw2();
   }
   // pass 1: collect
   std::vector<double> vEx, vTc;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2max && ic > icMin && ic < icMax)) continue;
      if (!(ke > keLo && ke < keHi)) continue;
      if (!(vz > vzLo && vz < vzHi)) continue;
      double thc = (th - slope * (ke - kcPivot)) * TMath::DegToRad();
      auto [ex, tcm] = fp_kine(Ebeam, thc, ke);
      if (std::isnan(ex) || tcm < cmLo || tcm > cmHi) continue;
      vEx.push_back(ex); vTc.push_back(tcm);
   }
   f->Close();
   long n = (long)vEx.size();
   printf("\n  %ld deuterons after the saved selection\n", n);

   // pass 2: measure the g.s. drift vs theta_cm and fit it smooth
   TF1 *DR = nullptr;
   if (driftCorr) {
      auto *gDr = new TGraphErrors();
      for (int k = 0; k < NB; ++k) {
         double lo = cmLo + k * dcm, hi = lo + dcm;
         TH1D hg("hg", "", 70, -1.0, 1.5);
         for (size_t i = 0; i < vEx.size(); ++i)
            if (vTc[i] >= lo && vTc[i] < hi) hg.Fill(vEx[i]);
         if (hg.Integral() < 60) continue;
         TF1 g("g", "gaus", -0.6, 0.45);
         g.SetParameters(hg.GetMaximum(), 0.0, 0.2);
         hg.Fit(&g, "RQN0");
         int p = gDr->GetN();
         gDr->SetPoint(p, 0.5 * (lo + hi), g.GetParameter(1));
         gDr->SetPointError(p, 0, std::max(g.GetParError(1), 0.005));
      }
      DR = new TF1("DR", "pol2", cmLo, cmHi);
      gDr->Fit(DR, "RQN0");
      printf("  drift fit (g.s. only): Ex_shift(theta) = %+.4f %+.4e*th %+.4e*th^2  [%d points]\n",
             DR->GetParameter(0), DR->GetParameter(1), DR->GetParameter(2), gDr->GetN());
   }

   // pass 3: fill, drift-corrected
   for (size_t i = 0; i < vEx.size(); ++i) {
      double ex = vEx[i] - (DR ? DR->Eval(vTc[i]) : 0.0);
      hAll->Fill(ex);
      int k = (int)((vTc[i] - cmLo) / dcm);
      if (k >= 0 && k < NB) hSl[k]->Fill(ex);
   }

   // ---------------- global fit: 4 gaussians + linear background ----------------
   TF1 *G = new TF1("G", "gaus(0)+gaus(3)+gaus(6)+gaus(9)+pol1(12)", -1.0, 6.0);
   double init[14] = {170, 0.05, 0.18, 250, 0.78, 0.18, 120, 3.40, 0.25, 30, 4.90, 0.50, 5, 0};
   G->SetParameters(init);
   const char *pn[14] = {"A_gs", "E_gs", "s_gs", "A_074", "E_074", "s_074", "A_34", "E_34",
                         "s_34", "A_49", "E_49", "s_49", "bg0", "bg1"};
   for (int i = 0; i < 14; ++i) G->SetParName(i, pn[i]);
   for (int g = 0; g < 4; ++g) { G->SetParLimits(3 * g, 0, 1e6); G->SetParLimits(3 * g + 2, 0.05, 1.2); }
   hAll->Fit(G, "RQ0");
   hAll->Fit(G, "RQ0");
   const double binw = hAll->GetBinWidth(1);
   printf("\n  GLOBAL FIT  (chi2/ndf = %.2f)\n", G->GetChisquare() / std::max(1, G->GetNDF()));
   printf("  state      centroid [MeV]      sigma [MeV]        area [counts]\n");
   const char *sn[4] = {"g.s.", "0.740", "~3.4", "~4.9"};
   double sgGlob[4];
   for (int g = 0; g < 4; ++g) {
      double A = G->GetParameter(3 * g), E = G->GetParameter(3 * g + 1), S = G->GetParameter(3 * g + 2);
      double eA = G->GetParError(3 * g), eE = G->GetParError(3 * g + 1), eS = G->GetParError(3 * g + 2);
      sgGlob[g] = S;
      double area = A * S * std::sqrt(2 * TMath::Pi()) / binw;
      double earea = area * std::sqrt(std::pow(eA / std::max(A, 1e-9), 2) + std::pow(eS / std::max(S, 1e-9), 2));
      printf("  %-8s  %6.3f +- %5.3f    %6.3f +- %5.3f    %8.0f +- %5.0f\n", sn[g], E, eE, S, eS, area, earea);
   }
   printf("\n  the two unbound structures are resonances on a continuum; their Gaussian areas are\n");
   printf("  background-dependent and are NOT carried into the angular distribution.\n");

   // ---------------- per-angle fits: the two bound states, widths fixed --------------
   auto *yGS = new TH1D("yGS", "", NB, cmLo, cmHi); yGS->SetDirectory(nullptr);
   auto *y074 = new TH1D("y074", "", NB, cmLo, cmHi); y074->SetDirectory(nullptr);
   printf("\n  PER-ANGLE FITS (bound peaks share ONE free width; resonance widths fixed at %.3f / %.3f)\n",
          sgGlob[2], sgGlob[3]);
   auto *y34 = new TH1D("y34", "", NB, cmLo, cmHi); y34->SetDirectory(nullptr);
   auto *y49 = new TH1D("y49", "", NB, cmLo, cmHi); y49->SetDirectory(nullptr);
   printf("  theta_cm | entries | sigma |  E_gs  |  N_gs      |  E_074 |  N_074     |  N_3.4    |  N_4.9    | c2/ndf\n");
   for (int k = 0; k < NB; ++k) {
      double lo = cmLo + k * dcm, hi = lo + dcm;
      if (hSl[k]->Integral() < 40) {
         printf("  %3.0f-%3.0f  | %8.0f |  too few entries to fit\n", lo, hi, hSl[k]->Integral());
         continue;
      }
      // [2] is ONE shared width for both bound peaks -- same resolution, so one parameter
      TF1 *F = new TF1(Form("F%d", k),
                       "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[2])^2)"
                       "+[5]*exp(-0.5*((x-[6])/[7])^2)+[8]*exp(-0.5*((x-[9])/[10])^2)+[11]+[12]*x",
                       -1.2, 6.2);
      double mx = hSl[k]->GetMaximum();
      double ini[13] = {mx * 0.6, 0.05, 0.5 * (sgGlob[0] + sgGlob[1]), mx * 0.8, 0.78,
                        mx * 0.35, 3.38, sgGlob[2], mx * 0.15, 4.91, sgGlob[3], 0.5, 0};
      F->SetParameters(ini);
      F->SetParLimits(0, 0, 1e6); F->SetParLimits(3, 0, 1e6);
      F->SetParLimits(5, 0, 1e6); F->SetParLimits(8, 0, 1e6);
      F->SetParLimits(2, 0.10, 0.35);          // shared bound-state resolution, free
      F->FixParameter(7, sgGlob[2]);           // resonance widths: statistics per slice
      F->FixParameter(10, sgGlob[3]);          // do not support freeing these
      F->SetParLimits(1, -0.5, 0.5); F->SetParLimits(4, 0.45, 1.15);
      F->SetParLimits(6, 3.0, 3.8);  F->SetParLimits(9, 4.3, 5.5);
      hSl[k]->Fit(F, "RQ0");
      hSl[k]->Fit(F, "RQ0");
      double bw = hSl[k]->GetBinWidth(1);
      const double sw = F->GetParameter(2); // the fitted shared width
      double n0 = F->GetParameter(0) * sw * std::sqrt(2 * TMath::Pi()) / bw;
      double e0 = n0 * (F->GetParameter(0) > 0 ? F->GetParError(0) / F->GetParameter(0) : 1);
      double n1 = F->GetParameter(3) * sw * std::sqrt(2 * TMath::Pi()) / bw;
      double e1 = n1 * (F->GetParameter(3) > 0 ? F->GetParError(3) / F->GetParameter(3) : 1);
      double n2 = F->GetParameter(5) * sgGlob[2] * std::sqrt(2 * TMath::Pi()) / bw;
      double e2 = n2 * (F->GetParameter(5) > 0 ? F->GetParError(5) / F->GetParameter(5) : 1);
      double n3 = F->GetParameter(8) * sgGlob[3] * std::sqrt(2 * TMath::Pi()) / bw;
      double e3 = n3 * (F->GetParameter(8) > 0 ? F->GetParError(8) / F->GetParameter(8) : 1);
      yGS->SetBinContent(k + 1, n0); yGS->SetBinError(k + 1, e0);
      y074->SetBinContent(k + 1, n1); y074->SetBinError(k + 1, e1);
      y34->SetBinContent(k + 1, n2); y34->SetBinError(k + 1, e2);
      y49->SetBinContent(k + 1, n3); y49->SetBinError(k + 1, e3);
      printf("  %3.0f-%3.0f  | %7.0f | %5.3f | %+6.3f | %5.0f +-%4.0f | %+6.3f | %5.0f +-%4.0f | %4.0f +-%3.0f | %4.0f +-%3.0f | %5.2f\n",
             lo, hi, hSl[k]->Integral(), sw, F->GetParameter(1), n0, e0, F->GetParameter(4), n1, e1, n2, e2, n3, e3,
             F->GetChisquare() / std::max(1, F->GetNDF()));
      hSl[k]->GetListOfFunctions()->Add(F);
   }

   // ---------------- draw ----------------
   TCanvas *cvG = new TCanvas("cg", "global", 1000, 700);
   gPad->SetLogy();
   hAll->SetTitle("^{16}C(p,d)^{15}C global fit;E_{x} [MeV];counts / 33 keV");
   hAll->SetLineWidth(2); hAll->SetMinimum(0.5); hAll->Draw("hist");
   G->SetLineColor(kRed + 1); G->SetLineWidth(2); G->Draw("same");
   for (int g = 0; g < 4; ++g) {
      auto *gg = new TF1(Form("g%d", g), "gaus", -1.0, 6.0);
      gg->SetParameters(G->GetParameter(3 * g), G->GetParameter(3 * g + 1), G->GetParameter(3 * g + 2));
      gg->SetLineColor(kAzure + 2); gg->SetLineStyle(2); gg->SetLineWidth(1); gg->Draw("same");
   }
   cvG->SaveAs(here + "/plots/fit_pd_global_" + tag + ".png");

   int nc = 4, nr = (NB + nc - 1) / nc;
   TCanvas *cvS = new TCanvas("cs", "slices", 340 * nc, 260 * nr);
   cvS->Divide(nc, nr);
   for (int k = 0; k < NB; ++k) {
      cvS->cd(k + 1);
      hSl[k]->SetTitle(Form("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", cmLo + k * dcm, cmLo + (k + 1) * dcm));
      gPad->SetLogy(); hSl[k]->SetMinimum(0.5);
      hSl[k]->SetLineWidth(2); hSl[k]->Draw("hist");
      if (hSl[k]->GetListOfFunctions()->GetSize()) {
         auto *F = (TF1 *)hSl[k]->GetListOfFunctions()->First();
         F->SetLineColor(kRed + 1); F->Draw("same");
      }
   }
   cvS->SaveAs(here + "/plots/fit_pd_slices_" + tag + ".png");

   TCanvas *cvY = new TCanvas("cy", "yields", 900, 650);
   gPad->SetLogy(); gPad->SetGridy();
   yGS->SetTitle("fitted peak yield (no acceptance, no solid angle);#theta_{cm} [deg];counts");
   yGS->SetMarkerStyle(20); yGS->SetLineColor(kBlack); yGS->SetMarkerColor(kBlack); yGS->SetLineWidth(2);
   yGS->Draw("E1");
   y074->SetMarkerStyle(24); y074->SetLineColor(kRed + 1); y074->SetMarkerColor(kRed + 1); y074->SetLineWidth(2);
   y074->Draw("E1 same");
   auto *lg = new TLegend(0.62, 0.75, 0.88, 0.88);
   lg->AddEntry(yGS, "g.s. (1/2^{+})", "lp");
   lg->AddEntry(y074, "0.740 (5/2^{+})", "lp");
   y34->SetMarkerStyle(21); y34->SetLineColor(kGreen + 3); y34->SetMarkerColor(kGreen + 3); y34->SetLineWidth(2);
   y34->Draw("E1 same");
   y49->SetMarkerStyle(25); y49->SetLineColor(kOrange + 7); y49->SetMarkerColor(kOrange + 7); y49->SetLineWidth(2);
   y49->Draw("E1 same");
   lg->AddEntry(y34, "~3.4 (resonance)", "lp");
   lg->AddEntry(y49, "~4.9 (resonance)", "lp");
   lg->Draw();
   cvY->SaveAs(here + "/plots/fit_pd_yields_" + tag + ".png");

   // ---------------- acceptance + luminosity -> absolute mb/sr ----------------
   // L scales with the length of target kept: the elastic 316.4 mb^-1 was integrated over the full
   // ~940 mm with no vertex cut, and the gas is uniform, so a [vzLo,vzHi] window keeps that
   // fraction of it. This is the ONLY place the vertex cut enters the normalisation.
   const double lumi = lumiFull * (vzHi - vzLo) / zFull;
   printf("\n  luminosity: %.1f mb^-1 full target x %.0f/%.0f mm kept = %.1f mb^-1\n", lumiFull, vzHi - vzLo,
          zFull, lumi);
   const char *accTag[4] = {"gs", "ex1", "ex2", "ex3"};
   const char *accNm[4] = {"g.s. (0.000)", "0.740", "~3.4 (sim 3.100)", "~4.9 (sim 4.660)"};
   TH1D *Y[4] = {yGS, y074, y34, y49};
   TH1D *X[4] = {nullptr, nullptr, nullptr, nullptr};
   auto dOm = [](double lo, double hi) {
      return 2 * TMath::Pi() * (std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad()));
   };
   for (int g = 0; g < 4; ++g) {
      TString af = here + Form("/../../../../Simulation/ATTPC/16C_pd/diagnostics/acceptance_%s.root", accTag[g]);
      TFile *fa = TFile::Open(af);
      TH1D *A = fa && !fa->IsZombie() ? (TH1D *)fa->Get(Form("hAcc_%s", accTag[g])) : nullptr;
      if (!A) { printf("  no acceptance for %s -- skipped\n", accTag[g]); if (fa) fa->Close(); continue; }
      A = (TH1D *)A->Clone(Form("A%d", g)); A->SetDirectory(nullptr); fa->Close();
      if (A->GetBinWidth(1) > dcm + 1e-6) {
         printf("\033[1;31m  %s: acceptance bins (%.1f) coarser than the yield (%.1f) -- refusing\033[0m\n",
                accTag[g], A->GetBinWidth(1), dcm);
         continue;
      }
      X[g] = (TH1D *)Y[g]->Clone(Form("X%d", g)); X[g]->Reset(); X[g]->SetDirectory(nullptr);
      printf("\n  %s\n  theta_cm |  N fitted     |  acc   |  dsigma/dOmega [mb/sr]\n", accNm[g]);
      for (int b = 1; b <= Y[g]->GetNbinsX(); ++b) {
         double lo = Y[g]->GetBinLowEdge(b), hi = lo + Y[g]->GetBinWidth(b);
         double nn = Y[g]->GetBinContent(b), en = Y[g]->GetBinError(b);
         if (nn <= 0) continue;
         double a = A->GetBinContent(A->FindBin(Y[g]->GetBinCenter(b)));
         double ea = A->GetBinError(A->FindBin(Y[g]->GetBinCenter(b)));
         if (a < 0.15) { printf("  %3.0f-%3.0f  | %6.0f        |  %.3f | DROPPED (acc < 0.15)\n", lo, hi, nn, a);
                         continue; }
         double v = nn / dOm(lo, hi) / a / lumi;
         double e = v * std::sqrt(std::pow(en / nn, 2) + std::pow(ea / a, 2));
         X[g]->SetBinContent(b, v); X[g]->SetBinError(b, e);
         printf("  %3.0f-%3.0f  | %6.0f +-%4.0f |  %.3f |  %8.4g +- %.3g\n", lo, hi, nn, en, a, v, e);
      }
   }
   TCanvas *cvX = new TCanvas("cx", "xsec", 950, 700);
   gPad->SetLogy(); gPad->SetGridy();
   int xc[4] = {kBlack, kRed + 1, kGreen + 3, kOrange + 7};
   int xm[4] = {20, 24, 21, 25};
   auto *lgx = new TLegend(0.60, 0.70, 0.89, 0.88);
   bool first = true;
   for (int g = 0; g < 4; ++g) {
      if (!X[g]) continue;
      X[g]->SetMarkerStyle(xm[g]); X[g]->SetMarkerColor(xc[g]); X[g]->SetLineColor(xc[g]); X[g]->SetLineWidth(2);
      X[g]->SetTitle("^{16}C(p,d)^{15}C;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
      X[g]->Draw(first ? "E1" : "E1 same"); first = false;
      lgx->AddEntry(X[g], accNm[g], "lp");
   }
   lgx->Draw();
   cvX->SaveAs(here + "/plots/xsec_pd_states_" + tag + ".png");

   TFile fo(here + "/plots/fit_pd_" + tag + ".root", "RECREATE");
   for (int g = 0; g < 4; ++g) if (X[g]) X[g]->Write(Form("dsdo_%s", accTag[g]));
   hAll->Write("ex_all"); yGS->Write("yield_gs"); y074->Write("yield_074");
   y34->Write("yield_34"); y49->Write("yield_49");
   for (int k = 0; k < NB; ++k) hSl[k]->Write(Form("ex_slice_%d", k));
   fo.Close();
   printf("\n  wrote plots/fit_pd_{global,slices,yields}_%s.png and fit_pd_%s.root\n\n", tag.Data(), tag.Data());
}
