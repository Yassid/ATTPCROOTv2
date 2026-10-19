/// @file fit_vz_gs_dp.C -- 2026-10-19. TRUE g.s. yield vs vz, data vs far sim, Spyral kinematics (census_arc trees).
/// Per theta_cm bin (10-20, 20-30, 30-40): Gaussian mean/sigma FIXED from a fit of that angle bin over all vz 200-980
/// (gaussian + linear bg, Poisson L), then per vz slice only amplitude + linear bg float. 100 keV bins, Ex -3..2.
/// vz 980-1010 EXCLUDED: data spike there is Ex-flat background (far_ex_dp.C, 2026-10-18).
/// Sim (pure g.s.) yield = counts in the bound window -1.5..0.729. Normalisation: data/sim on 30-40, vz 200-900.
/// Run from D2_UKF:  root -l -b -q fit_vz_gs_dp.C   -> fit_vz_gs_dp.png, fit_vz_gs_dp_slices.png, stdout table
void fit_vz_gs_dp()
{
   gStyle->SetOptStat(0);
   TChain d("c"), s("c");
   d.Add("/mnt/f/a1975/census_arc/data_*.root");
   s.Add("/mnt/f/a1975_C16_dp_sim_far/far_s41*_census.root");
   TGraph G;
   {
      std::ifstream in("pid/proton_d2_dv1104.json");
      std::string str((std::istreambuf_iterator<char>(in)), {});
      auto p = str.find("vertices"); p = str.find('[', p) + 1;
      while (true) {
         auto a = str.find('[', p); if (a == std::string::npos) break;
         auto b = str.find(']', a); double x, y; sscanf(str.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &x, &y);
         G.AddPoint(x, y); p = b + 1; if (str[str.find_first_not_of(" \n\t,", p)] == ']') break;
      }
   }
   auto cg = new TCutG("gp", G.GetN(), G.GetX(), G.GetY()); cg->SetVarX("ssq"); cg->SetVarY("sbr");
   TString W = "sv&&sth>89&&sth<161&&ske>1.2";
   const int NT = 3, NZ = 9;
   double TB[NT + 1] = {10, 20, 30, 40};
   double ZB[NZ + 1] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 980};
   const double XL = -3, XH = 2; const int NB = 50; const double BW = (XH - XL) / NB;
   auto mkF = [&](const char *n) {
      auto f = new TF1(n, "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", XL, XH);
      f->SetParNames("amp", "mean", "sigma", "c0", "c1"); return f;
   };
   double Y[NZ][NT], EY[NZ][NT], S[NZ][NT], MEAN[NT], SIG[NT];
   auto cs = new TCanvas("cs", "", 1800, 900); cs->Divide(NZ, NT, 0.001, 0.001);
   for (int t = 0; t < NT; ++t) {
      TString ca = Form("&&stcm>=%g&&stcm<%g", TB[t], TB[t + 1]);
      // global fit for mean/sigma
      TH1D hg(Form("hg%d", t), "", NB, XL, XH);
      d.Draw(Form("sex>>hg%d", t), W + "&&gp" + ca + "&&svz>=200&&svz<980", "goff");
      auto fg = mkF(Form("fg%d", t));
      double pk = hg.GetBinContent(hg.FindBin(-0.75));
      fg->SetParameters(pk, -0.75, 0.3, hg.GetBinContent(hg.FindBin(1.5)), 0);
      fg->SetParLimits(1, -1.5, 0.5); fg->SetParLimits(2, 0.05, 1.0); fg->SetParLimits(0, 0, 1e6);
      auto rg = hg.Fit(fg, "RQ0LS");
      MEAN[t] = fg->GetParameter(1); SIG[t] = fg->GetParameter(2);
      printf("\n=== thcm %2.0f-%2.0f  GLOBAL (vz 200-980): status %d  -2lnL %.2f  N=%.0f\n", TB[t], TB[t + 1],
             (int)rg, rg->MinFcnValue(), hg.Integral());
      for (int i = 0; i < 5; ++i) {
         double lo, hi; fg->GetParLimits(i, lo, hi); double v = fg->GetParameter(i); bool rail = (lo < hi) && (fabs(v - lo) < 1e-4 * (1 + fabs(v)) || fabs(v - hi) < 1e-4 * (1 + fabs(v)));   // judge on the VALUE
         printf("   %-6s %9.4f +- %8.4f %s\n", fg->GetParName(i), fg->GetParameter(i), fg->GetParError(i), rail ? "  <<< RAILED" : "");
      }
      printf("  vz slice  |  N   amp      c0      c1     -2lnL  | Yfit +- err | free mean  free sig | SIM win\n");
      for (int z = 0; z < NZ; ++z) {
         TString cz = Form("&&svz>=%g&&svz<%g", ZB[z], ZB[z + 1]);
         auto h = new TH1D(Form("h%d_%d", t, z), Form("#theta_{cm} %g-%g  vz %g-%g", TB[t], TB[t + 1], ZB[z], ZB[z + 1]), NB, XL, XH);
         d.Draw(Form("sex>>h%d_%d", t, z), W + "&&gp" + ca + cz, "goff");
         S[z][t] = s.GetEntries(W + ca + cz + "&&sex>-1.5&&sex<0.729");
         auto f = mkF(Form("f%d_%d", t, z));
         f->SetParameters(std::max(1.0, h->GetBinContent(h->FindBin(MEAN[t]))), MEAN[t], SIG[t], std::max(0.1, h->Integral() / NB), 0);
         f->SetParLimits(0, 0, 1e5);           // limits FIRST, then fix (SetParLimits after Fix releases)
         f->FixParameter(1, MEAN[t]); f->FixParameter(2, SIG[t]);
         auto r = h->Fit(f, "RQ0LS");
         double yld = f->GetParameter(0) * SIG[t] * sqrt(2 * TMath::Pi()) / BW;
         double ey = f->GetParError(0) * SIG[t] * sqrt(2 * TMath::Pi()) / BW;
         ey = std::max(ey, sqrt(std::max(yld, 1.0)));   // floor at sqrt(N)
         Y[z][t] = yld; EY[z][t] = ey;
         // free diagnostic
         auto ff = mkF(Form("ff%d_%d", t, z));
         ff->SetParameters(f->GetParameter(0) + 1, MEAN[t], SIG[t], f->GetParameter(3), f->GetParameter(4));
         ff->SetParLimits(0, 0, 1e5); ff->SetParLimits(1, -1.5, 0.5); ff->SetParLimits(2, 0.05, 1.0);
         h->Fit(ff, "RQ0L");
         bool frail = fabs(ff->GetParameter(2) - 0.05) < 1e-3 || fabs(ff->GetParameter(2) - 1.0) < 1e-3 || fabs(ff->GetParameter(1) + 1.5) < 1e-3 || fabs(ff->GetParameter(1) - 0.5) < 1e-3;
         printf("  %4.0f-%4.0f | %4.0f %6.2f %7.3f %7.3f %8.2f | %6.1f +- %4.1f | %8.3f %8.3f%s | %5.0f\n", ZB[z], ZB[z + 1],
                h->Integral(), f->GetParameter(0), f->GetParameter(3), f->GetParameter(4), r->MinFcnValue(), yld, ey,
                ff->GetParameter(1), ff->GetParameter(2), frail ? "R" : " ", S[z][t]);
         cs->cd(t * NZ + z + 1); h->SetLineColor(kBlack); h->Draw("E"); f->SetLineColor(kRed); f->Draw("same");
      }
   }
   cs->SaveAs("fit_vz_gs_dp_slices.png");
   double dn = 0, sn = 0;
   for (int z = 1; z < 8; ++z) { dn += Y[z][2]; sn += S[z][2]; }
   double norm = dn / sn;
   printf("\nnormalisation (fitted data / sim, thcm 30-40, vz 200-900) = %.4f\n", norm);
   printf("\ndata / (sim x norm), fitted g.s. yield\n  vz slice  |     10-20          20-30          30-40\n");
   auto c = new TCanvas("c", "", 900, 600);
   auto mg = new TMultiGraph(); auto leg = new TLegend(0.15, 0.65, 0.45, 0.88);
   int col[NT] = {kRed + 1, kBlue + 1, kGreen + 2};
   for (int t = 0; t < NT; ++t) {
      auto g = new TGraphErrors();
      for (int z = 0; z < NZ; ++z) if (S[z][t] > 0) {
         double r = Y[z][t] / (S[z][t] * norm), er = r * sqrt(pow(EY[z][t] / std::max(Y[z][t], 1.0), 2) + 1.0 / S[z][t]);
         if (Y[z][t] < 1) er = EY[z][t] / (S[z][t] * norm);
         g->AddPoint(0.5 * (ZB[z] + ZB[z + 1]), r); g->SetPointError(g->GetN() - 1, 0.5 * (ZB[z + 1] - ZB[z]), er);
      }
      g->SetMarkerStyle(20); g->SetMarkerColor(col[t]); g->SetLineColor(col[t]); mg->Add(g, "P");
      leg->AddEntry(g, Form("#theta_{cm} %g-%g", TB[t], TB[t + 1]), "lp");
   }
   for (int z = 0; z < NZ; ++z) {
      printf("  %4.0f-%4.0f |", ZB[z], ZB[z + 1]);
      for (int t = 0; t < NT; ++t) {
         if (S[z][t] <= 0) { printf("       -       "); continue; }
         double r = Y[z][t] / (S[z][t] * norm), er = EY[z][t] / (S[z][t] * norm);
         printf("  %5.2f +- %4.2f ", r, er);
      }
      printf("\n");
   }
   mg->SetTitle("16C(d,p) g.s.: fitted data yield / (far sim #times norm);vertex z (mm);data / sim");
   mg->Draw("A"); mg->GetYaxis()->SetRangeUser(0, std::max(2.5, mg->GetYaxis()->GetXmax()));
   TLine l(100, 1, 980, 1); l.SetLineStyle(2); l.Draw(); leg->Draw();
   c->SaveAs("fit_vz_gs_dp.png");
}
