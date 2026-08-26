/// @file omp_verdict_C14.C
/// @brief Everything that bears on WHICH optical potential describes p+14C, in one figure.
///
/// Three tests, and they do not all agree, which is the point:
///   (a) the diffraction MINIMUM. Measured at theta_cm 58.0 deg (elastic_dip_C14.C). Becchetti-
///       Greenlees and Menet put it there exactly; KD03 is 5 deg out.
///   (b) the SHAPE away from the minimum, over 20-50 and 70-148 deg. Perey best, KD03 second --
///       this is the ONE test that favours the potential the analysis actually uses.
///   (c) the LUMINOSITY against the FRIBDAQ scalers. This is the only test with no optical model
///       in it at all: L = ic_sca * f_beam * livetime * n_target. CH89, Menet and Becchetti-
///       Greenlees agree with the hardware beam count to a few %; KD03 is ~28% low.
///
/// Inputs, all measured and all on disk:
///   plots/omp_luminosity.txt     <- elastic_dip_C14.C
///   plots/scaler_luminosity.txt  <- scaler_lumi_C14.C
///   plots/ic_fraction.txt        <- ic_fraction_C14.C
///
///   root -b -q 'omp_verdict_C14.C()'
void omp_verdict_C14()
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   std::vector<TString> nm; std::vector<double> dip, L, rms;
   { std::ifstream in((here + "/plots/omp_luminosity.txt").Data());
     if (!in) { printf("\033[1;31mrun elastic_dip_C14.C first\033[0m\n"); return; }
     std::string line;
     while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream is(line); std::string k, w; is >> k;
        std::vector<std::string> t; while (is >> w) t.push_back(w);
        if (t.size() < 3) continue;
        std::string n; for (size_t i = 0; i + 3 < t.size(); ++i) n += (i ? " " : "") + t[i];
        nm.push_back(n.c_str()); dip.push_back(std::atof(t[t.size()-3].c_str()));
        L.push_back(std::atof(t[t.size()-2].c_str())); rms.push_back(std::atof(t[t.size()-1].c_str()));
     } }
   const int NP = nm.size();
   double Lc = 0, Lt = 0, fB = 0, icsca = 0;
   { std::ifstream in((here + "/plots/scaler_luminosity.txt").Data());
     if (!in) { printf("\033[1;31mrun scaler_lumi_C14.C first\033[0m\n"); return; }
     std::string l1, l2, l3; std::getline(in, l1); std::getline(in, l2); std::getline(in, l3);
     double lc, lt, nt, Ln; in >> icsca >> lc >> lt >> nt >> Lc >> Lt >> Ln;
     TString p(l2.c_str()); int i = p.Index("fBeam="); if (i >= 0) fB = TString(p(i + 6, 12)).Atof(); }
   const double thDip = 57.95;   // measured, elastic_dip_C14.C

   printf("\n  scaler L: clock %.1f, trigger %.1f mb^-1  (f_beam %.3f, ic_sca %.4g)\n", Lc, Lt, fB, icsca);
   printf("\n  %-24s %6s %8s %8s %10s\n", "potential", "dip", "L", "L/Lsc", "rms");
   for (int i = 0; i < NP; ++i)
      printf("  %-24s %6.1f %8.1f %8.3f %10.3f\n", nm[i].Data(), dip[i], L[i], L[i] / Lc, rms[i]);

   int col[8] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};
   auto *c = new TCanvas("cov", "", 1500, 1000); c->Divide(2, 2);

   // (a) dip position
   c->cd(1); gPad->SetGridy(); gPad->SetBottomMargin(0.22);
   auto *h1 = new TH1D("h1", "(a) diffraction minimum;;#theta_{dip} [deg]", NP, 0, NP);
   for (int i = 0; i < NP; ++i) { h1->SetBinContent(i + 1, dip[i]); h1->GetXaxis()->SetBinLabel(i + 1, nm[i]); }
   h1->SetFillColor(kAzure - 4); h1->SetBarWidth(0.6); h1->SetBarOffset(0.2);
   h1->GetXaxis()->SetLabelSize(0.052); h1->GetXaxis()->LabelsOption("v");
   h1->SetMinimum(54); h1->SetMaximum(65); h1->Draw("bar");
   auto *lm = new TLine(0, thDip, NP, thDip); lm->SetLineColor(kRed + 1); lm->SetLineWidth(3); lm->Draw();
   TLatex t1; t1.SetNDC(); t1.SetTextSize(0.048); t1.SetTextColor(kRed + 1);
   t1.DrawLatex(0.45, 0.30, Form("measured %.1f#circ", thDip));

   // (b) shape away from the dip
   c->cd(2); gPad->SetGridy(); gPad->SetBottomMargin(0.22);
   auto *h2 = new TH1D("h2", "(b) shape away from the minimum;;rms ln(data/calc)", NP, 0, NP);
   for (int i = 0; i < NP; ++i) { h2->SetBinContent(i + 1, rms[i]); h2->GetXaxis()->SetBinLabel(i + 1, nm[i]); }
   h2->SetFillColor(kGreen - 6); h2->SetBarWidth(0.6); h2->SetBarOffset(0.2);
   h2->GetXaxis()->SetLabelSize(0.052); h2->GetXaxis()->LabelsOption("v");
   h2->SetMinimum(0); h2->SetMaximum(0.6); h2->Draw("bar");
   TLatex t2; t2.SetNDC(); t2.SetTextSize(0.044);
   t2.DrawLatex(0.32, 0.86, "lower is better -- the only test");
   t2.DrawLatex(0.32, 0.80, "that favours KD03");

   // (c) luminosity vs the scalers
   c->cd(3); gPad->SetGridy(); gPad->SetBottomMargin(0.22);
   auto *h3 = new TH1D("h3", "(c) luminosity vs the hardware beam count;;L [counts/mb]", NP, 0, NP);
   for (int i = 0; i < NP; ++i) { h3->SetBinContent(i + 1, L[i]); h3->GetXaxis()->SetBinLabel(i + 1, nm[i]); }
   h3->SetFillColor(kOrange - 3); h3->SetBarWidth(0.6); h3->SetBarOffset(0.2);
   h3->GetXaxis()->SetLabelSize(0.052); h3->GetXaxis()->LabelsOption("v");
   h3->SetMinimum(0); h3->SetMaximum(1.5 * std::max(Lc, *std::max_element(L.begin(), L.end())));
   h3->Draw("bar");
   // the scaler band: clock to trigger livetime
   auto *bx = new TBox(0, std::min(Lc, Lt), NP, std::max(Lc, Lt));
   bx->SetFillColorAlpha(kRed + 1, 0.20); bx->SetLineColor(kRed + 1); bx->Draw("l");
   auto *lc2 = new TLine(0, Lc, NP, Lc); lc2->SetLineColor(kRed + 1); lc2->SetLineWidth(3); lc2->Draw();
   TLatex t3; t3.SetNDC(); t3.SetTextSize(0.044); t3.SetTextColor(kRed + 1);
   t3.DrawLatex(0.13, 0.86, Form("scalers: %.1f - %.1f mb^{-1}", std::min(Lc, Lt), std::max(Lc, Lt)));
   t3.SetTextColor(kBlack);
   t3.DrawLatex(0.13, 0.80, Form("f_{beam} = %.3f measured, no optical model", fB));

   // (d) the two independent tests against each other
   c->cd(4); gPad->SetGridx(); gPad->SetGridy();
   auto *fr = gPad->DrawFrame(-1.5, 0.55, 6.5, 1.35);
   fr->SetTitle("(d) the two model-free tests;#theta_{dip}^{calc} - #theta_{dip}^{data} [deg];L_{elastic} / L_{scaler}");
   for (int i = 0; i < NP; ++i) {
      auto *g = new TGraph(1); g->SetPoint(0, dip[i] - thDip, L[i] / Lc);
      g->SetMarkerStyle(20); g->SetMarkerSize(2.2); g->SetMarkerColor(col[i]); g->Draw("P same");
      TLatex tl; tl.SetTextSize(0.036); tl.SetTextColor(col[i]);
      tl.DrawLatex(dip[i] - thDip + 0.18, L[i] / Lc + 0.015, nm[i]);
   }
   auto *u = new TLine(-1.5, 1, 6.5, 1); u->SetLineColor(kRed + 1); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
   auto *v = new TLine(0, 0.55, 0, 1.35); v->SetLineColor(kRed + 1); v->SetLineStyle(2); v->SetLineWidth(2); v->Draw();
   TLatex t4; t4.SetTextSize(0.038); t4.SetTextColor(kRed + 1);
   t4.DrawLatex(-1.3, 1.28, "target: dip right AND luminosity right");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "09_omp_verdict.png");
   printf("\n  wrote %s09_omp_verdict.png\n\n", out.Data());
}
