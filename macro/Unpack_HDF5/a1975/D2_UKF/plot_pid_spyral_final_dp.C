/// @file plot_pid_spyral_final_dp.C
/// @brief Spyral PID with the fake-triton veto AND the vetoed tracks added back from the circle (2026-11-01).
///        Kept tracks: Spyral (sqp_g, br). Vetoed tracks (ln(br/bA) > vetoLn where br < vetoMaxBrho or sqp_g > vetoHighSq):
///        slow elastic deuterons near 90 deg whose Spyral radius is ~1.6x too large (flat in z, Spyral's z-ordered
///        spline is ill-posed); they are re-added at the circle Brho bA with the 3D dE/dx scaled to the Spyral scale,
///        s1_g * scale, scale = median(sqp_g / s1_g) on kept tracks with br 0.12-0.25 and |ln(br/bA)| < 0.1.
///        Top: final plane, and the re-added tracks alone. Bottom: ln-mass per slice of the plotted Brho, final (blue),
///        kept only (grey), re-added only (magenta); the deuteron mode of the final sample in each slice is 0.
///   root -l -b -q 'plot_pid_spyral_final_dp.C("plots/pid_spyral_final.png")'
void plot_pid_spyral_final_dp(TString png, double vetoLn = 0.25, double vetoMaxBrho = 0.45, double vetoHighSq = 30,
                              TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root",
                              TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&sqp_g>0&&br>0&&bA>0&&s1_g>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   TString drop = Form("(log(br/bA)>=%g&&(br<%g||sqp_g>%g))", vetoLn, vetoMaxBrho, vetoHighSq);
   w->Draw("sqp_g/s1_g", sel + "&&!" + drop + "&&br>0.12&&br<0.25&&abs(log(br/bA))<0.1", "goff");
   std::vector<double> r(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(r.begin(), r.end());
   const double sc = r[r.size() / 2];
   // the plotted coordinates: Spyral for kept tracks, circle Brho + scaled 3D dE/dx for re-added ones
   TString B = Form("(%s?bA:br)", drop.Data()), SQ = Form("(%s?s1_g*%g:sqp_g)", drop.Data(), sc);
   w->Draw("sqp_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2], xhi = x[x.size() * 995 / 1000] * 1.1;
   printf("scale %.4f; re-added %lld of %lld tracks\n", sc, w->GetEntries(sel + "&&" + drop), w->GetEntries(sel));
   auto *c = new TCanvas("c", "c", 2000, 1300);
   auto *top = new TPad("top", "", 0, 0.45, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.45); bot->Draw();
   top->Divide(2, 1);
   const char *ttl[2] = {"Spyral PID + veto, vetoed deuterons re-added from the circle", "re-added tracks (circle B#rho, scaled 3D dE/dx)"};
   TString cuts[2] = {sel, sel + "&&" + drop};
   for (int p = 0; p < 2; ++p) {
      top->cd(p + 1); gPad->SetLogz(); gPad->SetRightMargin(0.12);
      auto *h = new TH2D(Form("p%d", p), Form("%s;#sqrt{dE/dx} (Spyral scale);B#rho [T m]", ttl[p]), 300, 0, xhi, 300, 0, 1.0);
      w->Draw(Form("%s:%s>>p%d", B.Data(), SQ.Data(), p), cuts[p], "goff");
      h->Draw("colz");
   }
   const double bl[] = {0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.80};
   bot->Divide(6, 1);
   TString L = Form("log(%s/%g*%s)", SQ.Data(), md, B.Data());
   for (int b = 0; b < 6; ++b) {
      bot->cd(b + 1); gPad->SetLogy(); gPad->SetGridx();
      TString sl = Form("&&%s>%g&&%s<%g", B.Data(), bl[b], B.Data(), bl[b + 1]);
      // deuteron mode of the final sample: the highest smoothed bin within +-0.3 of the mode in the 0.15-0.25 slice
      static double dref = 0;
      if (b == 0) {
         TH1D hr("hr", "", 600, -4, 2);
         w->Draw(L + ">>hr", sel + Form("&&%s>0.15&&%s<0.25", B.Data(), B.Data()), "goff");
         hr.Smooth(2); dref = hr.GetBinCenter(hr.GetMaximumBin());
      }
      TH1D tmp(Form("t%d", b), "", 600, -4, 2);
      w->Draw(L + Form(">>t%d", b), sel + sl, "goff");
      tmp.Smooth(2);
      int i0 = tmp.FindBin(dref - 0.3), i1 = tmp.FindBin(dref + 0.3), bd = i0;
      for (int i = i0; i <= i1; ++i) if (tmp.GetBinContent(i) > tmp.GetBinContent(bd)) bd = i;
      const double d = tmp.GetBinCenter(bd);
      TString cs[3] = {sel, sel + "&&!" + drop, sel + "&&" + drop};
      const int col[3] = {kBlue + 2, kGray + 1, kMagenta + 1};
      TH1D *h[3];
      for (int k = 0; k < 3; ++k) {
         h[k] = new TH1D(Form("h%d%d", b, k), Form("B#rho %.2f-%.2f T m;ln mass - ln m_{d};tracks", bl[b], bl[b + 1]), 160, -1.2, 1.2);
         w->Draw(Form("%s-(%g)>>h%d%d", L.Data(), d, b, k), cs[k] + sl, "goff");
         h[k]->SetLineWidth(2); h[k]->SetLineColor(col[k]);
      }
      h[0]->SetMinimum(0.8); h[0]->SetMaximum(3 * h[0]->GetMaximum());
      h[0]->Draw("hist"); h[1]->Draw("hist same"); h[2]->Draw("hist same");
      for (double v : {std::log(938.272 / 1875.613), std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)}) {
         auto *l = new TLine(v, 0.8, v, h[0]->GetMaximum()); l->SetLineStyle(3); l->Draw();
      }
      if (b == 0) {
         auto *lg = new TLegend(0.12, 0.74, 0.75, 0.89);
         lg->AddEntry(h[0], "final", "l"); lg->AddEntry(h[1], "Spyral, kept", "l"); lg->AddEntry(h[2], "re-added (circle)", "l"); lg->Draw();
      }
   }
   c->SaveAs(png);
}
