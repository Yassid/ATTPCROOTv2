/// @file plot_pid_spyral_veto_dp.C
/// @brief Spyral PID (run + pad gain sqrt(dE/dx) sqp_g vs Spyral Brho br) with the FAKE-TRITON VETO (2026-11-01):
///        a track is dropped when Spyral's Brho exceeds the 3D circle Brho by more than vetoLn in log,
///        ln(br / bA) > vetoLn (default 0.25 = 28 %). Slow elastic deuterons whose Spyral Brho is ~1.5x too high
///        (Spyral's polar angle too steep) form the blob at the low-Brho end of the deuteron band; in the d + t sims
///        the cut keeps 97.5-100 % of tritons above circle Brho 0.25 and removes the sim deuterons' own plume.
///        Top: plane before, after, and the removed tracks. Bottom: ln-mass per Spyral-Brho slice, before (grey) and
///        after (blue), each shifted so its deuteron mode is 0; lines at p, t, 4He.
///        vetoMaxBrho / vetoHighSq: apply the veto only below this Spyral Brho, or above it only for sqrt(dE/dx) >
///        vetoHighSq (adopted 2026-11-01: 0.45 and 30 -- the blob plus the upper streak, normal bands untouched).
///   root -l -b -q 'plot_pid_spyral_veto_dp.C("plots/pid_spyral_veto.png",0.25,0.45,30)'
void plot_pid_spyral_veto_dp(TString png, double vetoLn = 0.25, double vetoMaxBrho = 99, double vetoHighSq = 1e9, TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root",
                             TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&sqp_g>0&&br>0&&bA>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   // the veto acts where Spyral's Brho is below vetoMaxBrho (the blob region) OR, above it, where sqrt(dE/dx) exceeds
   // vetoHighSq (the upper part of the fake-triton streak, which continues to Brho ~0.7 at high dE/dx)
   TString drop = Form("(log(br/bA)>=%g&&(br<%g||sqp_g>%g))", vetoLn, vetoMaxBrho, vetoHighSq), keep = "!" + drop;
   w->Draw("sqp_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2], xhi = x[x.size() * 995 / 1000] * 1.1;
   printf("veto removes %lld of %lld tracks\n", w->GetEntries(sel + "&&" + drop), w->GetEntries(sel));
   auto *c = new TCanvas("c", "c", 2000, 1300);
   auto *top = new TPad("top", "", 0, 0.45, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.45); bot->Draw();
   top->Divide(3, 1);
   TString cuts[3] = {sel, sel + "&&" + keep, sel + "&&" + drop};
   const char *ttl[3] = {"Spyral, run + pad gain: all", vetoMaxBrho < 50 ? Form("veto ln(B#rho_{S}/B#rho_{circle}) > %.2f where B#rho_{S} < %.2f or #sqrt{dE/dx} > %.0f", vetoLn, vetoMaxBrho, vetoHighSq) : Form("after the veto ln(B#rho_{Spyral}/B#rho_{circle}) < %.2f", vetoLn), "removed by the veto"};
   for (int p = 0; p < 3; ++p) {
      top->cd(p + 1); gPad->SetLogz(); gPad->SetRightMargin(0.12);
      auto *h = new TH2D(Form("p%d", p), Form("%s;#sqrt{dE/dx} (Spyral, run+pad gain);B#rho_{Spyral} [T m]", ttl[p]), 300, 0, xhi, 300, 0, 1.0);
      w->Draw(Form("br:sqp_g>>p%d", p), cuts[p], "goff");
      h->Draw("colz");
   }
   const double bl[] = {0.15, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80};
   bot->Divide(6, 1);
   TString L = Form("log(sqp_g/%g*br)", md);
   for (int b = 0; b < 6; ++b) {
      bot->cd(b + 1); gPad->SetLogy(); gPad->SetGridx();
      TString sl = Form("&&br>%g&&br<%g", bl[b], bl[b + 1]);
      TH1D tmp(Form("t%d", b), "", 600, -4, 2);
      w->Draw(L + Form(">>t%d", b), sel + sl, "goff");
      tmp.Smooth(2);
      double d = tmp.GetBinCenter(tmp.GetMaximumBin());
      TH1D *h[2];
      for (int k = 0; k < 2; ++k) {
         h[k] = new TH1D(Form("h%d%d", b, k), Form("B#rho_{Spyral} %.2f-%.2f T m;ln mass - ln m_{d};tracks", bl[b], bl[b + 1]), 160, -1.2, 1.2);
         w->Draw(Form("%s-(%g)>>h%d%d", L.Data(), d, b, k), (k ? sel + "&&" + keep : sel) + sl, "goff");
         h[k]->SetLineWidth(2); h[k]->SetLineColor(k ? kBlue + 2 : kGray + 1);
      }
      h[0]->SetMinimum(0.8); h[0]->SetMaximum(3 * h[0]->GetMaximum());
      h[0]->Draw("hist"); h[1]->Draw("hist same");
      for (double v : {std::log(938.272 / 1875.613), std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)}) {
         auto *l = new TLine(v, 0.8, v, h[0]->GetMaximum()); l->SetLineStyle(3); l->Draw();
      }
      if (b == 0) {
         auto *lg = new TLegend(0.12, 0.78, 0.7, 0.89);
         lg->AddEntry(h[0], "all", "l"); lg->AddEntry(h[1], "after veto", "l"); lg->Draw();
      }
   }
   c->SaveAs(png);
}
