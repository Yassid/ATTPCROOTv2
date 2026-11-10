#include "degen.C"
// What limits the 5/2+ amplitude ratio: statistics, or the assumed g.s. lineshape?
// Scan the FIXED g.s. right-hand width, refit everything else, and plot the ratio against it
// together with Delta(-2lnL). The likelihood is steep -- so the width is pinned and the ratio is
// NOT internally degenerate -- but the ratio slides fast with the width, so a lineshape error of
// a few per cent moves it by tens of per cent.
void runScan(const char* nm, double kb, double fw, double lo, int mode,
             std::vector<double>& W, std::vector<double>& R, std::vector<double>& D, Out& o)
{
   o = one(nm, kb, fw, lo, mode, false);
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   double best = 1e30; std::vector<double> L;
   for (double fwR = 0.40; fwR <= 1.251; fwR += 0.025) {
      TF1 g(Form("fg_%s_%d", nm, (int)(fwR * 1000)), fpk2, lo, SPLIT, 5);
      g.SetParameter(0, h->GetMaximum()); g.SetParLimits(0, 0, 1e9);
      g.SetParameter(1, o.mu);            g.SetParLimits(1, -0.8, 0.8);
      g.SetParameter(2, o.fwL / 2.355);   g.SetParLimits(2, 0.05, 1.2);
      g.FixParameter(3, fwR / 2.355);
      g.SetParameter(4, o.r52);           g.SetParLimits(4, 0, 0.9);
      TFitResultPtr r = h->Fit(&g, "RQ0NLS");
      double v = 2.0 * r->MinFcnValue();
      W.push_back(fwR); R.push_back(g.GetParameter(4)); L.push_back(v);
      best = std::min(best, v);
   }
   for (double v : L) D.push_back(v - best);
}
void panel(const char* nm, double kb, double fw, double lo, int mode, double other, const char* othernm)
{
   std::vector<double> W, R, D; Out o;
   runScan(nm, kb, fw, lo, mode, W, R, D, o);
   int n = W.size();
   TGraph* gr = new TGraph(n, W.data(), R.data());
   gr->SetTitle(Form("%s  --  5/2^{+} ratio vs the ASSUMED g.s. right-hand width;"
                     "g.s. FWHM_{R} (MeV, fixed);5/2^{+} / g.s. amplitude ratio", nm));
   gr->SetLineColor(kViolet + 1); gr->SetLineWidth(3);
   gr->GetYaxis()->SetRangeUser(0, 0.20);
   gr->GetXaxis()->SetLimits(0.40, 1.25);
   gr->Draw("AL");

   // the +-1 sigma band on the width, from Delta(-2lnL) = 1, INTERPOLATED between grid points
   // (the likelihood is so steep that only one 25 keV grid point has Delta <= 1, which would
   //  otherwise print a zero-width band)
   int ib = 0; for (int i = 1; i < n; i++) if (D[i] < D[ib]) ib = i;
   double wlo = W[0], whi = W[n - 1];
   for (int i = ib; i > 0; i--)
      if (D[i - 1] >= 1.0 && D[i] < 1.0) { wlo = W[i - 1] + (W[i] - W[i - 1]) * (D[i - 1] - 1.0) / (D[i - 1] - D[i]); break; }
   for (int i = ib; i < n - 1; i++)
      if (D[i] < 1.0 && D[i + 1] >= 1.0) { whi = W[i] + (W[i + 1] - W[i]) * (1.0 - D[i]) / (D[i + 1] - D[i]); break; }
   TBox* b = new TBox(wlo, 0, whi, 0.20);
   b->SetFillColorAlpha(kAzure + 1, 0.25); b->SetLineColor(kAzure + 2); b->Draw("l");
   gr->Draw("L same");

   TMarker* m = new TMarker(o.fwR, o.r52, 20);
   m->SetMarkerColor(kRed); m->SetMarkerSize(1.7); m->Draw();
   TLine* lo2 = new TLine(0.40, other, 1.25, other);
   lo2->SetLineColor(kOrange + 7); lo2->SetLineStyle(2); lo2->SetLineWidth(3); lo2->Draw();

   // where would this code have to sit to reproduce the other code's ratio?
   double wneed = -1;
   for (int i = 1; i < n; i++)
      if ((R[i - 1] - other) * (R[i] - other) <= 0 && wneed < 0)
         wneed = W[i - 1] + (W[i] - W[i - 1]) * (other - R[i - 1]) / (R[i] - R[i - 1] + 1e-12);
   double s10 = 0;
   for (int i = 1; i < n; i++) if (W[i - 1] <= o.fwR && W[i] > o.fwR) s10 = (R[i] - R[i - 1]) / (W[i] - W[i - 1]);
   TLatex t; t.SetNDC(true); t.SetTextSize(0.031);
   t.SetTextColor(kRed);
   t.DrawLatex(0.17, 0.355, Form("best fit  FWHM_{R} = %.3f,  ratio = %.4f #pm %.4f (stat)", o.fwR, o.r52, o.dr52));
   t.SetTextColor(kAzure + 2);
   t.DrawLatex(0.17, 0.305, Form("#Delta(-2lnL)=1 band  %.3f #minus %.3f  #rightarrow width PINNED to #pm%.0f keV",
                                 wlo, whi, 1000 * 0.5 * (whi - wlo)));
   t.SetTextColor(kOrange + 7);
   if (wneed > 0)
      t.DrawLatex(0.17, 0.255, Form("%s's ratio (%.4f) needs FWHM_{R} = %.2f, i.e. %+.0f %%",
                                    othernm, other, wneed, 100.0 * (wneed - o.fwR) / o.fwR));
   else
      t.DrawLatex(0.17, 0.255, Form("%s's ratio (%.4f) is off this scale", othernm, other));
   t.SetTextColor(kBlack);
   t.DrawLatex(0.17, 0.205, Form("slope %.3f /MeV: a 10 %% width error #rightarrow %.0f %% on the ratio",
                                 s10, std::fabs(100.0 * s10 * 0.1 * o.fwR / o.r52)));
}
void degenfig()
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.042, "t");
   TCanvas* c = new TCanvas("cd", "", 1750, 680);
   c->Divide(2, 1);
   c->cd(1); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetGridy();
   panel("ATTPCROOT", 191.25, 1.012, -1.6, 0, 0.0353, "Spyral");
   c->cd(2); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetGridy();
   panel("Spyral", 195.00, 0.548, -1.2, 1, 0.0997, "ATTPCROOT");
   c->SaveAs("/home/yassid/Data/C15d_statefit/fivehalf_lineshape.png");
   printf("wrote fivehalf_lineshape.png\n");
}
