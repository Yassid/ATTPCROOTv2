#include "ddfinal.C"
// The (d,d') spectrum with the CORRECTED continuum: sequential 15C(d,d')15C*(M), M -> 14C+n,
// folded with the kinematic acceptance, plus the six (p,p') resonances as barrier BWs.
// TGraphs, not TF1s -- the templates are globals rebuilt per fit and TF1::Draw is lazy.
static TGraph* gc(std::function<double(double)> f, double lo, double hi, int N = 1800)
{ auto* g = new TGraph(N);
  for (int i = 0; i < N; i++) { double x = lo + (hi - lo) * i / (N - 1.0); g->SetPoint(i, x, f(x)); }
  return g; }
void pane(const char* nm, double kb, double fw, double lo, double keHi)
{
   R o = fitFinal(nm, kb, fw, lo, 14, 86, (keHi > 19 ? 0.0 : 0.5), keHi, false);
   TH1D* h = o.h;
   double cs = o.cs, FL = o.floor;
   TGraph* ACC = GACC;
   auto cont = [&](double X) { return FL + cs * std::max(0.0, ACC->Eval(X)); };
   auto gsf  = [&](double X) { double s = (X < o.gsMu) ? o.fwL / 2.355 : o.fwR / 2.355;
                               return o.gsAmp * std::exp(-0.5 * std::pow((X - o.gsMu) / s, 2)); };
   auto e5f  = [&](double X) { return o.r52 * o.gsAmp * std::exp(-0.5 * std::pow((X - o.gsMu - E52) / (o.fwR / 2.355), 2)); };
   std::vector<Tmpl> T(NBW); for (int i = 0; i < NBW; i++) T[i] = TB[i];
   auto tot = [&](double X) { double v = cont(X) + gsf(X) + e5f(X);
                              for (int i = 0; i < NBW; i++) v += o.A[i] * T[i].Eval(X); return v; };

   h->SetTitle(Form("%s  --  ^{15}C(d,d'), sequential continuum + (p,p') resonances"
                    ";E_{x} (MeV);counts / 100 keV", nm));
   h->SetMarkerStyle(20); h->SetMarkerSize(0.45); h->SetLineColor(kBlack);
   h->GetXaxis()->SetRangeUser(lo, 20.0);
   h->SetMinimum(25); h->SetMaximum(h->GetBinContent(h->GetMaximumBin()) * 3.5);
   h->Draw("E");
   // the continuum, and the floor alone, so the two pieces are visible
   TGraph* gf = gc([&](double) { return FL; }, lo, 20); gf->SetLineColor(kGray + 1);
   gf->SetLineStyle(3); gf->SetLineWidth(2); gf->Draw("L");
   TGraph* gk = gc(cont, lo, 20); gk->SetLineColor(kGray + 3); gk->SetLineStyle(2); gk->SetLineWidth(3); gk->Draw("L");
   TGraph* gg = gc(gsf, lo, 3.2); gg->SetLineColor(kAzure + 2); gg->SetLineStyle(7); gg->SetLineWidth(2); gg->Draw("L");
   TGraph* ge = gc(e5f, lo, 3.6); ge->SetLineColor(kViolet + 1); ge->SetLineWidth(3); ge->Draw("L");
   int col[NBW] = {kGreen + 2, kSpring - 6, kCyan + 2, kOrange + 7, kRed - 4, kMagenta + 2};
   for (int i = 0; i < NBW; i++) {
      TGraph* gb = gc([&, i](double X) { return cont(X) + o.A[i] * T[i].Eval(X); }, 2.0, 14.0);
      gb->SetLineColor(col[i]); gb->SetLineStyle(7); gb->SetLineWidth(2); gb->Draw("L"); }
   TGraph* gt = gc(tot, lo, 20); gt->SetLineColor(kRed); gt->SetLineWidth(3); gt->Draw("L");
   h->Draw("E same");

   TLine* ln = new TLine(Bar::SN + o.gsMu, 25, Bar::SN + o.gsMu, h->GetMaximum() * 0.18);
   ln->SetLineStyle(3); ln->Draw();
   TLatex sn; sn.SetTextSize(0.027); sn.DrawLatex(Bar::SN + o.gsMu + 0.15, h->GetMaximum() * 0.20, "S_{n}");

   TLatex t; t.SetNDC(true); t.SetTextSize(0.026);
   t.DrawLatex(0.45, 0.905, Form("#chi^{2}/ndf = %.2f", o.chi2ndf));
   for (int i = 0; i < NBW; i++) { t.SetTextColor(col[i]);
      t.DrawLatex(0.45, 0.868 - i * 0.033, Form("BW%d  #font[12]{l}=%d  E_{x}=%.2f   %.0f#sigma",
                  i + 1, LPP[i], o.mu[i], o.err[i] > 0 ? o.amp[i] / o.err[i] : 0.0)); }
   t.SetTextColor(kGray + 3);
   t.DrawLatex(0.45, 0.868 - NBW * 0.033, "continuum = floor + N#timesdens(E_{x})#timesacc(E_{x})");
   t.SetTextColor(kGray + 1);
   t.DrawLatex(0.45, 0.868 - (NBW + 1) * 0.033, Form("floor FIXED at %.0f (measured below S_{n})", FL));
}
void plotfinal()
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.038, "t");
   TCanvas* c = new TCanvas("cf", "", 1900, 780);
   c->Divide(2, 1);
   c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03);
   pane("ATTPCROOT", 191.25, 1.012, -1.6, 18.0);
   c->cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03);
   pane("Spyral", 195.00, 0.548, -1.2, 20.0);
   c->SaveAs("/home/yassid/Data/C15d_statefit/dd_final_seqps.png");
   printf("wrote dd_final_seqps.png\n");
}
