#include "ddbw.C"
// The (d,d') spectrum with the 15C(p,p') six-resonance scheme: barrier-corrected Breit-Wigners
// with Lane-Thomas penetrability, convolved with the detector resolution, on a phase-space
// continuum.  Everything is evaluated into TGraphs immediately -- TF1::Draw is lazy and would
// otherwise paint panel 1 with panel 2's globals (the templates TB[] are rebuilt per fit).
static TGraph* gcurve(std::function<double(double)> f, double lo, double hi, int N = 1600)
{
   auto* g = new TGraph(N);
   for (int i = 0; i < N; i++) { double x = lo + (hi - lo) * i / (N - 1.0); g->SetPoint(i, x, f(x)); }
   return g;
}
void panelBW(const char* nm, double fw, double lo)
{
   Res r = fit(nm, fw, lo, 0.0, true, false);
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   double s0 = fw / 2.355;

   // rebuild the fitted pieces as explicit lambdas over the CURRENT templates
   double gsA = 0, A[NBW];
   {  // recover amplitudes: r.amp[i] is the peak height, convert back to the template scale
      for (int i = 0; i < NBW; i++) A[i] = r.amp[i] / std::max(1e-12, TB[i].Eval(r.mu[i]));
   }
   // the g.s. amplitude is not returned by Res; refit quickly for it
   TF1 gq("gq2", "gaus(0)+pol1(3)", lo, 1.6);
   gq.SetParameters(h->GetBinContent(h->GetMaximumBin()), r.gsMu, s0, 10, 0);
   gq.SetParLimits(2, 0.05, 1.5); h->Fit(&gq, "RQ0N"); gsA = gq.GetParameter(0);

   auto cont = [&](double X) { return r.fl + r.ps * (HPS ? HPS->Interpolate(X) : 0.0); };
   auto gsf  = [&](double X) { double s = (X < r.gsMu) ? r.fwL / 2.355 : r.fwR / 2.355;
                               return gsA * std::exp(-0.5 * std::pow((X - r.gsMu) / s, 2)); };
   auto e5f  = [&](double X) { return r.r52 * gsA * std::exp(-0.5 * std::pow((X - r.gsMu - E52) / (r.fwR / 2.355), 2)); };
   auto tot  = [&](double X) { double v = cont(X) + gsf(X) + e5f(X);
                               for (int i = 0; i < NBW; i++) v += A[i] * TB[i].Eval(X); return v; };

   h->SetTitle(Form("%s  --  ^{15}C(d,d') with the ^{15}C(p,p') six-resonance scheme"
                    ";E_{x} (MeV);counts / 100 keV", nm));
   h->SetMarkerStyle(20); h->SetMarkerSize(0.45); h->SetLineColor(kBlack);
   h->GetXaxis()->SetRangeUser(lo, 14.0);
   h->SetMinimum(25); h->SetMaximum(h->GetBinContent(h->GetMaximumBin()) * 3.5);
   h->Draw("E");

   TGraph* gc = gcurve(cont, lo, 14); gc->SetLineColor(kGray + 2); gc->SetLineStyle(2); gc->SetLineWidth(2); gc->Draw("L");
   TGraph* gg = gcurve(gsf, lo, 3.2);  gg->SetLineColor(kAzure + 2); gg->SetLineStyle(7); gg->SetLineWidth(2); gg->Draw("L");
   TGraph* ge = gcurve(e5f, lo, 3.6);  ge->SetLineColor(kViolet + 1); ge->SetLineWidth(3); ge->Draw("L");
   int col[NBW] = {kGreen + 2, kSpring - 6, kCyan + 2, kOrange + 7, kRed - 4, kMagenta + 2};
   for (int i = 0; i < NBW; i++) {
      TGraph* gb = gcurve([&, i](double X) { return cont(X) + A[i] * TB[i].Eval(X); }, 2.0, 14.0);
      gb->SetLineColor(col[i]); gb->SetLineStyle(7); gb->SetLineWidth(2); gb->Draw("L");
   }
   TGraph* gt = gcurve(tot, lo, 14); gt->SetLineColor(kRed); gt->SetLineWidth(3); gt->Draw("L");
   h->Draw("E same");

   TLine* ln = new TLine(Bar::SN + r.gsMu, 25, Bar::SN + r.gsMu, h->GetMaximum() * 0.2);
   ln->SetLineStyle(3); ln->Draw();
   TLatex sn; sn.SetTextSize(0.028); sn.DrawLatex(Bar::SN + r.gsMu + 0.1, h->GetMaximum() * 0.22, "S_{n}");

   TLatex t; t.SetNDC(true); t.SetTextSize(0.027);
   t.DrawLatex(0.47, 0.90, Form("#chi^{2}/ndf = %.2f   (energies and #Gamma FIXED to (p,p'))", r.chi2ndf));
   t.SetTextColor(kGray + 2);

   for (int i = 0; i < NBW; i++) {
      t.SetTextColor(col[i]);
      t.DrawLatex(0.47, 0.862 - i * 0.034,
                  Form("BW%d  #font[12]{l}=%d  E_{x}=%.2f   %.0f#sigma", i + 1, LPP[i], r.mu[i],
                       r.err[i] > 0 ? r.amp[i] / r.err[i] : 0.0));
   }
   t.SetTextColor(kGray + 2);
   t.DrawLatex(0.47, 0.862 - NBW * 0.034, Form("continuum: PS#times%.3g + floor %.0f", r.ps, r.fl));
}
void plotbw()
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.040, "t");
   TCanvas* c = new TCanvas("cbw", "", 1900, 760);
   c->Divide(2, 1);
   c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03);
   panelBW("ATTPCROOT", 1.012, -1.6);
   c->cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03);
   panelBW("Spyral", 0.548, -1.2);
   c->SaveAs("/home/yassid/Data/C15d_statefit/dd_bw_ppscheme.png");
   printf("wrote dd_bw_ppscheme.png\n");
}
