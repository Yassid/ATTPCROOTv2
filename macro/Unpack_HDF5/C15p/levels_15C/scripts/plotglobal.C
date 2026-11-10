#include "globalfit.C"
// One figure, both codes, the complete excitation spectrum with every component drawn.
//
// TRAP (cost me one wrong figure): TF1::Draw is LAZY -- the function is evaluated when the canvas
// is painted, i.e. at SaveAs, not when Draw is called. bkgOf() branches on the global MODE, so the
// ATTPCROOT panel was being painted with Spyral's background mode and parameters. Every draw-only
// function below therefore carries the mode, Kb and sigma_res as EXPLICIT TF1 parameters; nothing
// here reads a global.
//
// Top row:    full range, log y -- the g.s. is 9000 counts/bin and the states ~300, so a linear
//             axis cannot show both.
// Bottom row: linear zoom on the g.s. and the 0.740 MeV 5/2+, which is where those two live.

//  pd[0] mode   pd[1] a   pd[2] b   pd[3] c   pd[4] Kb   pd[5] sigma_res
static double bkgX(double X, const double* pd)
{
   int m = (int)(pd[0] + 0.5);
   if (m == 0) return std::max(0.0, pd[1] + pd[2] * X);
   if (m == 1) return std::max(0.0, pd[1]) + std::max(0.0, pd[2]) * PS::dens(X, pd[4]) * std::exp(-pd[3] * X);
   double thr = 0.5 * (1.0 + std::erf((X - SNUC) / (M_SQRT2 * pd[5])));
   return std::max(0.0, pd[1]) + std::max(0.0, pd[2]) * thr * std::exp(-pd[3] * X);
}
//  0-4 peak | 5-7 bkg abc | 8-10 A | 11-13 B | 14 mode | 15 Kb | 16 sigma_res
static double dTOT(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   double pd[6] = {p[14], p[5], p[6], p[7], p[15], p[16]};
   return bkgX(X, pd)
        + p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
        + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2))
        + p[8] * std::exp(-0.5 * std::pow((X - p[9]) / p[10], 2))
        + p[11] * std::exp(-0.5 * std::pow((X - p[12]) / p[13], 2));
}
static double dBKG(double* x, double* p) { double pd[6] = {p[0], p[1], p[2], p[3], p[4], p[5]}; return bkgX(x[0], pd); }
static double dGS(double* x, double* p)
{ double X = x[0], s = (X < p[1]) ? p[2] : p[3]; return p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2)); }
static double d52(double* x, double* p)
{ return p[4] * p[0] * std::exp(-0.5 * std::pow((x[0] - p[1] - E52) / p[3], 2)); }
static double dG3(double* x, double* p) { return p[0] * std::exp(-0.5 * std::pow((x[0] - p[1]) / p[2], 2)); }

struct Kit { TF1 *tot, *bk, *gs, *e5, *A, *B; Out o; TH1D* h; };

Kit build(const char* nm, double kb, double fw, double lo, int mode)
{
   Kit k;
   // the background systematic on the 5/2+ ratio: half the spread over the three background forms
   double rmin = 1e9, rmax = -1e9;
   for (int m = 0; m < 3; m++) { Out t = one(nm, kb, fw, lo, m, false);
      rmin = std::min(rmin, t.r52); rmax = std::max(rmax, t.r52); }
   k.o = one(nm, kb, fw, lo, mode, false);
   k.o.r52sys = 0.5 * (rmax - rmin);
   double sres = fw / 2.355;
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   k.h = (TH1D*)f.Get(nm); k.h->SetDirectory(nullptr);

   k.tot = new TF1(Form("T_%s", nm), dTOT, lo, 20.0, 17);
   for (int i = 0; i < 14; i++) k.tot->SetParameter(i, k.o.par[i]);
   k.tot->SetParameter(14, mode); k.tot->SetParameter(15, kb); k.tot->SetParameter(16, sres);
   k.tot->SetLineColor(kRed); k.tot->SetLineWidth(2); k.tot->SetNpx(2000);

   k.bk = new TF1(Form("K_%s", nm), dBKG, lo, 20.0, 6);
   k.bk->SetParameters(mode, k.o.par[5], k.o.par[6], k.o.par[7], kb, sres);
   k.bk->SetLineColor(kGray + 2); k.bk->SetLineStyle(2); k.bk->SetLineWidth(2); k.bk->SetNpx(1500);

   k.gs = new TF1(Form("G_%s", nm), dGS, lo, 3.2, 5);
   k.e5 = new TF1(Form("E_%s", nm), d52, lo, 3.6, 5);
   for (int i = 0; i < 5; i++) { k.gs->SetParameter(i, k.o.par[i]); k.e5->SetParameter(i, k.o.par[i]); }
   k.gs->SetLineColor(kAzure + 2);  k.gs->SetLineStyle(7); k.gs->SetLineWidth(2); k.gs->SetNpx(1200);
   k.e5->SetLineColor(kViolet + 1); k.e5->SetLineStyle(1); k.e5->SetLineWidth(3); k.e5->SetNpx(1200);

   k.A = new TF1(Form("A_%s", nm), dG3, 2.0, 20.0, 3);
   k.A->SetParameters(k.o.par[8], k.o.par[9], k.o.par[10]);
   k.A->SetLineColor(kGreen + 2); k.A->SetLineStyle(7); k.A->SetLineWidth(2); k.A->SetNpx(800);
   k.B = new TF1(Form("B_%s", nm), dG3, 2.0, 20.0, 3);
   k.B->SetParameters(k.o.par[11], k.o.par[12], k.o.par[13]);
   k.B->SetLineColor(kOrange + 7); k.B->SetLineStyle(7); k.B->SetLineWidth(2); k.B->SetNpx(800);
   return k;
}

void wide(Kit& k, const char* nm, double lo, int mode)
{
   TH1D* h = k.h;
   h->SetTitle(Form("%s  #font[12]{E}_{x}(^{15}C) from ^{15}C(d,d');E_{x} (MeV);counts / 100 keV", nm));
   h->SetMarkerStyle(20); h->SetMarkerSize(0.45); h->SetLineColor(kBlack);
   h->GetXaxis()->SetRangeUser(lo, 20.0);
   h->SetMinimum(8); h->SetMaximum(h->GetMaximum() * 4);
   h->Draw("E");
   k.bk->Draw("same"); k.gs->Draw("same"); k.e5->Draw("same");
   k.A->Draw("same"); k.B->Draw("same"); k.tot->Draw("same"); h->Draw("E same");

   TLine* ln = new TLine(SNUC, 8, SNUC, h->GetMaximum() * 0.25); ln->SetLineStyle(3); ln->Draw();
   TLatex l; l.SetTextSize(0.030); l.DrawLatex(SNUC + 0.15, h->GetMaximum() * 0.28, "S_{n}");
   Out& o = k.o;
   TLatex t; t.SetNDC(true); t.SetTextSize(0.032);
   t.DrawLatex(0.52, 0.885, Form("#chi^{2}/ndf = %.2f", o.chi2ndf));
   t.SetTextColor(kAzure + 2);  t.DrawLatex(0.52, 0.848, Form("g.s. 1/2^{+}  #mu = %+.3f  FWHM %.2f/%.2f", o.mu, o.fwL, o.fwR));
   t.SetTextColor(kViolet + 1); t.DrawLatex(0.52, 0.811, Form("5/2^{+} (0.740)  amp/g.s. = %.3f", o.r52));
   t.SetTextColor(kGreen + 2);  t.DrawLatex(0.52, 0.774, Form("A  %.3f #pm %.3f  FWHM %.2f", o.Amu, o.dAmu, o.Afw));
   t.SetTextColor(kOrange + 7); t.DrawLatex(0.52, 0.737, Form("B  %.3f #pm %.3f  FWHM %.2f", o.Bmu, o.dBmu, o.Bfw));
   t.SetTextColor(kGray + 2);   t.DrawLatex(0.52, 0.700, Form("bkg  %s", mode == 0 ? "c_{0} + c_{1}E_{x}  (empirical)"
                                     : (mode == 1 ? "floor + phase space #times free atten." : "floor + threshold")));
}

// The 5/2+ at 0.740 MeV is a SHOULDER on the g.s., not a separate peak -- on a linear axis next to
// a 9000-count g.s. it is invisible. So subtract the fitted g.s. and the background and plot the
// residual: the first excited state then appears as the bump it is, with the fitted 5/2+ Gaussian
// (which was NOT fitted to this residual -- it comes from the full fit) drawn over it.
void zoom(Kit& k, const char* nm)
{
   Out& o = k.o;
   double lo = o.mu - 1.4, hi = o.mu + 2.6;
   TH1D* hr = (TH1D*)k.h->Clone(Form("hr_%s", nm));
   hr->Reset();
   hr->SetTitle(Form("%s  --  data #minus (g.s. + background):  the 0.740 MeV 5/2^{+};"
                     "E_{x} (MeV);counts / 100 keV  (g.s. and bkg subtracted)", nm));
   double top = 0;
   for (int i = 1; i <= k.h->GetNbinsX(); i++) {
      double x = k.h->GetBinCenter(i);
      if (x < lo || x > hi) continue;
      double y = k.h->GetBinContent(i) - k.gs->Eval(x) - k.bk->Eval(x);
      hr->SetBinContent(i, y); hr->SetBinError(i, k.h->GetBinError(i));
      top = std::max(top, y + k.h->GetBinError(i));
   }
   hr->SetMarkerStyle(20); hr->SetMarkerSize(0.7); hr->SetLineColor(kBlack);
   hr->GetXaxis()->SetRangeUser(lo, hi);
   hr->SetMinimum(-0.50 * top); hr->SetMaximum(1.45 * top);
   hr->Draw("E");

   TF1* e5f = (TF1*)k.e5->Clone(Form("e5f_%s", nm));
   e5f->SetRange(lo, hi);
   e5f->SetFillColorAlpha(kViolet + 1, 0.35); e5f->SetFillStyle(1001);
   e5f->Draw("same FC");
   TLine* z = new TLine(lo, 0, hi, 0); z->SetLineStyle(2); z->SetLineColor(kGray + 2); z->Draw();
   hr->Draw("E same");

   TArrow* ar = new TArrow(o.mu + E52, top * 1.12, o.mu + E52, top * 1.02, 0.022, "|>");
   ar->SetLineColor(kViolet + 1); ar->SetFillColor(kViolet + 1); ar->SetLineWidth(3); ar->Draw();
   TLatex li; li.SetTextSize(0.040); li.SetTextColor(kViolet + 1);
   li.DrawLatex(o.mu + E52 - 0.55, top * 1.18, "5/2^{+},  0.740 MeV above the g.s.");
   li.SetTextColor(kBlack); li.SetTextSize(0.034);
   li.DrawLatex(lo + 0.08, top * 1.30, Form("g.s. #mu = %+.3f MeV   amp ratio %.4f #pm %.4f (stat) #pm %.4f (bkg)",
                                            o.mu, o.r52, o.dr52, o.r52sys));
   li.SetTextColor(kGray + 2); li.SetTextSize(0.030);
   li.DrawLatex(lo + 0.08, -0.40 * top, "the 5/2^{+} curve is from the FULL fit, not refitted to this residual");
}

void plotglobal()
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.045, "t");
   // LIN adopted for BOTH. PSATT fits better but its advantage comes from a free exp(-k*Ex)
   // attenuation that the 2026-11-09 study showed cannot be efficiency -- it is a shape
   // systematic, not the adopted model. See LEVELS_15C.md.
   Kit ka = build("ATTPCROOT", 191.25, 1.012, -1.6, 0);
   Kit ks = build("Spyral",    195.00, 0.548, -1.2, 0);

   TCanvas* c = new TCanvas("c", "", 1900, 1250);
   c->Divide(2, 2);
   c->cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03); wide(ka, "ATTPCROOT", -1.6, 0);
   c->cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.105); gPad->SetRightMargin(0.03); wide(ks, "Spyral", -1.2, 0);
   c->cd(3); gPad->SetLeftMargin(0.115); gPad->SetRightMargin(0.03); zoom(ka, "ATTPCROOT");
   c->cd(4); gPad->SetLeftMargin(0.115); gPad->SetRightMargin(0.03); zoom(ks, "Spyral");
   c->SaveAs("/home/yassid/Data/C15d_statefit/global_ex_fit.png");
   printf("wrote global_ex_fit.png\n");
}
