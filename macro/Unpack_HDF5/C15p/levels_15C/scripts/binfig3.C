#include "binfit.C"
// Per-theta_CM-bin fits: tailed response + six (p,p') barrier BWs + sequential phase-space
// continuum.  TGraphs (the templates are globals rebuilt per fit; TF1::Draw is lazy).
static TGraph* gcv(std::function<double(double)> f, double lo, double hi, int N = 1200)
{ auto* g = new TGraph(N);
  for (int i = 0; i < N; i++) { double x = lo + (hi - lo) * i / (N - 1.0); g->SetPoint(i, x, f(x)); }
  return g; }
void binfig3()
{
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1); gStyle->SetTitleSize(0.075,"t");
   int col[NB] = {kGreen + 2, kSpring - 6, kCyan + 2, kOrange + 7, kRed - 4, kMagenta + 2};
   TCanvas* c = new TCanvas("cb", "", 1800, 900); c->Divide(5, 2, 0.001, 0.001);
   for (int r = 0; r < 2; r++) {
      Cfg& k = CFG[r];
      REFTOT = 0; for (int i = 0; i < NB; i++) REF[i] = 0;
      { double sumA[NB] = {0}, tot = 0;
        for (int ib = 0; ib < 6; ib++) { BR o = fitBin(k, ib, -2.0, 9.0, false, false);
          if (o.N < 400) continue; for (int i = 0; i < NB; i++) sumA[i] += o.A[i]; tot += o.N; }
        for (int i = 0; i < NB; i++) REF[i] = sumA[i]; REFTOT = tot; }
      for (int ib = 0; ib < 5; ib++) {
         c->cd(r * 5 + ib + 1); gPad->SetLogy();
         gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.02);
         gPad->SetTopMargin(0.10); gPad->SetBottomMargin(0.14);
         USE_LEBKG = true; FREE_TURNON = true; TURNON = 2.5;
         BR o = fitBin(k, ib, -2.0, 9.0, false, true);
         if (o.N < 400 || !o.h) continue;
         // capture the fitted state
         double fl = o.floor, cs = o.cs, gsA = o.gsA, e5A = o.e5A;
         double gmu = GSMU, sl = SGL, sr = SGR, al = ALF_L, ar = ALF_R;
         std::vector<T2> T(NB); for (int i = 0; i < NB; i++) T[i] = TT[i];
         TGraph* GCl = GC;
         double bB = o.bgB, bW = o.bgW, bK = o.bgK;
         auto lebkg = [&](double X) { if (!USE_LEBKG) return 0.0;
            double c0 = FREE_TURNON ? TURNON : B2::SN;
            double thr = 0.5 * (1 + std::erf((X - c0) / (M_SQRT2 * std::max(0.05, bW))));
            return std::max(0.0, bB) * thr * std::exp(-std::max(0.0, bK) * X); };
         auto cont = [&](double X) { return fl + cs * std::max(0.0, GCl->Eval(X)) + lebkg(X); };
         auto ps   = [&](double X) { return fl + cs * std::max(0.0, GCl->Eval(X)); };
         auto tot  = [&](double X) { double v = cont(X)
                       + gsA * shapeGT(X, gmu, sl, sr, al, ar)
                       + e5A * shapeGT(X, gmu + E52d, sl, sr, al, ar);
                       for (int i = 0; i < NB; i++) v += o.A[i] * T[i].Eval(X); return v; };
         TH1D* h = o.h;
         h->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ%s", k.nm, o.cmLo, o.cmHi, ib==0?"  (degenerate)":""));
         h->SetMarkerStyle(20); h->SetMarkerSize(0.4); h->SetLineColor(kBlack);
         h->GetXaxis()->SetRangeUser(-2, 9);
         h->SetMinimum(0.5); h->SetMaximum(h->GetBinContent(h->GetMaximumBin()) * 4);
         h->GetXaxis()->SetTitle("E_{x} (MeV)"); h->GetYaxis()->SetTitle("counts / 100 keV");
         h->GetXaxis()->SetTitleSize(0.06); h->GetYaxis()->SetTitleSize(0.06);
         h->GetXaxis()->SetLabelSize(0.055); h->GetYaxis()->SetLabelSize(0.055);
         h->Draw("E");
         TGraph* gp = gcv(ps, -2, 9); gp->SetLineColor(kBlue + 1); gp->SetLineStyle(3);
         gp->SetLineWidth(2); gp->Draw("L");                       // phase space + floor alone
         TGraph* gk = gcv(cont, -2, 9); gk->SetLineColor(kGray + 2); gk->SetLineStyle(2);
         gk->SetLineWidth(3); gk->Draw("L");                       // + low-energy background
         for (int i = 0; i < NB; i++) { if (o.A[i] <= 0) continue;
            TGraph* gb = gcv([&, i](double X) { return cont(X) + o.A[i] * T[i].Eval(X); }, 2, 9);
            gb->SetLineColor(col[i]); gb->SetLineStyle(7); gb->SetLineWidth(2); gb->Draw("L"); }
         TGraph* gt = gcv(tot, -2, 9); gt->SetLineColor(kRed); gt->SetLineWidth(2); gt->Draw("L");
         h->Draw("E same");
         TLatex t; t.SetNDC(true); t.SetTextSize(0.055);
         t.DrawLatex(0.50, 0.87, Form("N %.0f", o.N));
         t.DrawLatex(0.50, 0.80, Form("#chi^{2}/ndf %.2f", o.chi2ndf));
         if (o.bgB > 1) { t.SetTextColor(kGray + 2);
            t.DrawLatex(0.50, 0.73, Form("bkg B %.0f", o.bgB)); t.SetTextColor(kBlack); }
      }
   }
   c->SaveAs("/home/yassid/Data/C15d_statefit/dd_cmbin_fits_ex9.png");
   printf("wrote dd_cmbin_fits_ex9.png\n");
}
