#include "binfit.C"
// Is the strength at Ex ~ 3.4 a missing LEVEL or missing CONTINUUM?
// Compare, in the same bins and window, adding:
//   (A) one barrier BW at 3.4, Gamma free       -> a level
//   (B) one barrier BW at 3.4, Gamma FIXED 0.2  -> a narrow level
//   (C) a smooth extra threshold term  erf((Ex-Sn)/w) * exp(-k Ex)  -> more continuum
// A level gives a narrow Gamma and beats the smooth term; continuum does the opposite.
static T2 TX2; static double SMW = 0.5, SMK = 0.3;
double fA(double* x, double* p)   // BW with free Gamma (rebuilt outside), amplitude p[9]
{ return fBin(x, p) + std::max(0.0, p[9]) * TX2.Eval(x[0]); }
double fC(double* x, double* p)   // smooth threshold term: amp p[9], width p[10], slope p[11]
{ double X = x[0];
  double thr = 0.5 * (1 + std::erf((X - B2::SN) / (M_SQRT2 * std::max(0.05, p[10]))));
  return fBin(x, p) + std::max(0.0, p[9]) * thr * std::exp(-std::max(0.0, p[11]) * X); }
void missing()
{
   for (auto& k : CFG) {
      REFTOT = 0; for (int i = 0; i < NB; i++) REF[i] = 0;
      { double sa[NB] = {0}, tot = 0;
        for (int ib = 0; ib < 6; ib++) { BR o = fitBin(k, ib, -2.0, 12.0, false, false);
          if (o.N < 400) continue; for (int i = 0; i < NB; i++) sa[i] += o.A[i]; tot += o.N; }
        for (int i = 0; i < NB; i++) REF[i] = sa[i]; REFTOT = tot; }
      printf("\n================ %s ================\n", k.nm);
      printf("   bin   base      (A) BW free Gamma        (B) BW Gamma=0.2     (C) smooth continuum\n");
      for (int ib = 0; ib < 2; ib++) {
         BR o = fitBin(k, ib, -2.0, 12.0, false, true);
         double base = o.nll;
         TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
         TH1D* h0 = (TH1D*)f.Get(Form("%s_cm%d", k.nm, ib));
         auto seed = [&](TF1& g) {
            g.SetParameter(0, o.gsA); g.SetParLimits(0, 0, 1e9);
            g.SetParameter(1, o.e5A); g.SetParLimits(1, 0, 1e9);
            g.SetParameter(2, o.cs);  g.SetParLimits(2, 0, 1e7);
            for (int i = 0; i < NB; i++) { double s2 = REF[i] * o.N / REFTOT;
               g.SetParameter(3 + i, s2); g.SetParLimits(3 + i, 0, std::max(3.0 * s2, 0.4 * REF[i])); } };
         // (A) scan Gamma, keep the best
         double bestG = 0, bestA = 1e30;
         for (double G : {0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 2.5, 4.0}) {
            TX2.build(3.40 + GSMU, G, 1, SGR, -3.0, 23.0);
            TH1D* h = (TH1D*)h0->Clone(Form("a%s%d%d", k.nm, ib, (int)(G * 100))); h->SetDirectory(nullptr);
            TF1 g(Form("ga%s%d%d", k.nm, ib, (int)(G * 100)), fA, -2.0, 12.0, 10);
            seed(g); double t = TX2.Eval(3.40 + GSMU);
            g.SetParameter(9, t > 0 ? 10.0 / t : 1); g.SetParLimits(9, 0, 1e9);
            TFitResultPtr r = h->Fit(&g, "RQ0NLS");
            double v = 2 * r->MinFcnValue(); if (v < bestA) { bestA = v; bestG = G; } }
         // (B) Gamma fixed 0.2
         TX2.build(3.40 + GSMU, 0.20, 1, SGR, -3.0, 23.0);
         double vB; { TH1D* h = (TH1D*)h0->Clone(Form("b%s%d", k.nm, ib)); h->SetDirectory(nullptr);
            TF1 g(Form("gb%s%d", k.nm, ib), fA, -2.0, 12.0, 10);
            seed(g); double t = TX2.Eval(3.40 + GSMU);
            g.SetParameter(9, t > 0 ? 10.0 / t : 1); g.SetParLimits(9, 0, 1e9);
            TFitResultPtr r = h->Fit(&g, "RQ0NLS"); vB = 2 * r->MinFcnValue(); }
         // (C) smooth extra continuum
         double vC; { TH1D* h = (TH1D*)h0->Clone(Form("c%s%d", k.nm, ib)); h->SetDirectory(nullptr);
            TF1 g(Form("gc%s%d", k.nm, ib), fC, -2.0, 12.0, 12);
            seed(g);
            g.SetParameter(9, 20); g.SetParLimits(9, 0, 1e6);
            g.SetParameter(10, 0.5); g.SetParLimits(10, 0.05, 3.0);
            g.SetParameter(11, 0.2); g.SetParLimits(11, 0.0, 2.0);
            TFitResultPtr r = h->Fit(&g, "RQ0NLS"); vC = 2 * r->MinFcnValue(); }
         printf("   %2d  %8.1f   d %7.1f (Gamma %.2f)   d %7.1f        d %7.1f\n",
                ib, base, base - bestA, bestG, base - vB, base - vC);
      }
   }
}
