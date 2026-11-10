#include "exlib.h"
// Is the non-Gaussian g.s. lineshape intrinsic, or the SUM over theta slices of different
// resolution?  Fit each slice with the SAME model used in the global fit -- bifurcated Gaussian
// + the 5/2+ held at mu+0.740 + a linear background -- so the comparison is like for like.
// (A first attempt without the 5/2+ term gave chi2/bin 10-49 and was meaningless: the 5/2+ sits
//  inside the fit window, so omitting it inflates chi2 precisely where the statistics are.)
static const double E52 = 0.740;
double fpk(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   return p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
        + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2))
        + std::max(0.0, p[5] + p[6] * X);
}
void thSlice(const char* file, const Pars& p, const char* nm, double lo, double hi)
{
   TFile* f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("%s MISSING\n", file); return; }
   TTree* t = (TTree*)f->Get("dd");
   double ke, th, vz, c2, icd; float icf; int np;
   t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vz", &vz); t->SetBranchAddress("chi2ndf", &c2);
   bool icF = TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName()) == "Float_t";
   if (icF) t->SetBranchAddress("ic", &icf); else t->SetBranchAddress("ic", &icd);
   t->SetBranchAddress("npulse", &np);

   const int NS = 7; const double T0 = 64, DT = 4;          // 64-92 deg, where the g.s. lives
   TH1D* hs[NS]; TH1D* hall = new TH1D(Form("%s_all", nm), "", 280, -5, 23);
   hall->SetDirectory(nullptr);
   for (int i = 0; i < NS; i++) { hs[i] = new TH1D(Form("%s_s%d", nm, i), "", 280, -5, 23);
                                  hs[i]->SetDirectory(nullptr); }
   for (Long64_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      double icv = icF ? icf : icd, keu = ke + p.keOff;
      if (c2 > p.chi2 || keu < p.keLo || keu > p.keHi) continue;
      if (vz < p.vzLo || vz > p.vzHi || icv < p.icLo || icv > p.icHi) continue;
      if (np < p.npLo || np > p.npHi) continue;
      double thf = std::fabs(th); if (thf > 90) thf = 180 - thf;
      double thc = thCorr(p, thf, keu);
      if (thc < p.thLo || thc > p.thHi) continue;
      double e = exOf(p, ebAt(p, vz, true), thc, keu);
      if (!std::isfinite(e)) continue;
      hall->Fill(e);
      int k = (int)((thc - T0) / DT);
      if (k >= 0 && k < NS) hs[k]->Fill(e);
   }
   f->Close();

   printf("\n== %s   g.s. + 5/2+ fitted separately per theta slice, window [%.1f, %.1f]\n", nm, lo, hi);
   printf("   theta       N      mu       FWHM L/R      5/2+ ratio   chi2/bin\n");
   auto fit1 = [&](TH1D* h, const char* lab) {
      if (h->GetEntries() < 800) { printf("   %-9s %7.0f   (too few)\n", lab, h->GetEntries()); return std::make_pair(0.0, 0.0); }
      TF1 g(Form("g_%s_%s", nm, lab), fpk, lo, hi, 7);
      g.SetParameter(0, h->GetMaximum());  g.SetParLimits(0, 0, 1e9);
      g.SetParameter(1, 0);                g.SetParLimits(1, -1.0, 1.0);
      g.SetParameter(2, 0.3);              g.SetParLimits(2, 0.05, 1.5);
      g.SetParameter(3, 0.3);              g.SetParLimits(3, 0.05, 1.5);
      g.SetParameter(4, 0.08);             g.SetParLimits(4, 0, 0.8);
      g.SetParameter(5, 10); g.SetParameter(6, 0);
      h->Fit(&g, "RQ0NL");
      double c = 0; int n = 0;
      for (int i = 1; i <= h->GetNbinsX(); i++) {
         double x = h->GetBinCenter(i), y = h->GetBinContent(i);
         if (x < lo || x > hi || y < 10) continue;
         double m = g.Eval(x); c += (y - m) * (y - m) / y; n++; }
      printf("   %-9s %7.0f  %+.3f   %.3f/%.3f   %.4f       %6.2f\n", lab, h->GetEntries(),
             g.GetParameter(1), 2.355 * g.GetParameter(2), 2.355 * g.GetParameter(3),
             g.GetParameter(4), c / std::max(1, n));
      return std::make_pair(2.355 * g.GetParameter(3), c / std::max(1, n));
   };
   for (int i = 0; i < NS; i++) fit1(hs[i], Form("%.0f-%.0f", T0 + i * DT, T0 + (i + 1) * DT));
   printf("   ---------------------------------------------------------------\n");
   fit1(hall, "ALL");
}
void thetaslice2()
{
   Pars A; A.ebeam = 191.25; A.chi2 = 0.2; A.thLo = 0; A.thHi = 180; A.keLo = 0.5; A.keHi = 18;
   A.vzLo = 0; A.vzHi = 1100; A.icLo = 1074; A.icHi = 1201;
   A.kcDenom = 3913.0; A.kcPivot = -1.5; A.keOff = 0; A.kcOn = true;
   A.ebzSlope = 11.0; A.ebzPivot = 500; A.ebzOn = true;
   Pars S = A; S.ebeam = 195; S.chi2 = 5; S.keLo = 0; S.keHi = 20; S.vzLo = 0; S.vzHi = 400;
   S.icLo = -1; S.icHi = 99999; S.kcDenom = 2535.2; S.kcPivot = 14.5;
   thSlice("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root", A, "ATTPCROOT", -1.6, 2.2);
   thSlice("/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root", S, "Spyral", -1.2, 2.2);
}
