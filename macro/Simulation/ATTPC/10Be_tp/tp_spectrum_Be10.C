/// @file tp_spectrum_Be10.C
/// @brief The 12Be excitation-energy spectrum the AT-TPC would actually measure: the four bound
///        levels SUMMED with their relative populations, and the test of whether the 0+_2 at
///        2.251 MeV is recoverable from it.
///
///   root -b -q 'tp_spectrum_Be10.C("/mnt/f/Be10_tp","b285_attpc")'
///   root -b -q 'tp_spectrum_Be10.C("/mnt/f/Be10_tp","b700_2mm",0.2,2.0,45.0)'
///
/// POPULATIONS. 0+ g.s., 2+ 2.109 and 1- 2.715 enter with equal weight; the intruder 0+_2 at
/// 2.251 enters 5x weaker (w2251 = 0.2). The suppression is applied HERE and not in the
/// generation: every level was simulated at full statistics so its acceptance and its resolution
/// are each measured to the same precision, and a relative population is a multiplicative weight
/// on a spectrum, not a property of the detector. The weight actually applied is pop_i / nGen_i,
/// with nGen_i READ from the acceptance file rather than assumed equal, so each level's own
/// acceptance survives into the sum and the histogram is proportional to what a detector counts.
///
/// TWO THINGS THIS MACRO GOT WRONG THE FIRST TIME, both now fixed, because both are traps that
/// would silently produce a believable-looking wrong answer:
///
/// 1. FITTING ALL ANGLES AT ONCE. sigma(Ex) runs from 0.19 MeV at theta_cm 2-20 to 1.49 MeV at
///    theta_cm 90-180 (summary_Be10tp.C). Summed over all angles the lineshape is a superposition
///    of widths, not a gaussian, and a four-gaussian fit to it does not find four levels -- it
///    uses the spare gaussians as background. It did exactly that: the 0+_2 component came back
///    with sigma 1.77 MeV and a RAILED centroid at 1.25, the 1- amplitude went to zero, and the
///    likelihood "preferred" the four-peak model by 161 units for entirely spurious reasons.
///    FIT IN theta_cm SLICES.
///
/// 2. LETTING EVERY CENTROID AND WIDTH FLOAT. An experiment knows the 12Be level scheme; what it
///    does not know is how much of each level it produced. So the model here is the physical one:
///    the four level energies are FIXED, one common energy-scale shift and one common resolution
///    width are fitted, and the four AREAS are free and non-negative. Six parameters, not twelve.
///
/// THE TEST. Two fits, differing by one parameter:
///     FIT A  all four levels, the 0+_2 area free
///     FIT B  the 0+_2 removed entirely -- the null hypothesis
/// compared on the LIKELIHOOD (option "L" minimises -2lnL; chi2 is not the objective of an "L"
/// fit). They are run on the ASIMOV dataset -- the expected spectrum with no Poisson fluctuation
/// -- so what comes out is the EXPECTED significance of the 0+_2 for a run of nTot counts, not
/// the significance of one lucky realisation. sqrt(Delta(-2lnL)) is that significance in sigma,
/// and it scales as sqrt(nTot), so a result can be read off for any beam time.
///
/// Ex is computed BOTH ways: with one constant beam energy for every vertex, and with the beam
/// energy at the RECONSTRUCTED vertex, whose z profile is solved out of MC truth here. The second
/// costs nothing but software and it is worth a factor 2.6 in this channel (dEx/dE_beam at the
/// transfer peak is 0.12-0.14, against 0.043-0.053 for 14C(d,p)).

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

static const int NL = 4;
static const char *LVL[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
static const double LEX[NL] = {0.0, 2.109, 2.251, 2.715};
static const char *LJP[NL] = {"0+", "2+", "0+_2", "1-"};

// 10Be(t,p)12Be
static const double SU = 931.49401;
static const double SM1 = 10.0135341 * SU, SM2 = 3.0160493 * SU, SM3 = 1.007825 * SU, SM4 = 12.0269221 * SU;

static double sp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double sp_ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + SM1, Et3 = Ke + SM3;
   double s = SM1 * SM1 + SM2 * SM2 + 2 * SM2 * Et1;
   double uu = SM2 * SM2 + SM3 * SM3 - 2 * SM2 * Et3;
   double a = (std::cos(th) * sp_om2(s, SM1 * SM1, SM2 * SM2) * sp_om2(uu, SM2 * SM2, SM3 * SM3) -
               (s - SM1 * SM1 - SM2 * SM2) * (SM2 * SM2 + SM3 * SM3 - uu)) / (2 * SM2 * SM2) + s + uu - SM2 * SM2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
static TString sp_find(const TString &dir, const TString &cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg.Data(), pre, lvl, cfg.Data()));
   return f.Strip(TString::kBoth);
}
static double sp_ngen(const TString &dir, const TString &cfg, const char *lvl)
{
   TString fa = sp_find(dir, cfg, lvl, "acceptance");
   if (fa.IsNull()) return -1;
   TFile *f = TFile::Open(fa);
   if (!f || f->IsZombie()) return -1;
   double n = -1;
   TIter nx(f->GetListOfKeys());
   while (auto *k = (TKey *)nx()) {
      TString nm = k->GetName();
      if (nm.BeginsWith("hGen_")) { TH1D *h = (TH1D *)f->Get(nm); if (h) n = h->Integral(); }
   }
   f->Close();
   return n;
}

/// THE PEAK SHAPE, and it is not a gaussian.
///
/// Measured on a single level with the vertex beam energy (theta_cm 2-45, 2.85 T AT-TPC): 65 % of
/// the events land within +-0.10 MeV of the centroid and 81 % within +-0.20, but ~19 % sit in a
/// tail reaching past +-0.5. The robust core width is 0.093 MeV while a single gaussian fitted to
/// the whole thing returns 0.21 -- it splits the difference between core and tail and describes
/// neither. Using that single gaussian is what made the first version of this macro report a
/// resolution twice too poor and then, correctly for the wrong reason, find nothing at 2.251.
///
/// So the shape is a CORE plus a TAIL gaussian, both centred on the same (fixed) level energy,
/// with the shift, the two widths and the tail fraction COMMON to all four levels -- because they
/// are properties of the detector, not of the level -- and only the four areas free. That is also
/// what an experiment can do: it knows the 12Be level scheme and it can measure its own lineshape
/// on the strong, isolated ground state.
///
/// p = {shift, sigmaCore, sigmaTail, fTail, A0 .. A(n-1)}
static int gNPk = 4;
static double gPkEx[NL];
static const int NSHAPE = 4;
/// index of the peak that is allowed to move and broaden (fit C's harder null), or -1 for none.
/// When set, two extra parameters follow the amplitudes: a centroid offset and an extra width
/// added in quadrature, applied to that peak alone.
static int gFreePk = -1;
static double sp_model(double *x, double *p)
{
   double sc = p[1], st = p[2], ft = p[3], v = 0;
   if (sc <= 0 || st <= 0) return 0;
   for (int i = 0; i < gNPk; ++i) {
      double mu = gPkEx[i] + p[0], sci = sc, sti = st;
      if (i == gFreePk) {
         double dmu = p[NSHAPE + gNPk], extra = p[NSHAPE + gNPk + 1];
         mu += dmu;
         sci = std::sqrt(sc * sc + extra * extra);
         sti = std::sqrt(st * st + extra * extra);
      }
      double d0 = (x[0] - mu) / sci, d1 = (x[0] - mu) / sti;
      v += p[NSHAPE + i] * ((1 - ft) * std::exp(-0.5 * d0 * d0) + ft * (sci / sti) * std::exp(-0.5 * d1 * d1));
   }
   return v;
}

void tp_spectrum_Be10(TString root = "/mnt/f/Be10_tp", TString cfg = "b285_attpc", Double_t w2251 = 0.2,
                      Double_t cmLo = 2.0, Double_t cmHi = 45.0, Double_t nTot = 20000, Bool_t vertexEbeam = kTRUE,
                      Double_t Ebeam = 112.20, TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull()) outDir = gSystem->pwd();
   const double pop[NL] = {1.0, 1.0, w2251, 1.0};
   // the fit window: wide enough to hold all four levels and their tails, narrow enough that the
   // fit is not describing empty axis
   const double exLo = -1.5, exHi = 4.5;
   const int nBins = 120;

   printf("\n=========== 10Be(t,p)12Be spectrum : %s / %s ===========\n", root.Data(), cfg.Data());
   printf("theta_cm %.0f-%.0f deg | populations 1 : 1 : %.2f : 1 | Ex from %s beam energy | %.0f counts\n", cmLo, cmHi,
          w2251, vertexEbeam ? "the RECONSTRUCTED-VERTEX" : "a CONSTANT", nTot);

   TH1D *hLev[NL], *hSum = new TH1D("hSum", ";E_{x}(^{12}Be) [MeV];counts / bin", nBins, exLo, exHi);
   hSum->SetDirectory(nullptr);
   double nGen[NL], nUse[NL];
   for (int l = 0; l < NL; ++l) {
      hLev[l] = new TH1D(Form("hL%d", l), "", nBins, exLo, exHi);
      hLev[l]->SetDirectory(nullptr);
      nGen[l] = sp_ngen(root, cfg, LVL[l]);
      nUse[l] = 0;
      TString f = sp_find(root, cfg, LVL[l], "exres");
      if (f.IsNull() || nGen[l] <= 0) {
         printf("\033[1;31mMISSING level %s in %s -- nothing to fit\033[0m\n", LVL[l], cfg.Data());
         return;
      }
      TFile *fr = TFile::Open(f);
      TTree *t = (TTree *)fr->Get("res");
      double exReco, cmTrue, thTrue, thReco, keTrue, keReco, zTrue, zReco;
      t->SetBranchAddress("exReco", &exReco);
      t->SetBranchAddress("cmTrue", &cmTrue);
      t->SetBranchAddress("thTrue", &thTrue);
      t->SetBranchAddress("thReco", &thReco);
      t->SetBranchAddress("keTrue", &keTrue);
      t->SetBranchAddress("keReco", &keReco);
      t->SetBranchAddress("zTrue", &zTrue);
      t->SetBranchAddress("zReco", &zReco);

      // E_beam(z), solved from TRUTH over all angles -- it is a property of the gas and the vertex
      // distribution, not of the angular slice, so fitting it per slice would only add noise.
      TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
      if (vertexEbeam) {
         std::vector<double> eb, zz;
         const double mres = SM4 + LEX[l];
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = 80., hi = 130.;
            double flo = sp_ex(mres, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = sp_ex(mres, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double m = 0.5 * (lo + hi), fm = sp_ex(mres, m, thTrue * TMath::DegToRad(), keTrue);
               if (std::isnan(fm)) break;
               if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
            }
            double e = 0.5 * (lo + hi);
            if (e > 85 && e < 125) { eb.push_back(e); zz.push_back(zTrue); }
         }
         if (eb.size() < 100) { printf("\033[1;31mtoo few events to fit E_beam(z) for %s\033[0m\n", LVL[l]); return; }
         TGraph g((int)eb.size(), zz.data(), eb.data());
         fEb.SetParameters(Ebeam, -0.006, 0.);
         g.Fit(&fEb, "QN");
      }

      const double w = pop[l] / nGen[l];
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (cmTrue < cmLo || cmTrue >= cmHi) continue;
         double ex = exReco;
         if (vertexEbeam) {
            ex = sp_ex(SM4, fEb.Eval(zReco), thReco * TMath::DegToRad(), keReco);
            if (std::isnan(ex)) continue;
         }
         hLev[l]->Fill(ex, w);
         nUse[l] += 1;
      }
      fr->Close();
      hSum->Add(hLev[l]);
   }

   printf("\n  %-8s %-6s %8s %8s %12s %14s\n", "level", "J^pi", "nGen", "nFit", "acceptance", "yield fraction");
   double tot = hSum->Integral();
   for (int l = 0; l < NL; ++l)
      printf("  %-8s %-6s %8.0f %8.0f %12.3f %14.4f\n", LVL[l], LJP[l], nGen[l], nUse[l], nUse[l] / nGen[l],
             tot > 0 ? hLev[l]->Integral() / tot : 0.0);
   if (!(tot > 0)) { printf("empty slice\n"); return; }
   // the TRUE injected ratio in this slice -- computed, not assumed, because acceptance is
   // level-dependent and the fit will be compared against this and not against w2251
   double rTrue = hLev[1]->Integral() > 0 ? hLev[2]->Integral() / hLev[1]->Integral() : NAN;

   // ASIMOV dataset: the expected spectrum, scaled to nTot counts, with NO Poisson fluctuation.
   // Fitting this returns the EXPECTED significance rather than that of one realisation.
   hSum->Scale(nTot / tot);
   for (int l = 0; l < NL; ++l) hLev[l]->Scale(nTot / tot);

   // THE GLOBALS ARE THE TRAP. sp_model reads gNPk / gPkEx / gFreePk, so a TF1 built for one fit
   // evaluates with WHATEVER configuration was set last. Fit C runs after fit A and leaves gNPk = 3
   // with a free peak, so drawing fA afterwards drew a curve that was not fit A at all -- it showed
   // the total falling BELOW one of its own components, which is what gave the bug away. The fitted
   // numbers were never affected (they are captured at fit time), only anything evaluated later.
   // setShape() restores a named configuration; call it before every evaluation of a stored TF1.
   auto setShape = [&](bool with2251, int freePk) {
      gNPk = 0;
      int idx2109 = -1;
      for (int l = 0; l < NL; ++l) {
         if (l == 2 && !with2251) continue;
         if (l == 1) idx2109 = gNPk;
         gPkEx[gNPk++] = LEX[l];
      }
      gFreePk = (freePk >= 0) ? idx2109 : -1;
   };
   auto runFit = [&](TH1D *h, bool with2251, const char *tag, int freePk = -1) {
      setShape(with2251, freePk);
      int npar = NSHAPE + gNPk + (gFreePk >= 0 ? 2 : 0);
      TF1 *f = new TF1(Form("f%s", tag), sp_model, exLo, exHi, npar);
      f->SetParName(0, "shift");
      f->SetParName(1, "sigmaCore");
      f->SetParName(2, "sigmaTail");
      f->SetParName(3, "fTail");
      f->SetParameter(0, 0.0);
      f->SetParLimits(0, -0.8, 0.8);
      f->SetParameter(1, 0.10);
      f->SetParLimits(1, 0.02, 0.60);
      f->SetParameter(2, 0.45);
      f->SetParLimits(2, 0.10, 3.0);
      f->SetParameter(3, 0.25);
      f->SetParLimits(3, 0.0, 0.90);
      for (int i = 0; i < gNPk; ++i) {
         f->SetParName(NSHAPE + i, Form("A_%.3f", gPkEx[i]));
         f->SetParameter(NSHAPE + i, h->GetMaximum() * 0.6);
         f->SetParLimits(NSHAPE + i, 0, 1e8);
      }
      if (gFreePk >= 0) {
         f->SetParName(NSHAPE + gNPk, "dmu_2109");
         f->SetParameter(NSHAPE + gNPk, 0.0);
         f->SetParLimits(NSHAPE + gNPk, -0.30, 0.30);
         f->SetParName(NSHAPE + gNPk + 1, "extraSig_2109");
         f->SetParameter(NSHAPE + gNPk + 1, 0.03);
         f->SetParLimits(NSHAPE + gNPk + 1, 0.0, 0.50);
      }
      f->SetNpx(1200);
      auto r = h->Fit(f, "LSQ0", "", exLo, exHi);
      return std::make_pair(f, r);
   };

   auto A = runFit(hSum, true, "A");
   auto B = runFit(hSum, false, "B");
   // FIT C : the HARDER null. Three levels, but the 2.109 is allowed to slide and to broaden --
   // exactly the freedom an unresolved doublet would exploit. If A still wins against C, the
   // spectrum is telling us there are two states there and not one badly-modelled one. Without
   // this test, "A beats B" only says the 2.109 peak is not a perfect single lineshape, which a
   // slightly wrong calibration would also produce.
   auto C = runFit(hSum, false, "C", 1);
   TF1 *fA = A.first;
   double lA = A.second.Get() ? A.second->MinFcnValue() : NAN;
   double lB = B.second.Get() ? B.second->MinFcnValue() : NAN;
   double lC = C.second.Get() ? C.second->MinFcnValue() : NAN;

   printf("\n---- FIT A : four levels at FIXED energies, common shift and width, free areas ----\n");
   printf("  %-12s %12s %12s %10s\n", "parameter", "value", "error", "");
   for (int p = 0; p < fA->GetNpar(); ++p) {
      double lo, hi;
      fA->GetParLimits(p, lo, hi);
      double v = fA->GetParameter(p);
      // Flag on the VALUE, not on the range. Scaling the tolerance by (hi-lo) made every
      // amplitude "railed" because their upper bound is 1e8: a perfectly healthy A = 611 sat
      // within 1e-4*1e8 = 10^4 of zero. Compare against the parameter's own error instead, which
      // is what "sitting on a bound" actually means.
      double tol = std::max(1e-9, 0.05 * std::fabs(fA->GetParError(p)));
      bool rail = (hi > lo) && (std::fabs(v - lo) < tol || std::fabs(v - hi) < tol);
      printf("  %-12s %12.5f %12.5f %10s\n", fA->GetParName(p), v, fA->GetParError(p), rail ? "<-- RAILED" : "");
   }
   // AREA. The shape integrates to sigmaCore*sqrt(2pi) per unit A whatever the tail fraction is,
   // because the tail term carries the factor (sc/st) precisely so that its integral is also
   // sc*sqrt(2pi): (1-f)*sc + f*(sc/st)*st = sc. So the area is A * sc * sqrt(2pi) / binwidth and
   // the four areas stay comparable to each other however the tail comes out.
   double bw = hSum->GetBinWidth(1), sg = fA->GetParameter(1);
   double area[NL];
   for (int l = 0; l < NL; ++l) area[l] = fA->GetParameter(NSHAPE + l) * sg * std::sqrt(2 * TMath::Pi()) / bw;
   // GUARD. A ratio is only meaningful while its denominator is a real peak. In a slice where the
   // fit collapses, A_2109 goes to ~0 and area/area(2+) prints something like 7e8 -- a number that
   // looks like a result and is not one. Require the 2+ to hold at least 1 % of the spectrum.
   const bool refOK = area[1] > 0.01 * nTot;
   printf("\n  %-8s %10s %12s %14s %14s\n", "level", "Ex [MeV]", "area", "area/area(2+)", "input ratio");
   for (int l = 0; l < NL; ++l) {
      if (refOK)
         printf("  %-8s %10.3f %12.1f %14.3f %14s\n", LVL[l], LEX[l], area[l], area[l] / area[1],
                l == 2 ? Form("%.3f", rTrue) : "-");
      else
         printf("  %-8s %10.3f %12.1f %14s %14s\n", LVL[l], LEX[l], area[l], "n/a", l == 2 ? Form("%.3f", rTrue) : "-");
   }
   if (!refOK)
      printf("  \033[1;31mTHE 2+ REFERENCE PEAK COLLAPSED (area %.1f of %.0f counts) -- every ratio below is\n"
             "  meaningless and this slice must not be quoted.\033[0m\n", area[1], nTot);
   printf("  fitted lineshape: core sigma %.3f MeV, tail sigma %.3f MeV with %.1f %% of the yield,\n"
          "                    energy-scale shift %+.3f MeV\n",
          sg, fA->GetParameter(2), 100 * fA->GetParameter(3), fA->GetParameter(0));

   double dl = lB - lA;
   double nsig = dl > 0 ? std::sqrt(dl) : 0.0;
   printf("\n---- IS THE 0+_2 THERE? FIT A vs FIT B (identical but for the 2.251 area) ----\n");
   printf("  A (with 0+_2, %d par) : -2lnL = %.3f\n", fA->GetNpar(), lA);
   printf("  B (no   0+_2, %d par) : -2lnL = %.3f\n", B.first->GetNpar(), lB);
   printf("  \033[1;3%dm Delta(-2lnL) = %.2f for ONE extra parameter  ->  expected significance %.1f sigma\033[0m\n",
          nsig >= 3 ? 2 : (nsig >= 2 ? 3 : 1), dl, nsig);
   printf("  at %.0f counts in this slice. It scales as sqrt(N): %.1f sigma at 10k, %.1f at 100k.\n", nTot,
          nsig * std::sqrt(10000. / nTot), nsig * std::sqrt(100000. / nTot));
   double dlC = lC - lA;
   double nsigC = dlC > 0 ? std::sqrt(dlC) : 0.0;
   printf("\n  the HARDER null, FIT C: no 0+_2, but the 2.109 free to slide (%+.3f MeV) and broaden\n",
          C.first->GetParameter(NSHAPE + C.first->GetNpar() - NSHAPE - 2));
   printf("  (extra sigma %.3f MeV) -- the freedom an unresolved doublet would exploit\n",
          C.first->GetParameter(C.first->GetNpar() - 1));
   printf("  C (%d par) : -2lnL = %.3f   ->  A beats it by %.2f  =  \033[1;3%dm%.1f sigma\033[0m\n",
          C.first->GetNpar(), lC, dlC, nsigC >= 3 ? 2 : (nsigC >= 2 ? 3 : 1), nsigC);
   printf("  THIS is the number to quote: A vs B only says the 2.109 is not a single clean peak,\n"
          "  which a slightly wrong calibration would also produce. A vs C says there are TWO states.\n");
   if (refOK)
      printf("  fitted 0+_2 / 2+ area ratio %.3f against the %.3f actually injected%s\n", area[2] / area[1], rTrue,
             (std::fabs(area[2] / area[1] - rTrue) < 0.35 * rTrue) ? "  (recovered)" : "  (NOT recovered)");
   else
      printf("  no usable area ratio: the 2+ reference peak collapsed in this slice\n");

   // ---- figure, on ONE Poisson realisation so it looks like data -----------------------------
   gRandom->SetSeed(20260901);
   TH1D *hObs = (TH1D *)hSum->Clone("hObs");
   hObs->SetDirectory(nullptr);
   for (int b = 1; b <= hObs->GetNbinsX(); ++b) hObs->SetBinContent(b, gRandom->Poisson(hSum->GetBinContent(b)));
   hObs->SetTitle(Form("^{10}Be(t,p)^{12}Be  %s  #theta_{cm} %.0f-%.0f#circ  (%s E_{beam})", cfg.Data(), cmLo, cmHi,
                       vertexEbeam ? "vertex" : "const"));
   TCanvas *c = new TCanvas("cSpec", "spectrum", 1100, 750);
   hObs->SetLineColor(kBlack);
   hObs->SetMarkerStyle(20);
   hObs->SetMarkerSize(0.6);
   double ymax = hObs->GetMaximum() * 1.35;
   hObs->GetYaxis()->SetRangeUser(0, ymax);
   hObs->Draw("E");
   // Sample each fit into a TGraph WITH ITS OWN SHAPE RESTORED. Drawing the live TF1s would
   // re-enter sp_model later, under whatever globals happen to be set at paint time.
   auto sample = [&](TF1 *f, bool with2251, int freePk, int color, int style) {
      setShape(with2251, freePk);
      const int N = 1200;
      auto *g = new TGraph(N);
      for (int i = 0; i < N; ++i) {
         double x = exLo + (exHi - exLo) * i / (N - 1.0);
         g->SetPoint(i, x, f->EvalPar(&x, f->GetParameters()));
      }
      g->SetLineColor(color);
      g->SetLineStyle(style);
      g->SetLineWidth(2);
      return g;
   };
   TGraph *gA = sample(fA, true, -1, kRed, 1);
   TGraph *gB = sample(B.first, false, -1, kGray + 2, 3);
   gA->Draw("L SAME");
   gB->Draw("L SAME");
   int col[NL] = {kBlue, kGreen + 2, kMagenta, kOrange + 7};
   TF1 *comp[NL];
   for (int l = 0; l < NL; ++l) {
      // The FULL level -- core PLUS tail -- so the dashed curves really do sum to the red total.
      // A plain "gaus" here would draw only the core, which is a different curve from the one the
      // fit contains and would misread as the component being taller than the total.
      comp[l] = new TF1(Form("c%d", l), "[0]*((1-[3])*exp(-0.5*((x-[1])/[2])^2)+[3]*([2]/[4])*exp(-0.5*((x-[1])/[4])^2))",
                        exLo, exHi);
      comp[l]->SetParameter(0, fA->GetParameter(NSHAPE + l));
      comp[l]->SetParameter(1, LEX[l] + fA->GetParameter(0));
      comp[l]->SetParameter(2, sg);
      comp[l]->SetParameter(3, fA->GetParameter(3));
      comp[l]->SetParameter(4, fA->GetParameter(2));
      comp[l]->SetLineColor(col[l]);
      comp[l]->SetLineStyle(2);
      comp[l]->SetNpx(1200);
      comp[l]->Draw("same");
   }
   TLegend *leg = new TLegend(0.60, 0.58, 0.98, 0.90);
   leg->AddEntry(hObs, Form("sum, %.0f counts", nTot), "lep");
   leg->AddEntry(gA, Form("fit A: 4 levels (#sigma_{core} = %.3f)", sg), "l");
   leg->AddEntry(gB, Form("fit B: no 0^{+}_{2}  (%.1f#sigma worse)", nsig), "l");
   for (int l = 0; l < NL; ++l) leg->AddEntry(comp[l], Form("%s  %.3f MeV  (w %.2f)", LJP[l], LEX[l], pop[l]), "l");
   leg->Draw();
   TString png = TString::Format("%s/spectrum_Be10tp_%s_cm%03.0f-%03.0f_%s.png", outDir.Data(), cfg.Data(), cmLo, cmHi,
                                 vertexEbeam ? "vtxE" : "constE");
   c->SaveAs(png);
   printf("\nwrote %s\n", png.Data());
   printf("spectrum done\n\n");
}
