/// @file tp_common_Be10.h
/// @brief Shared kinematics, level table and file lookup for the 10Be(t,p)12Be plot macros.
///
/// Factored out because the (d,p) campaign carried four private copies of the same two-body
/// expressions across its plot macros, and a copy is a place for them to drift apart.

#ifndef TP_COMMON_BE10_H
#define TP_COMMON_BE10_H

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TKey.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace tpc {

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const char *CFGL[NC] = {"2.85 T, AT-TPC", "2.85 T, 2 mm", "4 T, AT-TPC",
                               "4 T, 2 mm",      "7 T, AT-TPC",  "7 T, 2 mm"};
static const int CCOL[NC] = {kBlack, kGray + 2, kBlue, kAzure + 7, kRed + 1, kOrange + 7};
static const int CSTY[NC] = {1, 2, 1, 2, 1, 2};

static const int NL = 4;
static const char *LVN[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
static const double LEX[NL] = {0.0, 2.109, 2.251, 2.715};
static const char *LJP[NL] = {"0^{+} g.s.", "2^{+} 2.109", "0^{+}_{2} 2.251", "1^{-} 2.715"};
static const int LCOL[NL] = {kBlue, kGreen + 2, kMagenta, kOrange + 7};

// 10Be(t,p)12Be, and 14C(d,p)15C for the comparison panels
static const double U = 931.49401;
static const double M1 = 10.0135341 * U, M2 = 3.0160493 * U, M3 = 1.007825 * U, M4 = 12.0269221 * U;
static const double EB = 112.20;

/// THE WINDOW THIS CHANNEL IS ANALYSED IN. theta_cm 2-45 is the transfer peak, i.e. the backward
/// lab proton. Past ~45 the proton is fast and forward and sigma(KE)/KE runs 3-6 %, which no
/// correction reaches -- see the error-budget comment in kin_gui_Be10.C.
static const double PEAK_LO = 2.0, PEAK_HI = 45.0;

inline double om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
/// Ex from (theta_lab [rad], KE) at beam energy K, for residual ground-state mass m4
inline double ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + M1, Et3 = Ke + M3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double a = (std::cos(th) * om2(s, M1 * M1, M2 * M2) * om2(uu, M2 * M2, M3 * M3) -
               (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) / (2 * M2 * M2) + s + uu - M2 * M2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
/// theta_cm (DWBA convention) -> (theta_lab [rad], KE) for residual mass m4
inline bool fwd(double m4, double thcmA, double &thlab, double &Ke)
{
   const double thcm = TMath::Pi() - thcmA;
   double E1 = EB + M1, s = M1 * M1 + M2 * M2 + 2 * M2 * E1, rs = std::sqrt(s);
   if (rs < M3 + m4) return false;
   double pcm = om2(s, M3 * M3, m4 * m4) / (2 * rs), Ecm3 = std::sqrt(pcm * pcm + M3 * M3);
   double plab = std::sqrt(E1 * E1 - M1 * M1), beta = plab / (E1 + M2), gam = 1 / std::sqrt(1 - beta * beta);
   Ke = gam * (Ecm3 + beta * pcm * std::cos(thcm)) - M3;
   thlab = std::atan2(pcm * std::sin(thcm), gam * (pcm * std::cos(thcm) + beta * Ecm3));
   return true;
}
inline double quant(std::vector<double> v, double p)
{
   if (v.size() < 20) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
inline TString find(const TString &dir, const char *cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg, pre, lvl, cfg));
   return f.Strip(TString::kBoth);
}
/// The E_beam(z) profile solved from truth, as tp_spectrum_Be10.C and the GUI both do. Returns
/// false if there is not enough to fit, so a caller can fall back rather than use a wild profile.
inline bool ebeamProfile(TTree *t, double lvlEx, TF1 &fEb)
{
   double thT, keT, zT;
   t->SetBranchAddress("thTrue", &thT);
   t->SetBranchAddress("keTrue", &keT);
   t->SetBranchAddress("zTrue", &zT);
   std::vector<double> eb, zz;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      double lo = 80., hi = 130., thr = thT * TMath::DegToRad();
      double flo = ex(M4, lo, thr, keT) - lvlEx, fhi = ex(M4, hi, thr, keT) - lvlEx;
      if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
      for (int it = 0; it < 60; ++it) {
         double m = 0.5 * (lo + hi), fm = ex(M4, m, thr, keT) - lvlEx;
         if (std::isnan(fm)) break;
         if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
      }
      double e = 0.5 * (lo + hi);
      if (e > 85 && e < 125) { eb.push_back(e); zz.push_back(zT); }
   }
   fEb.SetParameters(EB, -0.006, 0.);
   if (eb.size() < 100) return false;
   TGraph g((int)eb.size(), zz.data(), eb.data());
   g.Fit(&fEb, "QN");
   return true;
}
/// gen/rec histograms out of an acceptance file
inline bool acceptanceHists(const TString &dir, const char *cfg, const char *lvl, TH1D *&hg, TH1D *&hr, TFile *&f)
{
   hg = nullptr; hr = nullptr;
   TString fa = find(dir, cfg, lvl, "acceptance");
   if (fa.IsNull()) return false;
   f = TFile::Open(fa);
   if (!f || f->IsZombie()) return false;
   TIter nx(f->GetListOfKeys());
   while (auto *k = (TKey *)nx()) {
      TString n = k->GetName();
      if (n.BeginsWith("hGen_")) hg = (TH1D *)f->Get(n);
      if (n.BeginsWith("hRec_")) hr = (TH1D *)f->Get(n);
   }
   return hg && hr;
}

} // namespace tpc
#endif
