/// @file dsdo_C17.C
/// @brief What the two MEASURED angular distributions would actually look like: the supplied DWBA
///        curves for the 1/2+ (217 keV) and 5/2+ (332 keV), with the expected data points and
///        error bars this experiment would put on them.
///
///   root -b -q 'dsdo_C17.C'
///   root -b -q 'dsdo_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b400")'
///   root -b -q 'dsdo_C17.C("...","pp_b285", 30, 90)'      // widen the window to show the cost
///
/// The ground state / elastic is left ARBITRARY here (R = N_el/N_217 = 10, flat in angle, as in
/// decompose_C17.C) -- it only enters through the overlap penalty on the two error bars, and no
/// elastic calculation exists to do better. The two inelastic curves are the supplied FRESCO
/// calculations and nothing here rescales them.
///
/// THE FIRST THING THE CURVES SAY, before any detector enters: the two shapes are nearly
/// PROPORTIONAL. dsigma(217)/dsigma(332) runs only 1.8 -> 2.6 across theta_cm 2-150 deg, i.e. the
/// 5/2+ is the 1/2+ divided by ~2.3 with the same structure -- first maximum at theta_cm ~46 deg,
/// diffraction minimum at ~92, secondary maximum at ~132. So the ANGULAR SHAPE carries almost no
/// information separating the two states; what the measurement extracts is two normalisations, and
/// the shape's only job is to confirm the assumed L transfer. That is worth knowing because it
/// means a window that covers one maximum well beats a wider window that smears both.
///
/// theta_cm is the PROJECTILE angle, so theta_lab = (180 - theta_cm)/2 for the light recoil. The
/// first DWBA maximum at theta_cm 46 deg is theta_lab 67 deg -- inside the usable resolution window
/// on the proton day. That is a piece of luck the proposal does not currently claim.
///
/// Points are drawn ON the curve (no Poisson realisation) with the error bars the ledger predicts,
/// which is the standard expected-sensitivity figure: it shows what can be resolved, not one
/// particular outcome. Inner bar = pure statistics, outer bar = with the three-component overlap
/// penalty measured by the toy fit.

#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace
{
double dsOmega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

double dsEx(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * dsOmega2(s, m1 * m1, m2 * m2) * dsOmega2(uu, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                         (2 * m2 * m2) +
                      s + uu - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

bool dsReadFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
{
   std::ifstream in(path);
   if (!in.good())
      return false;
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line[0] == '@')
         continue;
      double a, b;
      if (sscanf(line.c_str(), "%lf %lf", &a, &b) == 2) {
         th.push_back(a);
         xs.push_back(b);
      }
   }
   return th.size() > 10;
}

double dsInterp(const std::vector<double> &th, const std::vector<double> &xs, double x)
{
   if (th.empty() || x < th.front() || x > th.back())
      return 0;
   const size_t i = std::min((size_t)(x - th.front()), th.size() - 2);
   const double f = (x - th[i]) / (th[i + 1] - th[i]);
   return xs[i] * (1 - f) + xs[i + 1] * f;
}

/// sigma integrated over [lo,hi) in theta_cm, and the sigma-weighted mean angle in that range.
double dsSigmaIn(const std::vector<double> &th, const std::vector<double> &xs, double lo, double hi,
                 double *meanTh = nullptr)
{
   double s = 0, sw = 0;
   for (size_t i = 0; i < th.size(); ++i)
      if (th[i] >= lo && th[i] < hi) {
         const double d = xs[i] * 2 * TMath::Pi() * std::sin(th[i] * TMath::DegToRad()) * (TMath::Pi() / 180.0);
         s += d;
         sw += d * th[i];
      }
   if (meanTh)
      *meanTh = (s > 0) ? sw / s : 0.5 * (lo + hi);
   return s;
}

bool dsSolve3(double A[3][3], const double b[3], double a[3])
{
   double M[3][4];
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j)
         M[i][j] = A[i][j];
      M[i][3] = b[i];
   }
   for (int col = 0; col < 3; ++col) {
      int piv = col;
      for (int r = col + 1; r < 3; ++r)
         if (std::fabs(M[r][col]) > std::fabs(M[piv][col]))
            piv = r;
      if (std::fabs(M[piv][col]) < 1e-30)
         return false;
      if (piv != col)
         for (int j = 0; j < 4; ++j)
            std::swap(M[col][j], M[piv][j]);
      const double d = M[col][col];
      for (int j = 0; j < 4; ++j)
         M[col][j] /= d;
      for (int r = 0; r < 3; ++r) {
         if (r == col)
            continue;
         const double f = M[r][col];
         for (int j = 0; j < 4; ++j)
            M[r][j] -= f * M[col][j];
      }
   }
   for (int i = 0; i < 3; ++i)
      a[i] = M[i][3];
   return true;
}

void dsToy(TH1D *tpl[3], const double Ntrue[3], int nToy, double &d217, double &d332)
{
   d217 = d332 = -1;
   const int nB = tpl[0]->GetNbinsX();
   std::vector<double> mu(nB, 0.0);
   for (int b = 1; b <= nB; ++b)
      for (int l = 0; l < 3; ++l)
         mu[b - 1] += Ntrue[l] * tpl[l]->GetBinContent(b);
   double s1 = 0, s11 = 0, s2 = 0, s22 = 0;
   int nOk = 0;
   for (int it = 0; it < nToy; ++it) {
      std::vector<double> n(nB);
      for (int b = 0; b < nB; ++b)
         n[b] = gRandom->Poisson(mu[b]);
      double A[3][3] = {{0}}, rhs[3] = {0, 0, 0};
      for (int b = 0; b < nB; ++b) {
         const double v = std::max(n[b], 1.0);
         double T[3];
         for (int l = 0; l < 3; ++l)
            T[l] = tpl[l]->GetBinContent(b + 1);
         for (int i = 0; i < 3; ++i) {
            rhs[i] += T[i] * n[b] / v;
            for (int j = 0; j < 3; ++j)
               A[i][j] += T[i] * T[j] / v;
         }
      }
      double a[3];
      if (!dsSolve3(A, rhs, a))
         continue;
      ++nOk;
      s1 += a[1];
      s11 += a[1] * a[1];
      s2 += a[2];
      s22 += a[2] * a[2];
   }
   if (nOk < nToy / 4 || Ntrue[1] <= 0 || Ntrue[2] <= 0)
      return;
   const double m1 = s1 / nOk, m2 = s2 / nOk;
   d217 = std::sqrt(std::max(0.0, s11 / nOk - m1 * m1)) / Ntrue[1];
   d332 = std::sqrt(std::max(0.0, s22 / nOk - m2 * m2)) / Ntrue[2];
}
} // namespace

/// @param thLo,thHi  theta_lab window to place points in; <0 = the usable window from
///                   kine_lines_C17.C. Widen it deliberately to see what the extra coverage costs.
/// @param dTh        theta_lab bin width [deg]. 5 gives the finest binning the templates support.
void dsdo_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString cfgTag = "pp_b285",
              Double_t thLo = -1, Double_t thHi = -1, Double_t dTh = 5.0, Double_t duty = 0.70,
              Double_t purity = 1.00, Double_t cleanFrac = 0.930, Double_t days = 1.0, Double_t R = 10.0,
              Int_t nToy = 400, Double_t chi2Cut = 5.0, Double_t zFidLo = 50.0, Double_t zFidHi = 950.0,
              TString frescoDir = "./fresco/", TString outDir = "./plots/")
{
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   gRandom->SetSeed(20260903);

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};

   TString chan(cfgTag);
   chan.Remove(2);
   TString btag(cfgTag);
   btag.Remove(0, 3);

   if (thLo < 0 || thHi < 0) {
      if (cfgTag == "pp_b285") {
         thLo = 50;
         thHi = 70;
      } else if (cfgTag == "pp_b400") {
         thLo = 40;
         thHi = 70;
      } else if (cfgTag == "dd_b285") {
         thLo = 70;
         thHi = 80;
      } else {
         thLo = 60;
         thHi = 80;
      }
   }

   const double uAmu = 931.49401;
   const double mBeam = 17.0225787 * uAmu, mRes = mBeam;
   const double mLight = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu;
   const double EbeamConst = 135.0;
   double ebz_a = 0, ebz_b = 0;
   bool haveEbz = false;

   printf("\n\033[1;36m##########################################################################\033[0m\n");
   printf("\033[1;36m 17C EXPECTED ANGULAR DISTRIBUTIONS -- %s, theta_lab %.0f-%.0f in %.0f deg bins\033[0m\n",
          cfgTag.Data(), thLo, thHi, dTh);
   printf("\033[1;36m##########################################################################\033[0m\n");

   std::vector<double> fth[2], fxs[2];
   if (!dsReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) ||
       !dsReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1])) {
      printf("\033[1;31m  no FRESCO tables in %s. Aborting.\033[0m\n", frescoDir.Data());
      return;
   }
   double sig4pi[2];
   for (int k = 0; k < 2; ++k)
      sig4pi[k] = dsSigmaIn(fth[k], fxs[k], 0, 181);

   // ---- what the curves say on their own -------------------------------------------------------
   printf("\n  \033[1mTHE TWO DWBA SHAPES ARE NEARLY PROPORTIONAL\033[0m\n");
   {
      double rmin = 1e9, rmax = -1e9, pk[2] = {0, 0}, pkv[2] = {-1, -1}, mn = 0, mnv = 1e9;
      for (size_t i = 0; i < fth[0].size(); ++i) {
         if (fth[0][i] < 2 || fth[0][i] > 150)
            continue;
         const double r = fxs[0][i] / std::max(1e-12, fxs[1][i]);
         rmin = std::min(rmin, r);
         rmax = std::max(rmax, r);
      }
      for (int k = 0; k < 2; ++k)
         for (size_t i = 0; i < fth[k].size(); ++i)
            if (fth[k][i] < 100 && fxs[k][i] > pkv[k]) {
               pkv[k] = fxs[k][i];
               pk[k] = fth[k][i];
            }
      for (size_t i = 0; i < fth[0].size(); ++i)
         if (fth[0][i] > 60 && fth[0][i] < 120 && fxs[0][i] < mnv) {
            mnv = fxs[0][i];
            mn = fth[0][i];
         }
      printf("    dsigma(217)/dsigma(332) over theta_cm 2-150 deg:  %.2f to %.2f  (a factor ~%.1f, no shape)\n", rmin,
             rmax, 0.5 * (rmin + rmax));
      printf("    first maximum:  theta_cm %.0f / %.0f deg  ->  theta_lab %.0f / %.0f deg\n", pk[0], pk[1],
             (180 - pk[0]) / 2, (180 - pk[1]) / 2);
      printf("    diffraction minimum: theta_cm %.0f deg  ->  theta_lab %.0f deg\n", mn, (180 - mn) / 2);
      printf("    => the shape does NOT separate the two states. The measurement is two\n");
      printf("       NORMALISATIONS; the shape only confirms the assumed L transfer.\n");
      printf("    window theta_lab %.0f-%.0f  <->  theta_cm %.0f-%.0f: %s the first maximum, %s the minimum\n", thLo,
             thHi, 180 - 2 * thHi, 180 - 2 * thLo, (pk[0] >= 180 - 2 * thHi && pk[0] <= 180 - 2 * thLo) ? "COVERS" : "misses",
             (mn >= 180 - 2 * thHi && mn <= 180 - 2 * thLo) ? "REACHES" : "misses");
   }

   const double nTgt = 3.308e-5 / 1.007 * 6.02214076e23 * 100.0;
   const double lumi = 940.0 * 86400.0 * days * duty * purity * cleanFrac * nTgt * 1e-27; // per mb

   // ---- campaign: per theta_lab bin templates and acceptance -----------------------------------
   const double cmLo = 10.0, cmHi = 178.0;
   const int nAB = 24;
   const double exLo = -1.5, exHi = 2.0;
   const int nB = 70;
   const int nW = std::max(1, (int)std::lround((thHi - thLo) / dTh));
   TH1D *tplW[nL][16];
   long nAccW[nL][16] = {{0}};
   std::vector<double> cmAccW[nL][16];
   double nGen[nL] = {0, 0, 0};

   for (int l = 0; l < nL; ++l) {
      for (int w = 0; w < nW; ++w) {
         tplW[l][w] = new TH1D(Form("ds_%s_%d", stTag[l], w), "", nB, exLo, exHi);
         tplW[l][w]->SetDirectory(nullptr);
      }
      TString dirq = TString("\"") + root + "/" + cfgTag + "\"";
      TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
      TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-10s MISSING (%s)\033[0m\n", stTag[l], pat.Data());
         return;
      }
      TString accLog = dirq + "/" + chan + "_" + stTag[l] + "_" + btag + "_s*_acc.log";
      nGen[l] = atof(gSystem
                        ->GetFromPipe("grep -h 'generated reactions' " + accLog +
                                      " 2>/dev/null | head -1 | awk '{print $3}'")
                        .Strip(TString::kBoth)
                        .Data());
      TFile *f = TFile::Open(found);
      TTree *t = f ? (TTree *)f->Get("res") : nullptr;
      if (!t) {
         printf("\033[1;31m  no res tree in %s\033[0m\n", found.Data());
         return;
      }
      double exR, cmT, c2n, thT, thR, keT, keR, zT, zR;
      t->SetBranchAddress("exReco", &exR);
      t->SetBranchAddress("cmTrue", &cmT);
      t->SetBranchAddress("chi2ndf", &c2n);
      t->SetBranchAddress("thTrue", &thT);
      t->SetBranchAddress("thReco", &thR);
      t->SetBranchAddress("keTrue", &keT);
      t->SetBranchAddress("keReco", &keR);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);

      if (!haveEbz) {
         TGraph g;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = EbeamConst - 40, hi = EbeamConst + 40;
            if ((dsEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                   (dsEx(mBeam, mLight, mLight, mRes, hi, thT, keT) - ExGen[l]) >
                0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((dsEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                      (dsEx(mBeam, mLight, mLight, mRes, mid, thT, keT) - ExGen[l]) <=
                   0)
                  hi = mid;
               else
                  lo = mid;
            }
            g.SetPoint(g.GetN(), zT, 0.5 * (lo + hi));
         }
         if (g.GetN() > 100) {
            TF1 lin("lin", "pol1");
            g.Fit(&lin, "QN");
            ebz_a = lin.GetParameter(0);
            ebz_b = lin.GetParameter(1);
            haveEbz = true;
         }
      }

      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (c2n >= chi2Cut || thR < thLo || thR >= thHi || zR < zFidLo || zR > zFidHi)
            continue;
         if (haveEbz) {
            const double ev = dsEx(mBeam, mLight, mLight, mRes, ebz_a + ebz_b * zR, thR, keR);
            if (ev < -1e8)
               continue;
            exR = ev;
         }
         double w = 1.0;
         if (l == 1)
            w = dsInterp(fth[0], fxs[0], cmT);
         else if (l == 2)
            w = dsInterp(fth[1], fxs[1], cmT);
         if (w <= 0)
            continue;
         const int wb = std::min(nW - 1, (int)((thR - thLo) / (thHi - thLo) * nW));
         if (wb >= 0) {
            tplW[l][wb]->Fill(exR, w);
            cmAccW[l][wb].push_back(cmT);
            ++nAccW[l][wb];
         }
      }
      f->Close();
   }

   // acceptance x cuts, per level, restricted to one theta_lab bin's accepted events
   auto detFrac = [&](int l, const std::vector<double> &cms) {
      if (nGen[l] <= 0)
         return -1.0;
      std::vector<double> accN(nAB, 0.0);
      for (double th : cms) {
         const int b = (int)((th - cmLo) / (cmHi - cmLo) * nAB);
         if (b >= 0 && b < nAB)
            accN[b] += 1.0;
      }
      const double dcosTot = std::cos(cmLo * TMath::DegToRad()) - std::cos(cmHi * TMath::DegToRad());
      const std::vector<double> &fth_l = (l == 2) ? fth[1] : fth[0];
      const std::vector<double> &fxs_l = (l == 2) ? fxs[1] : fxs[0];
      const double sig_l = (l == 2) ? sig4pi[1] : sig4pi[0];
      double d = 0;
      for (int b = 0; b < nAB; ++b) {
         const double lo = cmLo + (cmHi - cmLo) * b / nAB;
         const double hi = cmLo + (cmHi - cmLo) * (b + 1) / nAB;
         const double dcos = std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad());
         const double gen = nGen[l] * dcos / dcosTot;
         if (gen <= 0)
            continue;
         d += (dsSigmaIn(fth_l, fxs_l, lo, hi) / sig_l) * (accN[b] / gen);
      }
      return d;
   };

   // ---- the points -----------------------------------------------------------------------------
   printf("\n  \033[1mEXPECTED DATA POINTS  (%.1f d, duty %.2f, R = %.0f)\033[0m\n", days, duty, R);
   printf("    %-13s %-13s %8s %8s %9s %9s %9s %9s\n", "theta_lab", "theta_cm", "N(217)", "N(332)", "ds217", "+-",
          "ds332", "+-");
   printf("    %-13s %-13s %8s %8s %9s %9s %9s %9s\n", "[deg]", "[deg]", "", "", "[mb/sr]", "[%]", "[mb/sr]", "[%]");

   std::vector<double> x, xe, y1, ye1, y1s, y2, ye2, y2s;
   for (int w = 0; w < nW; ++w) {
      const double wlo = thLo + (thHi - thLo) * w / nW, whi = thLo + (thHi - thLo) * (w + 1) / nW;
      const double cLo = 180 - 2 * whi, cHi = 180 - 2 * wlo;
      double mth1 = 0, mth2 = 0;
      const double s1 = dsSigmaIn(fth[0], fxs[0], cLo, cHi, &mth1);
      const double s2 = dsSigmaIn(fth[1], fxs[1], cLo, cHi, &mth2);
      const double dOm = 2 * TMath::Pi() *
                         std::fabs(std::cos(cLo * TMath::DegToRad()) - std::cos(cHi * TMath::DegToRad()));
      double N[nL];
      for (int l = 1; l < nL; ++l)
         N[l] = lumi * sig4pi[l - 1] * detFrac(l, cmAccW[l][w]);
      N[0] = R * lumi * sig4pi[0] * detFrac(0, cmAccW[0][w]);
      if (N[1] <= 0 || N[2] <= 0 || nAccW[1][w] < 50 || nAccW[2][w] < 50) {
         printf("    %5.0f-%-7.0f %5.0f-%-7.0f %8.0f %8.0f   template too thin\n", wlo, whi, cLo, cHi,
                std::max(0.0, N[1]), std::max(0.0, N[2]));
         continue;
      }
      TH1D *shw[nL];
      for (int l = 0; l < nL; ++l) {
         shw[l] = (TH1D *)tplW[l][w]->Clone(Form("shw%d_%d", l, w));
         shw[l]->SetDirectory(nullptr);
         if (shw[l]->Integral() > 0)
            shw[l]->Scale(1.0 / shw[l]->Integral());
      }
      double e1, e2;
      dsToy(shw, N, nToy, e1, e2);
      for (int l = 0; l < nL; ++l)
         delete shw[l];
      if (e1 <= 0)
         e1 = 1.0 / std::sqrt(N[1]);
      if (e2 <= 0)
         e2 = 1.0 / std::sqrt(N[2]);
      const double d1 = s1 / dOm, d2 = s2 / dOm;
      printf("    %5.0f-%-7.0f %5.0f-%-7.0f %8.0f %8.0f %9.4f %9.1f %9.4f %9.1f\n", wlo, whi, cLo, cHi, N[1], N[2], d1,
             100 * e1, d2, 100 * e2);
      x.push_back(0.5 * (mth1 + mth2));
      xe.push_back(0.5 * (cHi - cLo));
      y1.push_back(d1);
      ye1.push_back(d1 * e1);
      y1s.push_back(d1 / std::sqrt(N[1]));
      y2.push_back(d2);
      ye2.push_back(d2 * e2);
      y2s.push_back(d2 / std::sqrt(N[2]));
   }
   if (x.size() < 2) {
      printf("\033[1;31m  fewer than two usable points -- nothing to draw\033[0m\n");
      return;
   }

   // ---- the figure -----------------------------------------------------------------------------
   gSystem->mkdir(outDir, kTRUE);
   TCanvas c("cds", "", 1250, 620);
   c.Divide(2, 1);

   // -- panel 1: the two distributions, curve + expected points
   c.cd(1)->SetLogy();
   c.cd(1)->SetGridx();
   c.cd(1)->SetGridy();
   TGraph *g217 = new TGraph(), *g332 = new TGraph();
   for (size_t i = 0; i < fth[0].size(); ++i) {
      g217->SetPoint(g217->GetN(), fth[0][i], fxs[0][i]);
      g332->SetPoint(g332->GetN(), fth[1][i], fxs[1][i]);
   }
   g217->SetLineColor(kAzure + 2);
   g217->SetLineWidth(3);
   g332->SetLineColor(kRed + 1);
   g332->SetLineWidth(3);
   g217->SetTitle(Form("^{17}C(%s,%s') DWBA + expected data, %s, %.1f d;#theta_{cm} [deg]  "
                       "(#theta_{lab} = (180-#theta_{cm})/2);d#sigma/d#Omega [mb/sr]",
                       chan == "dd" ? "d" : "p", chan == "dd" ? "d" : "p", cfgTag.Data(), days));
   g217->GetXaxis()->SetLimits(0, 180);
   g217->SetMinimum(0.008);
   g217->SetMaximum(4.0);
   g217->Draw("AL");
   // the measurable window as a shaded band
   TBox *bx = new TBox(180 - 2 * thHi, 0.008, 180 - 2 * thLo, 4.0);
   bx->SetFillColorAlpha(kGreen - 9, 0.30);
   bx->Draw();
   g217->Draw("L");
   g332->Draw("L");
   TGraphErrors *p217 = new TGraphErrors(x.size(), &x[0], &y1[0], &xe[0], &ye1[0]);
   TGraphErrors *p332 = new TGraphErrors(x.size(), &x[0], &y2[0], &xe[0], &ye2[0]);
   TGraphErrors *q217 = new TGraphErrors(x.size(), &x[0], &y1[0], nullptr, &y1s[0]);
   TGraphErrors *q332 = new TGraphErrors(x.size(), &x[0], &y2[0], nullptr, &y2s[0]);
   for (auto *g : {p217, q217}) {
      g->SetMarkerStyle(20);
      g->SetMarkerSize(1.3);
      g->SetMarkerColor(kAzure + 2);
      g->SetLineColor(kAzure + 2);
   }
   for (auto *g : {p332, q332}) {
      g->SetMarkerStyle(21);
      g->SetMarkerSize(1.3);
      g->SetMarkerColor(kRed + 1);
      g->SetLineColor(kRed + 1);
   }
   p217->SetLineWidth(2);
   p332->SetLineWidth(2);
   q217->SetLineWidth(5);
   q332->SetLineWidth(5);
   q217->Draw("PZ");
   q332->Draw("PZ");
   p217->Draw("PZ");
   p332->Draw("PZ");
   TLegend *lg = new TLegend(0.40, 0.14, 0.89, 0.36);
   lg->SetFillColorAlpha(kWhite, 0.85);
   lg->AddEntry(g217, "1/2^{+} 217 keV  (FRESCO, 8.58 mb)", "l");
   lg->AddEntry(g332, "5/2^{+} 332 keV  (FRESCO, 3.60 mb)", "l");
   lg->AddEntry(p217, Form("expected points, %.0f%% duty, R = %.0f", 100 * duty, R), "pe");
   lg->AddEntry(bx, "usable resolution window", "f");
   lg->Draw();

   // -- panel 2: the ratio, i.e. is there any shape information at all
   c.cd(2);
   c.cd(2)->SetGridx();
   c.cd(2)->SetGridy();
   TGraph *gr = new TGraph();
   for (size_t i = 0; i < fth[0].size(); ++i)
      if (fxs[1][i] > 1e-9 && fth[0][i] <= 150)
         gr->SetPoint(gr->GetN(), fth[0][i], fxs[0][i] / fxs[1][i]);
   gr->SetLineColor(kBlack);
   gr->SetLineWidth(3);
   gr->SetTitle("ratio of the two DWBA curves;#theta_{cm} [deg];d#sigma(217) / d#sigma(332)");
   gr->GetXaxis()->SetLimits(0, 180);
   gr->SetMinimum(0);
   gr->SetMaximum(5);
   gr->Draw("AL");
   TBox *bx2 = new TBox(180 - 2 * thHi, 0, 180 - 2 * thLo, 5);
   bx2->SetFillColorAlpha(kGreen - 9, 0.30);
   bx2->Draw();
   gr->Draw("L");
   std::vector<double> rv, rve;
   for (size_t i = 0; i < x.size(); ++i) {
      rv.push_back(y1[i] / y2[i]);
      rve.push_back((y1[i] / y2[i]) *
                    std::sqrt(std::pow(ye1[i] / y1[i], 2) + std::pow(ye2[i] / y2[i], 2)));
   }
   TGraphErrors *pr = new TGraphErrors(x.size(), &x[0], &rv[0], &xe[0], &rve[0]);
   pr->SetMarkerStyle(20);
   pr->SetMarkerSize(1.3);
   pr->SetLineWidth(2);
   pr->Draw("PZ");
   TLatex tx;
   tx.SetNDC();
   tx.SetTextSize(0.030);
   tx.DrawLatex(0.14, 0.86, "flat to #pm15% over #theta_{cm} 2-150#circ:");
   tx.DrawLatex(0.14, 0.82, "the two shapes carry no separating power.");
   tx.DrawLatex(0.14, 0.78, "The measurement is two NORMALISATIONS.");

   TString stem = TString::Format("dsdo_%s_%.0f-%.0f_d%.0f", cfgTag.Data(), thLo, thHi, dTh);
   c.SaveAs(outDir + stem + ".png");
   printf("\n    wrote %s%s.png\n\n", outDir.Data(), stem.Data());
}
