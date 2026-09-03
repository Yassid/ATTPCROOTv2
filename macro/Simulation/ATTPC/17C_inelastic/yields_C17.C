/// @file yields_C17.C
/// @brief The COUNT LEDGER for the 17C M_n/M_p proposal: what does "~1000 counts" actually become
///        once every factor between a 4pi cross section and a fitted differential point is written
///        down and either measured from the campaign or named as a parameter?
///
///   root -b -q 'yields_C17.C'
///   root -b -q 'yields_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b285")'
///   root -b -q 'yields_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","dd_b400",-1,-1,0.8)'
///
/// WHY THIS EXISTS. The proposal's expected statistics -- "~1000 counts of scattered protons and a
/// similar amount for deuterons" (p.4) -- is not derived from a cross section at all. It is an
/// analogy to Ref.[24], the 14C(p,p') work, scaled by a factor 2 in beam rate. decompose_C17.C
/// replaced that analogy with arithmetic from the supplied FRESCO calculations and got 1379 (1/2+)
/// and 579 (5/2+) at 4pi, 1163 and 487 with acceptance folded, which happens to bracket the quoted
/// number. But those two numbers are still NOT what the experiment records, for three reasons that
/// all push the same way:
///
///   1. ONE DAY OF BEAM TIME IS NOT 86400 s OF BEAM ON TARGET. Every yield in this study, and in
///      the proposal, multiplies 940 pps by a full 86400 s. Nothing in the simulation can know the
///      duty factor, so it is a named parameter here with an explicitly unsourced default. It is
///      the single largest correction in the ledger and the cheapest to get right -- it needs a
///      number from ReA6 operations, not a calculation.
///   2. THE COUNTS OUTSIDE THE USABLE ANGULAR WINDOW ARE NOT COUNTS. inel_summary_C17.C and
///      kine_lines_C17.C established that sigma(Ex) varies 20x across theta_lab and that only
///      50-70 deg (p,p' at 2.85 T), 40-70 (p,p' at 4 T), 70-80 (d,d' at 2.85 T) and 60-80 (d,d' at
///      4 T) are usable. The FRESCO distributions put only 43 % of the (p,p') yield in the 2.85 T
///      window and 68 % in the 4 T one. The "1163 and 487 detected" figure integrates over angles
///      whose Ex spectrum is a featureless bump, so it counts events that cannot enter the fit.
///   3. THE DELIVERABLE IS A DIFFERENTIAL CROSS SECTION, NOT A YIELD. M_n/M_p comes from the
///      deformation lengths, i.e. from the normalisation of a dsigma/dOmega CURVE. The proposal's
///      ~10 % is quoted per state; what has to hold is ~10 % per angular point, which is another
///      sqrt(n_bins) worse. This macro reports both.
///
/// So this is deliberately not a new simulation. Every detector effect here is already measured by
/// the existing campaign; what was missing was one place where the chain from mb to fitted points
/// is written out step by step, with each factor either MEASURED (marked [sim]) or NAMED (marked
/// [par], with its default's provenance). A factor that is a parameter is not a weakness as long
/// as it is visible; the failure mode this macro exists to prevent is a factor that is silently 1.
///
/// WHAT IS STILL NOT IN THE LEDGER, and cannot be until something is added upstream:
///   - The ELASTIC angular distribution. R = N_elastic/N_217 is still scanned, and per angular bin
///     it is assumed constant, which is certainly wrong -- elastic is forward-peaked, so the real R
///     rises steeply towards small theta_lab, i.e. towards the forward end of every usable window.
///     Getting the elastic curve remains the highest-value missing input, and the per-bin table
///     below is where it would land.
///   - The TRIGGER. The AT-TPC trigger counts PADS above threshold, and a backward-angle recoil in
///     a 2.85-4 T field is a tightly curled helix that revisits the same pads. The digitised events
///     exist in the campaign's sim files; nothing has yet counted their pad multiplicity against a
///     threshold. This bites hardest in exactly the backward window the resolution study selected.
///   - COMPETING CHANNELS on the deuterium day ((d,p), (d,t), (d,3He), breakup) and the 17O beam
///     contaminant, neither of which exists in the simulation at any rate.
///
/// METHOD. Identical machinery to decompose_C17.C -- same trees, same vertex correction, same
/// flat-generator acceptance folding, same fixed-position three-component toy fit -- run once
/// globally and then once per angular bin inside the usable window, with the ledger factors
/// applied to the normalisation.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
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
double ylOmega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

/// Two-body inversion: measured ejectile (theta_lab, KE) plus a beam energy -> residual excitation
/// energy. Same form as decompose_C17.C:dcEx and inel_summary_C17.C:s17_ex.
double ylEx(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * ylOmega2(s, m1 * m1, m2 * m2) * ylOmega2(uu, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                         (2 * m2 * m2) +
                      s + uu - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

bool ylReadFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
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

double ylInterp(const std::vector<double> &th, const std::vector<double> &xs, double x)
{
   if (th.empty() || x < th.front() || x > th.back())
      return 0;
   const size_t i = std::min((size_t)(x - th.front()), th.size() - 2);
   const double f = (x - th[i]) / (th[i + 1] - th[i]);
   return xs[i] * (1 - f) + xs[i + 1] * f;
}

/// dsigma/dOmega integrated over [lo,hi] in theta_cm, on the table's own 1 deg grid.
double ylSigmaIn(const std::vector<double> &th, const std::vector<double> &xs, double lo, double hi)
{
   double s = 0;
   for (size_t i = 0; i < th.size(); ++i)
      if (th[i] >= lo && th[i] < hi)
         s += xs[i] * 2 * TMath::Pi() * std::sin(th[i] * TMath::DegToRad()) * (TMath::Pi() / 180.0);
   return s;
}

bool ylSolve3(double A[3][3], const double b[3], double a[3])
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

/// One fixed-position three-component toy experiment set. Returns the fractional width of the
/// recovered 1/2+ and 5/2+ yields, or (-1,-1) if the templates are too thin to fit.
void ylToy(TH1D *tpl[3], const double Ntrue[3], int nToy, double &d217, double &d332)
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
      if (!ylSolve3(A, rhs, a))
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

/// @param root       campaign directory (INEL_ROOT)
/// @param cfgTag     pp_b285 | pp_b400 | dd_b285 | dd_b400
/// @param thLo,thHi  usable theta_lab window; <0 takes the per-configuration default established by
///                   kine_lines_C17.C (pp 2.85 T: 50-70, pp 4 T: 40-70, dd 2.85 T: 70-80,
///                   dd 4 T: 60-80)
/// @param duty       [par] beam-on-target fraction of the requested day. DEFAULT IS A GUESS: no
///                   ReA6 number was supplied with the proposal. Set it from operations experience.
/// @param purity     [par] 17C fraction of the delivered cocktail. The proposal states the only
///                   expected contaminant is 17O and that its amount "was not able to be
///                   estimated", so the default is 1.0 -- an upper limit, not a measurement.
/// @param cleanFrac  [par] fraction of events free of beam pile-up inside one drift window.
///                   Default is computed, not guessed: 100 cm at 1.30 cm/us = 76.9 us, 940 pps
///                   -> exp(-0.072) = 0.930. Pile-up may be recoverable in analysis rather than
///                   lost, so this is the pessimistic end.
/// @param days       number of days on this gas.
/// @param R          [par] N_elastic/N_217, still unmeasured. Assumed CONSTANT in angle, which is
///                   wrong in a known direction: elastic is forward-peaked.
void yields_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString cfgTag = "pp_b285",
                Double_t thLo = -1, Double_t thHi = -1, Double_t duty = 0.70, Double_t purity = 1.00,
                Double_t cleanFrac = 0.930, Double_t days = 1.0, Double_t R = 10.0, Int_t nToy = 400,
                Double_t chi2Cut = 5.0, Double_t zFidLo = 50.0, Double_t zFidHi = 950.0,
                TString frescoDir = "./fresco/", TString outDir = "./plots/")
{
   gStyle->SetOptStat(0);
   gRandom->SetSeed(20260903);

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};

   TString chan(cfgTag);
   chan.Remove(2);
   TString btag(cfgTag);
   btag.Remove(0, 3);

   // Usable windows, from kine_lines_C17.C's per-window separation table.
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
   const double mBeam = 17.0225787 * uAmu;
   const double mRes = mBeam;
   const double mLight = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu;
   const double EbeamConst = 135.0;
   double ebz_a = 0, ebz_b = 0;
   bool haveEbz = false;

   printf("\n\033[1;36m##########################################################################\033[0m\n");
   printf("\033[1;36m 17C COUNT LEDGER -- %s, theta_lab %.0f-%.0f, %.1f day(s)\033[0m\n", cfgTag.Data(), thLo, thHi,
          days);
   printf("\033[1;36m##########################################################################\033[0m\n");

   // ---- cross sections -------------------------------------------------------------------------
   std::vector<double> fth[2], fxs[2];
   if (!ylReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) ||
       !ylReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1])) {
      printf("\033[1;31m  no FRESCO tables in %s -- cannot normalise. Aborting.\033[0m\n", frescoDir.Data());
      return;
   }
   double sig4pi[2];
   for (int k = 0; k < 2; ++k)
      sig4pi[k] = ylSigmaIn(fth[k], fxs[k], 0, 181);

   // 300 torr H2 (media.geo H_300torr_RT: atomic H, A = 1.007, rho = 3.308e-5) over 1 m; D2 at the
   // same pressure has the same molecular number density to 0.8 %.
   const double nTgt = 3.308e-5 / 1.007 * 6.02214076e23 * 100.0;
   const double beamRaw = 940.0 * 86400.0 * days;
   const double beamEff = beamRaw * duty * purity;

   printf("\n  \033[1mBEAM AND TARGET\033[0m\n");
   printf("    %-46s %12.4g\n", "target areal density [nuclei/cm2]  [sim]", nTgt);
   printf("    %-46s %12.4g\n", "beam particles requested (940 pps x day)", beamRaw);
   printf("    x %-44s %12.3f\n", "beam-on-target duty            [par, GUESS]", duty);
   printf("    x %-44s %12.3f\n", "17C purity of the cocktail     [par, LIMIT]", purity);
   printf("    %-46s %12.4g\n", "= 17C on target", beamEff);

   if (duty >= 0.999)
      printf("    \033[1;33mduty = 1 means the ledger reproduces the proposal's own assumption, which is that a\n"
             "    requested day delivers 86400 s of beam. Nothing in this study supports that.\033[0m\n");

   // ---- the campaign: acceptance inside the window, per theta_cm bin ---------------------------
   // The generator is FLAT IN cos(theta_cm), so the number generated in a theta_cm bin is
   // nGen x Dcos(bin)/Dcos(range) and the acceptance in that bin is accepted/generated. Because the
   // theta_lab window, the chi2 cut and the z fiducial are all applied in the event loop below,
   // that single measured ratio carries all of them -- no factor is applied twice.
   const double cmLo = 10.0, cmHi = 178.0;
   const int nAB = 24;
   const double exLo = -1.5, exHi = 2.0;
   const int nB = 70; // 50 keV bins; coarser than decompose_C17.C because the per-bin templates
                      // are built from a fraction of the sample

   // per-angular-bin templates: bins of the usable theta_lab window
   const int nW = std::max(1, (int)std::lround((thHi - thLo) / 5.0));
   TH1D *tpl[nL];      // global, window only
   TH1D *tplW[nL][12]; // per 5 deg bin
   long nAccW[nL][12] = {{0}};
   std::vector<double> cmAcc[nL];
   std::vector<double> cmAccW[nL][12];
   double nGen[nL] = {0, 0, 0};
   long nAcc[nL] = {0, 0, 0}, nRead[nL] = {0, 0, 0};

   for (int l = 0; l < nL; ++l) {
      tpl[l] = new TH1D(Form("yl_%s", stTag[l]), "", nB, exLo, exHi);
      tpl[l]->SetDirectory(nullptr);
      for (int w = 0; w < nW; ++w) {
         tplW[l][w] = new TH1D(Form("ylw_%s_%d", stTag[l], w), "", nB, exLo, exHi);
         tplW[l][w]->SetDirectory(nullptr);
      }
      // Quote ONLY the directory -- quoting the whole pattern stops the "s*" wildcard expanding.
      TString dirq = TString("\"") + root + "/" + cfgTag + "\"";
      TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
      TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-10s MISSING (%s)\033[0m\n", stTag[l], pat.Data());
         return;
      }
      TString accLog = dirq + "/" + chan + "_" + stTag[l] + "_" + btag + "_s*_acc.log";
      TString ng = gSystem->GetFromPipe("grep -h 'generated reactions' " + accLog +
                                        " 2>/dev/null | head -1 | awk '{print $3}'");
      nGen[l] = atof(ng.Strip(TString::kBoth).Data());

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
      nRead[l] = t->GetEntries();

      // E_beam(z) from truth, on the first level present. The tree's exReco is the constant-Ebeam
      // reconstruction; the vertex correction has to be REDONE from the tracked (theta, KE).
      if (!haveEbz) {
         TGraph g;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = EbeamConst - 40, hi = EbeamConst + 40;
            if ((ylEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                   (ylEx(mBeam, mLight, mLight, mRes, hi, thT, keT) - ExGen[l]) >
                0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((ylEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                      (ylEx(mBeam, mLight, mLight, mRes, mid, thT, keT) - ExGen[l]) <=
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
         if (c2n >= chi2Cut)
            continue;
         if (thR < thLo || thR >= thHi)
            continue;
         if (zR < zFidLo || zR > zFidHi)
            continue;
         if (haveEbz) {
            const double ev = ylEx(mBeam, mLight, mLight, mRes, ebz_a + ebz_b * zR, thR, keR);
            if (ev < -1e8)
               continue;
            exR = ev;
         }
         // Reweight the flat generator to the real angular distribution. The elastic has no
         // calculation, so it stays flat.
         double w = 1.0;
         if (l == 1)
            w = ylInterp(fth[0], fxs[0], cmT);
         else if (l == 2)
            w = ylInterp(fth[1], fxs[1], cmT);
         if (w <= 0)
            continue;
         tpl[l]->Fill(exR, w);
         cmAcc[l].push_back(cmT);
         ++nAcc[l];
         const int wb = std::min(nW - 1, (int)((thR - thLo) / (thHi - thLo) * nW));
         if (wb >= 0) {
            tplW[l][wb]->Fill(exR, w);
            cmAccW[l][wb].push_back(cmT);
            ++nAccW[l][wb];
         }
      }
      f->Close();
      if (tpl[l]->Integral() <= 0) {
         printf("\033[1;31m  %-10s no accepted events in the window\033[0m\n", stTag[l]);
         return;
      }
   }
   if (haveEbz)
      printf("\n    E_beam(z) = %.3f %+.5f z[mm]   [sim, from truth]\n", ebz_a, ebz_b);

   // detected fraction of the 4pi yield, FRESCO-weighted acceptance, window+cuts folded
   auto detFracOf = [&](int l, const std::vector<double> &cms) {
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
         d += (ylSigmaIn(fth_l, fxs_l, lo, hi) / sig_l) * (accN[b] / gen);
      }
      return d;
   };

   // ---- the ledger -----------------------------------------------------------------------------
   // theta_cm of the window edges: theta_cm = 180 - 2 theta_lab (theta_cm is the PROJECTILE angle).
   const double cmWinLo = std::max(0.0, 180.0 - 2 * thHi), cmWinHi = std::min(180.0, 180.0 - 2 * thLo);
   const double fWin[2] = {ylSigmaIn(fth[0], fxs[0], cmWinLo, cmWinHi) / sig4pi[0],
                           ylSigmaIn(fth[1], fxs[1], cmWinLo, cmWinHi) / sig4pi[1]};
   const double detF[nL] = {detFracOf(0, cmAcc[0]), detFracOf(1, cmAcc[1]), detFracOf(2, cmAcc[2])};

   double Nfit[nL];
   printf("\n  \033[1mFROM CROSS SECTION TO ANALYSABLE COUNTS\033[0m\n");
   printf("    %-46s %12s %12s\n", "", "1/2+ 217", "5/2+ 332");
   printf("    %-46s %12.3f %12.3f\n", "sigma(4pi) [mb]                    [FRESCO]", sig4pi[0], sig4pi[1]);
   const double N4pi[2] = {beamRaw * nTgt * sig4pi[0] * 1e-27, beamRaw * nTgt * sig4pi[1] * 1e-27};
   printf("    %-46s %12.0f %12.0f\n", "x luminosity, day at face value", N4pi[0], N4pi[1]);
   printf("    x %-44s %12.0f %12.0f\n", "duty x purity                        [par]", N4pi[0] * duty * purity,
          N4pi[1] * duty * purity);
   printf("    x %-44s %12.0f %12.0f\n", "pile-up-free fraction                [par]",
          N4pi[0] * duty * purity * cleanFrac, N4pi[1] * duty * purity * cleanFrac);
   printf("    %-46s %12.3f %12.3f\n", "  (FRESCO fraction inside the window)", fWin[0], fWin[1]);
   printf("    x %-44s %12.0f %12.0f\n", "window x acceptance x chi2 x zfid    [sim]",
          N4pi[0] * duty * purity * cleanFrac * detF[1], N4pi[1] * duty * purity * cleanFrac * detF[2]);
   for (int l = 1; l < nL; ++l)
      Nfit[l] = N4pi[l - 1] * duty * purity * cleanFrac * detF[l];
   Nfit[0] = R * N4pi[0] * duty * purity * cleanFrac * detF[0];
   printf("    \033[1m%-46s %12.0f %12.0f\033[0m\n", "= COUNTS THAT ENTER THE FIT", Nfit[1], Nfit[2]);
   printf("    %-46s %12.0f\n", "  elastic under them, at R x N(217) [par]", Nfit[0]);
   printf("\n    for reference, the numbers this replaces:\n");
   printf("      proposal, by analogy with Ref.[24]                ~1000        ~1000\n");
   printf("      decompose_C17.C, 4pi x acc, ALL angles (pp285)     1163          487\n");
   printf("      same window, without duty/purity/pile-up      %8.0f %12.0f\n", N4pi[0] * detF[1],
          N4pi[1] * detF[2]);
   printf("      this ledger, in the usable window             %8.0f %12.0f   (x %.2f, x %.2f\n", Nfit[1], Nfit[2],
          Nfit[1] / 1000.0, Nfit[2] / 1000.0);
   printf("                                                                             of the quoted ~1000)\n");

   // ---- the fit, globally over the window ------------------------------------------------------
   TH1D *sh[nL];
   for (int l = 0; l < nL; ++l) {
      sh[l] = (TH1D *)tpl[l]->Clone(Form("sh_%d", l));
      sh[l]->SetDirectory(nullptr);
      sh[l]->Scale(1.0 / sh[l]->Integral());
   }
   double d217, d332;
   ylToy(sh, Nfit, nToy, d217, d332);
   printf("\n  \033[1mYIELD ERRORS, fixed-position three-component fit, R = %.0f\033[0m\n", R);
   printf("    %-30s %12s %12s\n", "", "d(217)/N [%]", "d(332)/N [%]");
   printf("    %-30s %12.1f %12.1f\n", "pure-statistics floor 1/sqrt(N)", 100 / std::sqrt(Nfit[1]),
          100 / std::sqrt(Nfit[2]));
   printf("    %-30s %12.1f %12.1f\n", "window-integrated, with overlap", 100 * d217, 100 * d332);
   if (d217 > 0)
      printf("    %-30s %12.1f %12.1f\n", "  overlap penalty (x floor)", d217 * std::sqrt(Nfit[1]),
             d332 * std::sqrt(Nfit[2]));

   // ---- per angular bin: the actual deliverable ------------------------------------------------
   printf("\n  \033[1mPER ANGULAR POINT -- what the differential cross section actually costs\033[0m\n");
   printf("    M_n/M_p comes from the normalisation of a dsigma/dOmega CURVE, so this table and not\n");
   printf("    the one above is the figure of merit the proposal's ~10 %% has to be compared with.\n\n");
   printf("    %-14s %8s %8s %10s %10s %10s %10s\n", "theta_lab", "N(217)", "N(332)", "floor217", "floor332",
          "d217 [%]", "d332 [%]");
   double best217 = 1e9, best332 = 1e9;
   std::vector<double> gx, gy217, gy332;
   for (int w = 0; w < nW; ++w) {
      const double wlo = thLo + (thHi - thLo) * w / nW, whi = thLo + (thHi - thLo) * (w + 1) / nW;
      double Nb[nL];
      for (int l = 1; l < nL; ++l)
         Nb[l] = N4pi[l - 1] * duty * purity * cleanFrac * detFracOf(l, cmAccW[l][w]);
      Nb[0] = R * N4pi[0] * duty * purity * cleanFrac * detFracOf(0, cmAccW[0][w]);
      if (Nb[1] <= 0 || Nb[2] <= 0 || nAccW[1][w] < 50 || nAccW[2][w] < 50) {
         printf("    %5.0f-%-8.0f %8.0f %8.0f %10s %10s %10s %10s   (template too thin)\n", wlo, whi,
                std::max(0.0, Nb[1]), std::max(0.0, Nb[2]), "-", "-", "-", "-");
         continue;
      }
      TH1D *shw[nL];
      for (int l = 0; l < nL; ++l) {
         shw[l] = (TH1D *)tplW[l][w]->Clone(Form("shw_%d_%d", l, w));
         shw[l]->SetDirectory(nullptr);
         if (shw[l]->Integral() > 0)
            shw[l]->Scale(1.0 / shw[l]->Integral());
      }
      double e1, e2;
      ylToy(shw, Nb, nToy, e1, e2);
      printf("    %5.0f-%-8.0f %8.0f %8.0f %10.1f %10.1f %10.1f %10.1f\n", wlo, whi, Nb[1], Nb[2],
             100 / std::sqrt(Nb[1]), 100 / std::sqrt(Nb[2]), 100 * e1, 100 * e2);
      if (e1 > 0) {
         best217 = std::min(best217, 100 * e1);
         best332 = std::min(best332, 100 * e2);
         gx.push_back(0.5 * (wlo + whi));
         gy217.push_back(100 * e1);
         gy332.push_back(100 * e2);
      }
      for (int l = 0; l < nL; ++l)
         delete shw[l];
   }

   // ---- what it would take ----------------------------------------------------------------------
   printf("\n  \033[1mWHAT IT WOULD TAKE\033[0m\n");
   if (best332 < 1e8) {
      const double need = (best332 / 10.0) * (best332 / 10.0);
      printf("    best angular point on the 5/2+ is %.1f %%. Reaching 10 %% THERE needs x%.1f the\n", best332, need);
      printf("    integrated luminosity: %.1f days at duty %.2f, or one day at duty %.2f.\n", days * need, duty,
             std::min(1.0, duty * need));
      if (duty * need > 1.0)
         printf("    \033[1;33mThat exceeds duty = 1: the requested beam time cannot deliver 10 %% per angular\n"
                "    point on the 5/2+ in this window, whatever the accelerator does. Either the\n"
                "    deliverable is the window-integrated yield (%.1f %%) rather than a per-point curve,\n"
                "    or the request needs more days, or the window needs widening -- which is what 4 T\n"
                "    buys.\033[0m\n",
                100 * d332);
   }
   printf("    Each ledger factor is worth, in days-equivalent: duty %.2f -> x%.2f, purity %.2f -> x%.2f,\n", duty,
          1 / duty, purity, 1 / purity);
   printf("    pile-up %.3f -> x%.2f, window %.3f -> x%.2f. The window is the biggest and the only one\n", cleanFrac,
          1 / cleanFrac, fWin[0], 1 / fWin[0]);
   printf("    a detector setting (the field) can move.\n");

   // ---- figure -----------------------------------------------------------------------------------
   if (outDir.Length() && gx.size() > 1) {
      gSystem->mkdir(outDir, kTRUE);
      TCanvas c("cy", "", 900, 650);
      c.SetGridy();
      TGraph g1(gx.size(), &gx[0], &gy217[0]), g2(gx.size(), &gx[0], &gy332[0]);
      g1.SetTitle(Form("%s, %.1f d, duty %.2f, R = %.0f;#theta_{lab} [deg];#delta N / N per angular point [%%]",
                       cfgTag.Data(), days, duty, R));
      g1.SetMarkerStyle(20);
      g1.SetMarkerColor(kBlue + 1);
      g1.SetLineColor(kBlue + 1);
      g1.SetLineWidth(2);
      g2.SetMarkerStyle(21);
      g2.SetMarkerColor(kRed + 1);
      g2.SetLineColor(kRed + 1);
      g2.SetLineWidth(2);
      g1.SetMinimum(0);
      g1.Draw("APL");
      g2.Draw("PL");
      TLegend leg(0.55, 0.72, 0.88, 0.88);
      leg.AddEntry(&g1, "1/2^{+} 217 keV", "pl");
      leg.AddEntry(&g2, "5/2^{+} 332 keV", "pl");
      leg.Draw();
      TF1 ten("ten", "10", thLo, thHi);
      ten.SetLineStyle(2);
      ten.SetLineColor(kBlack);
      ten.Draw("same");
      c.SaveAs(outDir + "yields_" + cfgTag + ".png");
      printf("\n    wrote %syields_%s.png\n", outDir.Data(), cfgTag.Data());
   }
   printf("\n");
}
