/// @file decompose_C17.C
/// @brief THE headline result for the 17C M_n/M_p proposal: with the 1/2+ and 5/2+ NOT resolved,
///        can a fixed-position fit still recover their two yields at one day of beam?
///
///   root -b -q 'decompose_C17.C'
///   root -b -q 'decompose_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b285")'
///
/// WHY THIS EXISTS. The proposal does not claim to RESOLVE the doublet. It says (p.4): "while the
/// separation between the two states is only 130 keV [it is 115], the excitation energies between
/// the two states is well known. Therefore, even with the 300 keV energy resolution achievable in
/// the AT-TPC, we expect to be able to fit a doublet in order to extract counts from the
/// individual states." That is a statement about a CONSTRAINED FIT, not about resolution, and no
/// resolution number confirms or refutes it. This macro tests it directly.
///
/// It is also harder than the proposal states, in a way worth being explicit about: it is not a
/// doublet. The ground state sits 217 keV below the 1/2+ and the ELASTIC channel is far stronger
/// than either inelastic, so the two peaks to be extracted sit on the flank of a much larger one
/// less than one resolution width away. The 14C(p,p') precedent the proposal cites (Ref.[24])
/// does not cover this: there the states of interest are at 6-7 MeV, some twenty resolution widths
/// from elastic.
///
/// METHOD.
///   1. Templates: the SIMULATED Ex distribution of each of the three levels, from the campaign.
///      Their positions and shapes are fixed -- exactly the constraint the proposal proposes to
///      use -- and only three amplitudes are free.
///   2. Normalisation: the FRESCO cross sections in fresco/ (8.582 mb for the 1/2+, 3.601 mb for
///      the 5/2+), 300 torr H2 over 1 m, 940 pps, one day. That is 1379 and 579 counts at 4pi
///      before acceptance -- the simulation's own acceptance is already folded in because the
///      templates are built from ACCEPTED events.
///   3. The elastic normalisation is SCANNED, not assumed: R = N_elastic / N_217 over a decade and
///      a half. No FRESCO elastic calculation was supplied with the proposal, and R is the single
///      number that decides whether this works, so it is reported as a dependence rather than
///      folded into one answer. Getting the elastic curve is the highest-value missing input.
///   4. Toys: Poisson-fluctuate the composite, refit the three amplitudes by linear least squares
///      with the templates fixed, repeat. The spread of the recovered yields is the answer.
///
/// The fit is a linear least squares with per-bin variance max(n,1) rather than a Poisson
/// likelihood. For a feasibility estimate at these counts the difference is small, and the linear
/// form makes the amplitude covariance -- which is the thing that actually degrades here -- exact
/// and inspectable rather than a numerical minimiser's output.
///
/// thLabMin CUTS ON RECONSTRUCTED LAB ANGLE, and on (d,d') it is the whole game. sigma(KE) for a
/// deuteron runs from 3.0 MeV at theta_lab 20-30 down to 0.04 MeV at 70-80 -- the recoil is fast and
/// stiff at forward angles and slow and tightly curled at backward ones. Integrated over everything
/// that gives sigma(Ex) = 0.47 MeV and a useless spectrum; above 60 deg it is 0.05-0.14 MeV, as good
/// as the proton channel. The acceptance folding handles the cut automatically: events removed by it
/// simply stop counting as accepted in their theta_cm bin, so the detected yield falls with the
/// resolution gain and the trade is priced rather than assumed.
///
/// ANGULAR WEIGHTING. The generator is flat in cos(theta_cm), so each event is reweighted by its
/// level's FRESCO dsigma/dOmega at its own theta_cm. The ELASTIC is left flat: no elastic
/// calculation exists yet, and since its normalisation is scanned anyway, the only thing the
/// missing shape changes is the elastic template's width. That is flagged in the output, not
/// hidden.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
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
double dcOmega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

/// Two-body inversion, identical in form to acceptance_C14.C:acc_kine and to
/// inel_summary_C17.C:s17_ex -- measured ejectile (theta_lab, KE) plus a beam energy -> residual
/// excitation energy. Needed here because the vertex correction has to be REDONE from the tracked
/// quantities; the tree's exReco is the CONSTANT-Ebeam reconstruction and nothing else.
double dcEx(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * dcOmega2(s, m1 * m1, m2 * m2) * dcOmega2(uu, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                         (2 * m2 * m2) +
                      s + uu - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

/// Read a two-column FRESCO xsec table (theta_cm [deg], dsigma/dOmega [mb/sr]).
bool dcReadFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
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

/// Linear interpolation on a 1 deg grid; 0 outside.
double dcInterp(const std::vector<double> &th, const std::vector<double> &xs, double x)
{
   if (th.empty() || x < th.front() || x > th.back())
      return 0;
   const size_t i = std::min((size_t)(x - th.front()), th.size() - 2);
   const double f = (x - th[i]) / (th[i + 1] - th[i]);
   return xs[i] * (1 - f) + xs[i + 1] * f;
}

/// Solve the 3x3 symmetric normal equations A a = b by Gaussian elimination with partial pivoting.
/// Returns false if singular. On success `cov` holds A^-1, whose diagonal is the amplitude
/// variance -- the quantity that blows up when two templates overlap.
bool dcSolve3(double A[3][3], const double b[3], double a[3], double cov[3][3])
{
   double M[3][6];
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         M[i][j] = A[i][j];
         M[i][3 + j] = (i == j) ? 1.0 : 0.0;
      }
   }
   for (int col = 0; col < 3; ++col) {
      int piv = col;
      for (int r = col + 1; r < 3; ++r)
         if (std::fabs(M[r][col]) > std::fabs(M[piv][col]))
            piv = r;
      if (std::fabs(M[piv][col]) < 1e-30)
         return false;
      if (piv != col)
         for (int j = 0; j < 6; ++j)
            std::swap(M[col][j], M[piv][j]);
      const double d = M[col][col];
      for (int j = 0; j < 6; ++j)
         M[col][j] /= d;
      for (int r = 0; r < 3; ++r) {
         if (r == col)
            continue;
         const double f = M[r][col];
         for (int j = 0; j < 6; ++j)
            M[r][j] -= f * M[col][j];
      }
   }
   for (int i = 0; i < 3; ++i) {
      a[i] = 0;
      for (int j = 0; j < 3; ++j) {
         cov[i][j] = M[i][3 + j];
         a[i] += M[i][3 + j] * b[j];
      }
   }
   return true;
}
} // namespace

void decompose_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString cfgTag = "pp_b285",
                   TString frescoDir = "./fresco/", Int_t nToy = 400, Double_t chi2Cut = 5.0,
                   Bool_t vertexCorrected = kTRUE, TString outDir = "./plots/", Double_t thLabMin = 0.0)
{
   gStyle->SetOptStat(0);
   gRandom->SetSeed(20260902);

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};
   const char *stName[nL] = {"3/2+ g.s.", "1/2+ 217", "5/2+ 332"};

   TString chan(cfgTag);
   chan.Remove(2);
   TString btag(cfgTag);
   btag.Remove(0, 3);

   const double uAmu = 931.49401;
   const double mBeam = 17.0225787 * uAmu;                       // 17C
   const double mRes = mBeam;                                    // scattering: residual IS the beam
   const double mLight = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu;
   const double EbeamConst = 135.0;
   double ebz_a = 0, ebz_b = 0;
   bool haveEbz = false;

   printf("\n\033[1;33m##########################################################################\033[0m\n");
   printf("\033[1;33m 17C triplet decomposition -- %s, %s Ex, theta_lab > %.0f\033[0m\n", cfgTag.Data(),
          vertexCorrected ? "vertex-corrected" : "constant-Ebeam", thLabMin);
   printf("\033[1;33m##########################################################################\033[0m\n");

   // ---- FRESCO cross sections and the one-day yields -------------------------------------------
   std::vector<double> fth[2], fxs[2];
   const bool haveFresco = dcReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) &&
                           dcReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1]);
   if (!haveFresco) {
      printf("\033[1;31m  no FRESCO tables in %s -- cannot normalise. Aborting.\033[0m\n", frescoDir.Data());
      return;
   }
   double sig4pi[2] = {0, 0};
   for (int k = 0; k < 2; ++k)
      for (size_t i = 0; i < fth[k].size(); ++i)
         sig4pi[k] += fxs[k][i] * 2 * TMath::Pi() * std::sin(fth[k][i] * TMath::DegToRad()) * (TMath::Pi() / 180.0);

   // 300 torr H2 (media.geo H_300torr_RT: atomic H, A = 1.007, rho = 3.308e-5) over 1 m; D2 at the
   // same pressure has the same molecular number density to 0.8 %, so the deuteron day sees the
   // same areal density of target nuclei.
   const double nTgt = 3.308e-5 / 1.007 * 6.02214076e23 * 100.0; // nuclei/cm2
   const double beamPerDay = 940.0 * 86400.0;
   const double N217 = beamPerDay * nTgt * sig4pi[0] * 1e-27;
   const double N332 = beamPerDay * nTgt * sig4pi[1] * 1e-27;
   printf("  FRESCO: sigma(217) = %.3f mb, sigma(332) = %.3f mb\n", sig4pi[0], sig4pi[1]);
   printf("  one day at 940 pps, 300 torr, 1 m:  N(217) = %.0f   N(332) = %.0f   (4pi, before acceptance)\n", N217,
          N332);
   if (chan == "dd")
      printf("  \033[1;33mNOTE: these are the PROTON cross sections. No (d,d') calculation was supplied;\n"
             "  the deuteron yields are assumed equal, which the proposal also assumes.\033[0m\n");

   // ---- templates from the campaign ------------------------------------------------------------
   const double exLo = -1.5, exHi = 2.0;
   const int nB = 140; // 25 keV bins
   TH1D *tpl[nL];
   long nAcc[nL] = {0, 0, 0};
   // accepted theta_cm (unweighted) and the number of reactions GENERATED, per level. Together
   // with the fact that the generator is flat in cos(theta_cm) over [CM_LO, CM_HI], these give the
   // acceptance as a function of angle -- which is what turns a 4pi cross section into counts that
   // actually land in the spectrum.
   std::vector<double> cmAcc[nL];
   double nGen[nL] = {0, 0, 0};
   for (int l = 0; l < nL; ++l) {
      tpl[l] = new TH1D(Form("tpl_%s", stTag[l]), "", nB, exLo, exHi);
      tpl[l]->SetDirectory(nullptr);
      // Quote ONLY the directory -- quoting the whole pattern stops the "s*" wildcard expanding.
      TString dirq = TString("\"") + root + "/" + cfgTag + "\"";
      TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
      TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-10s MISSING (%s)\033[0m\n", stTag[l], pat.Data());
         return;
      }
      // The generated-reaction count lives in the acceptance log, which the same sample wrote.
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

      // MEASURE E_beam(z) from truth on the first level present, by solving per event for the beam
      // energy that reproduces that level's GENERATED Ex. Same procedure as inel_summary_C17.C.
      // Without this the templates are the constant-Ebeam reconstruction, which is 2.3x wider --
      // and mislabelling that as vertex-corrected would overstate the difficulty by the same
      // factor.
      if (vertexCorrected && !haveEbz) {
         TGraph g;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = EbeamConst - 40, hi = EbeamConst + 40;
            if ((dcEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                   (dcEx(mBeam, mLight, mLight, mRes, hi, thT, keT) - ExGen[l]) >
                0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((dcEx(mBeam, mLight, mLight, mRes, lo, thT, keT) - ExGen[l]) *
                      (dcEx(mBeam, mLight, mLight, mRes, mid, thT, keT) - ExGen[l]) <=
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
            printf("  E_beam(z) measured from truth: %.4f %+.6f * z[mm]  (n = %d)\n", ebz_a, ebz_b, g.GetN());
         }
      }

      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (c2n >= chi2Cut || thR < thLabMin)
            continue;
         // The quantity that goes into the template. exReco is the CONSTANT-Ebeam value; the
         // vertex-corrected one has to be rebuilt from the tracked (theta, KE) at E_beam(z_reco).
         if (vertexCorrected && haveEbz) {
            const double ev = dcEx(mBeam, mLight, mLight, mRes, ebz_a + ebz_b * zR, thR, keR);
            if (ev < -1e8)
               continue;
            exR = ev;
         }
         // Reweight the flat generator to the real angular distribution. The elastic has no
         // calculation, so it stays flat -- see the header note.
         double w = 1.0;
         if (l == 1)
            w = dcInterp(fth[0], fxs[0], cmT);
         else if (l == 2)
            w = dcInterp(fth[1], fxs[1], cmT);
         if (w <= 0)
            continue;
         tpl[l]->Fill(exR, w);
         cmAcc[l].push_back(cmT);
         ++nAcc[l];
      }
      f->Close();
      if (tpl[l]->Integral() <= 0) {
         printf("\033[1;31m  %-10s empty template\033[0m\n", stTag[l]);
         return;
      }
      tpl[l]->Scale(1.0 / tpl[l]->Integral()); // unit-normalised shape
   }
   printf("  template statistics (accepted, chi2/ndf < %.1f):  gs %ld   217 %ld   332 %ld\n", chi2Cut, nAcc[0],
          nAcc[1], nAcc[2]);

   // ---- fold the acceptance in, per angle ------------------------------------------------------
   // The generator is FLAT IN cos(theta_cm) over [cmLo, cmHi], so the number generated in an angular
   // bin is nGen x Dcos(bin)/Dcos(range). The acceptance in that bin is accepted/generated, and the
   // detected fraction of a real 4pi cross section is the FRESCO weight of each bin times its
   // acceptance. Skipping this step treats every one of the 1379 counts as detected and understates
   // the errors by ~10 %.
   const double cmLo = 10.0, cmHi = 178.0;
   const int nAB = 24; // 7 deg bins
   auto binOf = [&](double th) { return (int)((th - cmLo) / (cmHi - cmLo) * nAB); };
   double detFrac[nL] = {0, 0, 0};
   for (int l = 0; l < nL; ++l) {
      if (nGen[l] <= 0) {
         printf("\033[1;31m  no generated-reaction count for %s; acceptance not folded\033[0m\n", stTag[l]);
         detFrac[l] = 1.0;
         continue;
      }
      std::vector<double> accN(nAB, 0.0);
      for (double th : cmAcc[l]) {
         const int b = binOf(th);
         if (b >= 0 && b < nAB)
            accN[b] += 1.0;
      }
      const double dcosTot = std::cos(cmLo * TMath::DegToRad()) - std::cos(cmHi * TMath::DegToRad());
      const std::vector<double> &fth_l = (l == 2) ? fth[1] : fth[0];
      const std::vector<double> &fxs_l = (l == 2) ? fxs[1] : fxs[0];
      const double sig_l = (l == 2) ? sig4pi[1] : sig4pi[0];
      for (int b = 0; b < nAB; ++b) {
         const double lo = cmLo + (cmHi - cmLo) * b / nAB;
         const double hi = cmLo + (cmHi - cmLo) * (b + 1) / nAB;
         const double dcos = std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad());
         const double gen = nGen[l] * dcos / dcosTot;
         if (gen <= 0)
            continue;
         const double A = accN[b] / gen;
         // FRESCO weight of this bin, as a fraction of 4pi
         double sb = 0;
         for (size_t i = 0; i < fth_l.size(); ++i)
            if (fth_l[i] >= lo && fth_l[i] < hi)
               sb += fxs_l[i] * 2 * TMath::Pi() * std::sin(fth_l[i] * TMath::DegToRad()) * (TMath::Pi() / 180.0);
         detFrac[l] += (sb / sig_l) * A;
      }
   }
   printf("  detected fraction of the 4pi yield (acceptance folded with the FRESCO shape): "
          "gs %.3f  217 %.3f  332 %.3f\n",
          detFrac[0], detFrac[1], detFrac[2]);

   // ---- the elastic scan -------------------------------------------------------------------------
   const int nR = 5;
   const double Rscan[nR] = {0.0, 3.0, 10.0, 30.0, 100.0};

   printf("\n  \033[1;33mFIXED-POSITION THREE-COMPONENT FIT, %d toys per point\033[0m\n", nToy);
   printf("  positions and shapes frozen at the simulated lineshapes; only the three yields float\n");
   printf("  R = N(elastic) / N(1/2+ 217). No elastic calculation was supplied, so it is scanned.\n\n");
   printf("  %8s %10s %12s %12s %12s %12s %10s\n", "R", "N_elastic", "d(217)/217", "d(332)/332", "bias(217)",
          "bias(332)", "corr");
   printf("  %8s %10s %12s %12s %12s %12s %10s\n", "", "[counts]", "[%]", "[%]", "[%]", "[%]", "217:332");

   std::vector<double> rr;
   std::vector<double> e217, e332;
   for (int ir = 0; ir < nR; ++ir) {
      const double N217d = N217 * detFrac[1], N332d = N332 * detFrac[2];
      const double Ngs = Rscan[ir] * N217 * detFrac[0];
      const double Ntrue[nL] = {Ngs, N217d, N332d};

      // expected composite
      std::vector<double> mu(nB, 0.0);
      for (int b = 1; b <= nB; ++b)
         for (int l = 0; l < nL; ++l)
            mu[b - 1] += Ntrue[l] * tpl[l]->GetBinContent(b);

      double s217 = 0, s2217 = 0, s332 = 0, s2332 = 0, s217332 = 0;
      int nOk = 0;
      for (int it = 0; it < nToy; ++it) {
         // Poisson toy
         std::vector<double> n(nB);
         for (int b = 0; b < nB; ++b)
            n[b] = gRandom->Poisson(mu[b]);

         // linear least squares with per-bin variance max(n,1)
         double A[3][3] = {{0}}, rhs[3] = {0, 0, 0};
         for (int b = 0; b < nB; ++b) {
            const double v = std::max(n[b], 1.0);
            double T[3];
            for (int l = 0; l < nL; ++l)
               T[l] = tpl[l]->GetBinContent(b + 1);
            for (int i = 0; i < 3; ++i) {
               rhs[i] += T[i] * n[b] / v;
               for (int j = 0; j < 3; ++j)
                  A[i][j] += T[i] * T[j] / v;
            }
         }
         double a[3], cov[3][3];
         if (!dcSolve3(A, rhs, a, cov))
            continue;
         ++nOk;
         s217 += a[1];
         s2217 += a[1] * a[1];
         s332 += a[2];
         s2332 += a[2] * a[2];
         s217332 += a[1] * a[2];
      }
      if (!nOk) {
         printf("  %8.0f %10.0f   singular\n", Rscan[ir], Ngs);
         continue;
      }
      const double m217 = s217 / nOk, m332 = s332 / nOk;
      const double v217 = std::max(0.0, s2217 / nOk - m217 * m217);
      const double v332 = std::max(0.0, s2332 / nOk - m332 * m332);
      const double c = (s217332 / nOk - m217 * m332) / std::sqrt(std::max(1e-30, v217 * v332));
      printf("  %8.0f %10.0f %12.1f %12.1f %12.1f %12.1f %10.2f\n", Rscan[ir], Ngs, 100 * std::sqrt(v217) / N217d,
             100 * std::sqrt(v332) / N332d, 100 * (m217 - N217d) / N217d, 100 * (m332 - N332d) / N332d, c);
      rr.push_back(Rscan[ir]);
      e217.push_back(100 * std::sqrt(v217) / N217d);
      e332.push_back(100 * std::sqrt(v332) / N332d);
   }

   printf("\n  Detected counts in one day: N(217) = %.0f, N(332) = %.0f. The statistical floor with\n"
          "  PERFECT separation would be 1/sqrt(N): %.1f %% (217) and %.1f %% (332).\n",
          N217 * detFrac[1], N332 * detFrac[2], 100 / std::sqrt(N217 * detFrac[1]),
          100 / std::sqrt(N332 * detFrac[2]));
   printf("  Anything above that is the price of the overlap. The proposal needs ~10 %% on the\n");
   printf("  deformation, which is ~10 %% on each yield, per angular bin.\n");

   // ---- the figure -------------------------------------------------------------------------------
   if (outDir.Length()) {
      gSystem->mkdir(outDir, kTRUE);
      TCanvas *cv = new TCanvas("cvdec", "decomposition", 1100, 750);
      const double Rshow = 10.0;
      TH1D *sum = (TH1D *)tpl[0]->Clone("sumshow");
      sum->Reset();
      const int col[nL] = {kGray + 2, kRed + 1, kBlue + 1};
      const double Nshow[nL] = {Rshow * N217 * detFrac[0], N217 * detFrac[1], N332 * detFrac[2]};
      for (int l = 0; l < nL; ++l) {
         TH1D *h = (TH1D *)tpl[l]->Clone(Form("show_%d", l));
         h->Scale(Nshow[l]);
         h->SetLineColor(col[l]);
         h->SetLineWidth(2);
         sum->Add(h);
         if (l == 0) {
            sum->SetTitle(Form("%s, one day, N_{elastic}/N_{217} = %.0f;"
                               "E_{x}(^{17}C) [MeV];counts / 25 keV",
                               cfgTag.Data(), Rshow));
         }
      }
      sum->SetLineColor(kBlack);
      sum->SetLineWidth(3);
      sum->Draw("HIST");
      TLegend *lg = new TLegend(0.60, 0.62, 0.89, 0.88);
      lg->SetBorderSize(0);
      lg->AddEntry(sum, "sum (what is measured)", "l");
      for (int l = 0; l < nL; ++l) {
         TH1D *h = (TH1D *)gDirectory->Get(Form("show_%d", l));
         if (h) {
            h->Draw("HIST SAME");
            lg->AddEntry(h, stName[l], "l");
         }
      }
      lg->Draw();
      cv->SaveAs(outDir + "decompose_" + cfgTag + ".png");
      printf("\n  wrote %sdecompose_%s.png\n", outDir.Data(), cfgTag.Data());
   }
   printf("\n  decompose done\n\n");
}
