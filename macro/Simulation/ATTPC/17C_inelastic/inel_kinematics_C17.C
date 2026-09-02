/// @file inel_kinematics_C17.C
/// @brief Analytic 17C(p,p') and 17C(d,d') kinematics for the M_n/M_p proposal
///        (C17p_FRIB_Proposal.pdf, this directory). No simulation input: this is the map the
///        simulation is measured against, and the GATE that decides whether the campaign is
///        worth running at all.
///
///   root -b -q 'inel_kinematics_C17.C'            // at the mean vertex beam energy, 135 MeV
///   root -b -q 'inel_kinematics_C17.C(142.29)'    // at the chamber entrance
///
/// WHAT THIS HAS TO ANSWER. The proposal measures d(sigma)/d(Omega) for the 1/2+ (217 keV) and
/// 5/2+ (332 keV) states of 17C. Those two are 115 keV apart -- NOT the 130 keV the proposal
/// text states -- and the ground state is 217 keV below the 1/2+. All three therefore sit inside
/// ONE 300 keV resolution width, with the elastic channel one to two orders of magnitude stronger
/// than either inelastic. The question is not "can 115 keV be resolved" but "how much proton
/// energy separates these states, and is that more or less than the tracking resolves".
///
/// WHAT THIS MACRO CAN AND CANNOT DO -- read before trusting the sigma(Ex) block below.
/// It computes the LEVERAGE exactly: dEx/dKE, dEx/dtheta and dEx/dEbeam are properties of the
/// kinematics and nothing else, and dEx/dKE turns out to be flat at -0.533 (p,p') / -0.563 (d,d')
/// across every angle, so the 115 keV level gap always appears as ~217 keV of recoil energy.
/// That part is solid and is the useful output.
///
/// It CANNOT predict sigma(Ex), and the block below that tries to is kept only as a worked
/// example of the propagation, with its input clearly wrong. It feeds in the single angle-averaged
/// sigma(KE) = 0.343 MeV that 14C(p,p') quotes for its elastic channel (14C_pp/highfield/
/// RESULTS.md) and gets a sigma(Ex) flat in angle at 0.18 MeV. The campaign then measured
/// sigma(KE) on this channel and it is NOT one number: it runs from 0.49 MeV at theta_lab 20-30
/// to 0.023 MeV at 60-70, a factor of twenty, because a slower recoil makes a tighter helix that
/// is measured far better. The real sigma(Ex) follows that, not the flat line.
///
/// The second thing the propagation misses is that the constant-Ebeam term DOMINATES until it is
/// corrected: measured sigma(Ex) = 0.225 MeV against a method floor of 0.206 MeV, so the detector
/// is nearly invisible and the field changes almost nothing. See RESULTS.md.
///
/// The requirement is stated instead: what sigma(KE) a given separation needs. Compare it with the
/// MEASURED per-slice sigma(KE) that inel_summary_C17.C prints.
///
/// ANGLE CONVENTION -- GET THIS WRONG AND EVERY TABLE IS SILENTLY MIRRORED. theta_cm is the
/// projectile (17C) centre-of-mass scattering angle. For a light recoil this is the supplement of
/// the recoil's cm polar angle about the beam, which is exactly what acceptance_C14.C computes
/// (theta_cm = pi - acos(...), acceptance_C14.C:44). The elastic check is that the recoil comes
/// out at theta_lab = (180 - theta_cm)/2; this macro asserts that numerically rather than
/// trusting it, because the (d,p) arm of this proposal shipped a reversed table for exactly this
/// reason (17C_dp/RESULTS.md, "theta_cm CONVENTION").
///
/// FRESCO's tabulated angle is the same one: its output is labelled "for projectile", partition 1.
/// So the distributions in fresco/ can weight these tables directly, and the yield-weighted rows
/// below say which lab angles actually carry the measurement.
///
/// @param Ebeam    17C kinetic energy in MeV. 142.29 at the chamber entrance (8.37 MeV/u, the
///                 proposal value); ~14.5 MeV is lost crossing the metre, so 135 is the mean over
///                 a uniform vertex distribution and 136 is what the FRESCO runs used.
/// @param outDir   where to write kinematics_C17inel.png. Empty disables plotting.
/// @param frescoDir directory holding c17pp_217keV.out / c17pp_332keV.out. Empty disables the
///                 yield weighting (the tables are still printed, unweighted).

#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TSystem.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace
{
const double kU = 931.49401;      // MeV per amu
const double kM17C = 17.02257865; // amu, AME2020
const double kMp = 1.0078250322;
const double kMd = 2.0141017778;

/// The three levels the proposal is about. 332 keV is the proposal's value for the 5/2+; ENSDF
/// adopts 331 keV. The 1 keV is irrelevant to everything below, but the separation it implies
/// (115 keV) is NOT the 130 keV the proposal quotes.
const int kNL = 3;
const double kEx[kNL] = {0.0, 0.217, 0.332};
const char *kName[kNL] = {"3/2+ g.s.", "1/2+ 217", "5/2+ 332"};

/// Read a two-column FRESCO xsec table (theta_cm [deg], dsigma/dOmega [mb/sr]). Header and
/// grace-format lines all start with '#' or '@'.
bool readFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
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

/// Integrate a tabulated dsigma/dOmega over 4pi, in mb.
double integrate4pi(const std::vector<double> &th, const std::vector<double> &xs)
{
   double s = 0;
   const double dth = TMath::Pi() / 180.0; // the tables are on a 1 deg grid
   for (size_t i = 0; i < th.size(); ++i)
      s += xs[i] * 2 * TMath::Pi() * std::sin(th[i] * TMath::DegToRad()) * dth;
   return s;
}
} // namespace

void inel_kinematics_C17(Double_t Ebeam = 135.0, TString outDir = "./plots/", TString frescoDir = "./fresco/")
{
   gStyle->SetOptStat(0);

   // ---------------------------------------------------------------------------------------
   // The two channels. In both, the BEAM is excited and the light target particle is what the
   // AT-TPC fits: ejectile mass = target mass, residual mass = beam mass.
   // ---------------------------------------------------------------------------------------
   struct Chan {
      const char *tag;
      const char *label;
      double mLight; // target = ejectile, amu
      int color;
   };
   const Chan chan[2] = {{"pp", "17C(p,p')", kMp, kRed + 1}, {"dd", "17C(d,d')", kMd, kBlue + 1}};

   printf("\n\033[1;33m===================================================================\033[0m\n");
   printf("\033[1;33m 17C inelastic kinematics for the M_n/M_p proposal\033[0m\n");
   printf("\033[1;33m Ebeam = %.2f MeV = %.3f MeV/u   (entrance 142.29 = 8.37 MeV/u,\033[0m\n", Ebeam, Ebeam / 17.0);
   printf("\033[1;33m ~14.5 MeV lost over the metre, FRESCO ran at 136)\033[0m\n");
   printf("\033[1;33m===================================================================\033[0m\n");
   printf("  levels: %s (0), %s (%.3f), %s (%.3f)   ->  1/2+ to 5/2+ separation = %.0f keV\n", kName[0], kName[1],
          kEx[1], kName[2], kEx[2], (kEx[2] - kEx[1]) * 1000);
   printf("  NOTE: the proposal text says 130 keV for that separation. 332 - 217 = 115.\n");

   // 14C(p,p') measured tracking resolution at this field and pad plane -- the transfer basis for
   // the sigma(Ex) prediction below (14C_pp/highfield/RESULTS.md, b285_attpc, elastic).
   const double sigKE_ref = 0.343;  // MeV
   const double sigTh_ref = 0.099;  // deg

   // FRESCO angular distributions, if present.
   std::vector<double> fth[2], fxs[2];
   bool haveFresco = false;
   if (frescoDir.Length()) {
      const bool ok1 = readFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]);
      const bool ok2 = readFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1]);
      haveFresco = ok1 && ok2;
      if (haveFresco) {
         const double s1 = integrate4pi(fth[0], fxs[0]), s2 = integrate4pi(fth[1], fxs[1]);
         printf("\n\033[1;36m--- FRESCO (p,p') angular distributions, %s ---\033[0m\n", frescoDir.Data());
         printf("  sigma(217) = %.3f mb    sigma(332) = %.3f mb    ratio = %.2f\n", s1, s2, s1 / s2);
         // Counts per day: H_300torr_RT is ATOMIC hydrogen in media.geo (A = 1.007, Z = 1,
         // rho = 3.308e-5 g/cm3), so n = rho/A * N_A directly, over the 1 m active length.
         const double nH = 3.308e-5 / 1.007 * 6.02214076e23 * 100.0; // atoms/cm2
         const double perDay = 940.0 * 86400.0;
         printf("  at 300 torr H2 x 1 m x 940 pps x 1 day, 4pi:  %.0f counts (217)   %.0f counts (332)\n",
                perDay * nH * s1 * 1e-27, perDay * nH * s2 * 1e-27);
         printf("  the proposal quotes ~1000 scattered protons: right for the 1/2+, ~2x optimistic\n");
         printf("  for the 5/2+, which is therefore the state that sets the statistical error.\n");
      } else {
         printf("\n  [no FRESCO tables in %s -- tables below are unweighted]\n", frescoDir.Data());
      }
   }

   // Graphs for the figure, [channel][level].
   TGraph *gThLab[2][kNL], *gKE[2][kNL], *gSigEx[2][kNL];
   TGraph *gSepKE[2];               // KE separation of the 1/2+ and 5/2+ vs theta_cm
   TGraph *gDerKE[2], *gDerTh[2];   // |dEx/dKE| and |dEx/dtheta| vs theta_lab, for the field panel

   for (int c = 0; c < 2; ++c) {
      const double mB = kM17C * kU;   // beam 17C
      const double mT = chan[c].mLight * kU; // target p or d
      const double mE = mT;           // ejectile: the same light particle, this is scattering
      const double mR0 = mB;          // residual: 17C, excited by Ex

      const double Eb = Ebeam + mB;
      const double pb = std::sqrt(Eb * Eb - mB * mB);
      TLorentzVector Lb(0, 0, pb, Eb), Lt(0, 0, 0, mT);
      TLorentzVector W = Lb + Lt;
      TVector3 bst = W.BoostVector();
      const double s = W.M2();

      /// Lab (theta, KE) of the light ejectile at projectile cm angle thcmDeg with the residual
      /// left at excitation ex.
      auto labOf = [&](double thcmDeg, double ex, double &thLab, double &ke) {
         const double mR = mR0 + ex;
         const double arg = (s - (mR + mE) * (mR + mE)) * (s - (mR - mE) * (mR - mE));
         if (arg <= 0) {
            thLab = ke = -1;
            return false;
         }
         const double pcm = std::sqrt(arg) / (2 * std::sqrt(s));
         // supplement: theta_cm is the PROJECTILE angle, the ejectile is the recoil
         const double th = (180.0 - thcmDeg) * TMath::DegToRad();
         TLorentzVector L(pcm * std::sin(th), 0, pcm * std::cos(th), std::sqrt(pcm * pcm + mE * mE));
         L.Boost(bst);
         thLab = L.Vect().Theta() * TMath::RadToDeg();
         ke = L.E() - mE;
         return true;
      };

      /// The inversion the analysis performs: excitation energy from the MEASURED (KE, theta_lab)
      /// of the light track, at an ASSUMED beam energy.
      auto exOf = [&](double keM, double thM, double EbeamAssumed) {
         const double E = keM + mE;
         const double p = std::sqrt(std::max(0., E * E - mE * mE));
         const double t = thM * TMath::DegToRad();
         TLorentzVector Le(p * std::sin(t), 0, p * std::cos(t), E);
         const double EbA = EbeamAssumed + mB;
         TLorentzVector LbA(0, 0, std::sqrt(EbA * EbA - mB * mB), EbA);
         return ((LbA + Lt) - Le).M() - mR0;
      };

      // --- convention assertion: elastic recoil must come out at (180 - theta_cm)/2 ----------
      {
         double th, ke;
         labOf(60.0, 0.0, th, ke);
         const double expect = (180.0 - 60.0) / 2.0;
         printf("\n\033[1;36m--- %s, Ebeam = %.2f MeV ---\033[0m\n", chan[c].label, Ebeam);
         printf("  convention check: elastic theta_cm = 60 deg  ->  theta_lab = %.2f deg "
                "(expect %.2f, (180-th_cm)/2)  %s\n",
                th, expect, std::fabs(th - expect) < 0.5 ? "\033[1;32mOK\033[0m" : "\033[1;31mMISMATCH\033[0m");
      }

      // dKE_gap is the quantity a detector actually has to resolve: how far apart the 1/2+ and the
      // 5/2+ recoils are IN ENERGY, at the SAME lab angle. Comparing KE at the same theta_cm is
      // the wrong comparison -- theta_lab moves with the level too, and the detector does not
      // measure theta_cm.
      printf("  %8s %10s %10s %12s %12s %12s %11s\n", "theta_cm", "theta_lab", "KE(gs)", "dEx/dKE", "dEx/dtheta",
             "dEx/dEbeam", "dKE(gap)");
      printf("  %8s %10s %10s %12s %12s %12s %11s\n", "[deg]", "[deg]", "[MeV]", "", "[MeV/deg]", "", "[keV]");
      for (double tcm : {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140., 160.}) {
         double th0, ke0;
         if (!labOf(tcm, kEx[0], th0, ke0))
            continue;
         const double dKE = 0.005, dTh = 0.01, dEb = 0.05;
         const double dExdKE = (exOf(ke0 + dKE, th0, Ebeam) - exOf(ke0 - dKE, th0, Ebeam)) / (2 * dKE);
         const double dExdTh = (exOf(ke0, th0 + dTh, Ebeam) - exOf(ke0, th0 - dTh, Ebeam)) / (2 * dTh);
         const double dExdEb = (exOf(ke0, th0, Ebeam + dEb) - exOf(ke0, th0, Ebeam - dEb)) / (2 * dEb);
         const double gapKE = (kEx[2] - kEx[1]) / std::fabs(dExdKE) * 1000.0; // keV
         printf("  %8.0f %10.1f %10.3f %12.4f %12.5f %12.4f %11.0f\n", tcm, th0, ke0, dExdKE, dExdTh, dExdEb, gapKE);
      }

      // --- the gate: sigma(Ex) predicted from the 14C-measured tracking resolution -----------
      printf("\n  \033[1;33mWORKED EXAMPLE ONLY -- propagating ONE angle-averaged sigma(KE) from 14C(p,p').\033[0m\n");
      printf("  \033[1;31m  Its flatness is an ARTEFACT of using a single sigma(KE): the campaign measured\n");
      printf("    sigma(KE) varying 20x across these slices. Use inel_summary_C17.C's measured table.\033[0m\n");
      printf("  (sigma(KE) = %.3f MeV, sigma(theta) = %.3f deg; the 1/2+ to 5/2+ gap is 0.115 MeV)\n", sigKE_ref,
             sigTh_ref);
      printf("  %8s %10s %12s %12s %12s %10s\n", "theta_cm", "theta_lab", "from KE", "from theta", "quadrature",
             "gap/sigma");
      for (double tcm : {20., 30., 40., 50., 60., 70., 80., 90., 100., 120., 140.}) {
         double th0, ke0;
         if (!labOf(tcm, kEx[0], th0, ke0))
            continue;
         const double dKE = 0.005, dTh = 0.01;
         const double dExdKE = (exOf(ke0 + dKE, th0, Ebeam) - exOf(ke0 - dKE, th0, Ebeam)) / (2 * dKE);
         const double dExdTh = (exOf(ke0, th0 + dTh, Ebeam) - exOf(ke0, th0 - dTh, Ebeam)) / (2 * dTh);
         const double sKE = std::fabs(dExdKE) * sigKE_ref;
         const double sTh = std::fabs(dExdTh) * sigTh_ref;
         const double tot = std::sqrt(sKE * sKE + sTh * sTh);
         printf("  %8.0f %10.1f %12.4f %12.4f %12.4f %10.2f\n", tcm, th0, sKE, sTh, tot, 0.115 / tot);
      }

      // --- what sigma(KE) would be needed ------------------------------------------------------
      // sigma(KE) is the whole detector budget here (the theta term is a tenth of it), so the
      // requirement can be stated cleanly. What CANNOT be done from this macro alone is predict
      // sigma(KE) itself -- see the warning below.
      {
         double th0, ke0;
         const double tRef = 78.0; // the FRESCO yield-weighted mean theta_cm
         labOf(tRef, kEx[1], th0, ke0);
         const double dKE = 0.005;
         const double dExdKE = (exOf(ke0 + dKE, th0, Ebeam) - exOf(ke0 - dKE, th0, Ebeam)) / (2 * dKE);
         printf("\n  \033[1;33mTHE REQUIREMENT (at theta_cm = %.0f, theta_lab = %.0f, where the yield is)\033[0m\n",
                tRef, th0);
         printf("  %14s %12s %12s\n", "separation", "sigma(Ex)", "sigma(KE)");
         for (double sep : {0.5, 1.0, 2.0}) {
            const double need = 0.115 / (2 * sep);
            printf("  %14.1f %12.3f %12.3f\n", sep, need, need / std::fabs(dExdKE));
         }
      }

      // --- build the graphs -----------------------------------------------------------------
      for (int l = 0; l < kNL; ++l) {
         gThLab[c][l] = new TGraph();
         gKE[c][l] = new TGraph();
         gSigEx[c][l] = new TGraph();
         int n = 0, m = 0;
         for (double tcm = 2; tcm <= 178; tcm += 1.0) {
            double th, ke;
            if (!labOf(tcm, kEx[l], th, ke) || ke <= 0)
               continue;
            gThLab[c][l]->SetPoint(n, tcm, th);
            gKE[c][l]->SetPoint(n, th, ke);
            ++n;
            const double dKE = 0.005, dTh = 0.01;
            const double dExdKE = (exOf(ke + dKE, th, Ebeam) - exOf(ke - dKE, th, Ebeam)) / (2 * dKE);
            const double dExdTh = (exOf(ke, th + dTh, Ebeam) - exOf(ke, th - dTh, Ebeam)) / (2 * dTh);
            const double tot =
               std::sqrt(std::pow(dExdKE * sigKE_ref, 2) + std::pow(dExdTh * sigTh_ref, 2));
            gSigEx[c][l]->SetPoint(m++, th, tot);
         }
      }

      // leverage vs lab angle, evaluated on the 1/2+ level (the one carrying the yield)
      gDerKE[c] = new TGraph();
      gDerTh[c] = new TGraph();
      {
         int n = 0;
         for (double tcm = 5; tcm <= 175; tcm += 1.0) {
            double th, ke;
            if (!labOf(tcm, kEx[1], th, ke) || ke <= 0.05)
               continue;
            const double dKE = 0.005, dTh = 0.01;
            gDerKE[c]->SetPoint(n, th,
                                std::fabs((exOf(ke + dKE, th, Ebeam) - exOf(ke - dKE, th, Ebeam)) / (2 * dKE)));
            gDerTh[c]->SetPoint(n, th,
                                std::fabs((exOf(ke, th + dTh, Ebeam) - exOf(ke, th - dTh, Ebeam)) / (2 * dTh)));
            ++n;
         }
      }

      // KE separation of the two proposal states, vs theta_cm
      gSepKE[c] = new TGraph();
      int k = 0;
      for (double tcm = 2; tcm <= 178; tcm += 1.0) {
         double t1, k1, t2, k2;
         if (!labOf(tcm, kEx[1], t1, k1) || !labOf(tcm, kEx[2], t2, k2))
            continue;
         gSepKE[c]->SetPoint(k++, tcm, (k1 - k2) * 1000.0); // keV
      }

      // --- yield-weighted summary -----------------------------------------------------------
      if (haveFresco && c == 0) {
         double sw = 0, swTh = 0, swThLab = 0, swSig = 0;
         for (size_t i = 0; i < fth[0].size(); ++i) {
            const double w = fxs[0][i] * std::sin(fth[0][i] * TMath::DegToRad());
            double th, ke;
            if (!labOf(fth[0][i], kEx[1], th, ke))
               continue;
            const double dKE = 0.005, dTh = 0.01;
            const double dExdKE = (exOf(ke + dKE, th, Ebeam) - exOf(ke - dKE, th, Ebeam)) / (2 * dKE);
            const double dExdTh = (exOf(ke, th + dTh, Ebeam) - exOf(ke, th - dTh, Ebeam)) / (2 * dTh);
            sw += w;
            swTh += w * fth[0][i];
            swThLab += w * th;
            swSig += w * std::sqrt(std::pow(dExdKE * sigKE_ref, 2) + std::pow(dExdTh * sigTh_ref, 2));
         }
         printf("\n  \033[1;33mFRESCO-weighted (1/2+ 217 keV distribution):\033[0m\n");
         printf("    mean theta_cm  = %.1f deg    mean theta_lab = %.1f deg\n", swTh / sw, swThLab / sw);
         printf("    yield-weighted predicted sigma(Ex) = %.3f MeV   ->  gap/sigma = %.2f\n", swSig / sw,
                0.115 / (swSig / sw));
      }
   }

   // ------------------------------------------------------------------------------------------
   // Figure
   // ------------------------------------------------------------------------------------------
   if (!outDir.Length())
      return;
   gSystem->mkdir(outDir, kTRUE);

   TCanvas *cv = new TCanvas("cv", "17C inelastic kinematics", 1400, 1000);
   cv->Divide(2, 2);

   const int lsty[kNL] = {1, 2, 3};

   cv->cd(1);
   TMultiGraph *mg1 = new TMultiGraph();
   for (int c = 0; c < 2; ++c)
      for (int l = 0; l < kNL; ++l) {
         gThLab[c][l]->SetLineColor(chan[c].color);
         gThLab[c][l]->SetLineStyle(lsty[l]);
         gThLab[c][l]->SetLineWidth(2);
         mg1->Add(gThLab[c][l], "L");
      }
   mg1->SetTitle("recoil lab angle;#theta_{cm} (projectile) [deg];#theta_{lab} of p / d [deg]");
   mg1->Draw("A");
   {
      TLegend *lg = new TLegend(0.55, 0.62, 0.88, 0.88);
      lg->SetBorderSize(0);
      for (int c = 0; c < 2; ++c)
         lg->AddEntry(gThLab[c][0], chan[c].label, "l");
      lg->Draw();
   }

   cv->cd(2);
   TMultiGraph *mg2 = new TMultiGraph();
   for (int c = 0; c < 2; ++c)
      for (int l = 0; l < kNL; ++l) {
         gKE[c][l]->SetLineColor(chan[c].color);
         gKE[c][l]->SetLineStyle(lsty[l]);
         gKE[c][l]->SetLineWidth(2);
         mg2->Add(gKE[c][l], "L");
      }
   mg2->SetTitle("recoil energy;#theta_{lab} [deg];KE of p / d [MeV]");
   mg2->Draw("A");

   cv->cd(3);
   TMultiGraph *mg3 = new TMultiGraph();
   for (int c = 0; c < 2; ++c) {
      gSepKE[c]->SetLineColor(chan[c].color);
      gSepKE[c]->SetLineWidth(3);
      mg3->Add(gSepKE[c], "L");
   }
   mg3->SetTitle("recoil-energy separation of the two proposal states;#theta_{cm} [deg];KE(217) - KE(332) [keV]");
   mg3->Draw("A");
   {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.035);
      tx.DrawLatex(0.15, 0.85, "the 115 keV level gap, as it appears in the MEASURED quantity");
   }

   // Panel 4: the LEVERAGE, which is the thing this macro can actually compute from kinematics
   // alone. |dEx/dKE| flat at ~0.53/0.56 is the whole reason a single number characterises the
   // requirement; |dEx/dtheta| is shown beside it to justify calling the angle term subdominant.
   // (An earlier version drew a predicted sigma(Ex) here from one angle-averaged sigma(KE). It was
   // flat, and that flatness was an artefact of the input -- see the header. Measured sigma(Ex)
   // belongs to inel_summary_C17.C, which has the per-slice sigma(KE) to build it from.)
   cv->cd(4);
   TMultiGraph *mg4 = new TMultiGraph();
   TLegend *lg4 = new TLegend(0.45, 0.66, 0.89, 0.89);
   lg4->SetBorderSize(0);
   for (int c = 0; c < 2; ++c) {
      gDerKE[c]->SetLineColor(chan[c].color);
      gDerKE[c]->SetLineWidth(3);
      gDerKE[c]->SetLineStyle(1);
      mg4->Add(gDerKE[c], "L");
      lg4->AddEntry(gDerKE[c], Form("%s  |dE_{x}/dKE|", chan[c].label), "l");
      gDerTh[c]->SetLineColor(chan[c].color);
      gDerTh[c]->SetLineWidth(2);
      gDerTh[c]->SetLineStyle(2);
      mg4->Add(gDerTh[c], "L");
      lg4->AddEntry(gDerTh[c], Form("%s  |dE_{x}/d#theta| [MeV/deg]", chan[c].label), "l");
   }
   mg4->SetTitle("leverage: how a tracking error becomes an E_{x} error;#theta_{lab} [deg];|dE_{x}/dKE|  and  |dE_{x}/d#theta|");
   mg4->SetMinimum(0.0);
   mg4->SetMaximum(0.75);
   mg4->Draw("A");
   lg4->Draw();
   {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.032);
      tx.DrawLatex(0.14, 0.20, "|dE_{x}/dKE| flat #Rightarrow the 115 keV gap is 217 keV of recoil energy,");
      tx.DrawLatex(0.14, 0.15, "at every angle. So #sigma(E_{x}) is set by #sigma(KE) alone.");
   }

   cv->SaveAs(outDir + "kinematics_C17inel.png");
   printf("\n  wrote %skinematics_C17inel.png\n\n", outDir.Data());
}
