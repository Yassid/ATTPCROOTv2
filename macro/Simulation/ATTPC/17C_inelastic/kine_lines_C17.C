/// @file kine_lines_C17.C
/// @brief The reconstructed kinematic lines: recoil KE vs theta_lab, as the experiment would see
///        it, with the three 17C levels drawn on top. Answers "are the loci separable in the
///        kinematics plot itself", not just in the derived Ex spectrum.
///
///   root -b -q 'kine_lines_C17.C'
///   root -b -q 'kine_lines_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","dd_b285")'
///
/// WHY BOTH PANELS. The raw reconstructed (theta_lab, KE) of a real event carries the beam energy
/// at ITS OWN vertex, and the beam loses 14.6 MeV crossing the chamber. Three loci 217 keV apart
/// therefore sit inside a band several MeV wide that is pure vertex position -- they are not
/// separable raw, and no amount of tracking resolution changes that.
///
/// The correction is the same one that dominates the Ex result: rebuild each event at a COMMON
/// reference beam energy. For a measured (theta_lab, KE) at vertex z, take the excitation it
/// implies at E_beam(z), then ask what KE that same excitation would give at E0 for the same
/// theta_lab. The vertex spread collapses and the lines appear.
///
///     KE_corr = KE such that  Ex(E0, theta_lab, KE_corr) = Ex(E_beam(z_reco), theta_lab, KE_reco)
///
/// This is not a cosmetic re-plot: it is exactly the analysis step the proposal needs to adopt,
/// shown in the variable an experimenter looks at first.
///
/// E_beam(z) is MEASURED from this campaign's own truth (same bisection as inel_summary_C17.C),
/// never assumed from a stopping-power table.
///
/// Events are weighted by the FRESCO angular distributions in fresco/, so the density on the plot
/// is the density the experiment would actually collect rather than the flat generator's.

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
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
double klOmega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}

/// Two-body inversion: measured ejectile (theta_lab, KE) plus a beam energy -> residual Ex.
/// Identical in form to acceptance_C14.C:acc_kine.
double klEx(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * klOmega2(s, m1 * m1, m2 * m2) * klOmega2(uu, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                         (2 * m2 * m2) +
                      s + uu - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}

/// Invert klEx in KE: the ejectile energy that leaves the residual at `ex`, at this beam energy
/// and lab angle. Ex is monotonic in KE (dEx/dKE = -0.53, flat), so bisection is safe and exact.
/// Returns -1 if the state is not kinematically reachable at this angle.
double klKE(double m1, double m2, double m3, double m4, double Ebeam, double thLabDeg, double ex, double keMax)
{
   // SCAN FIRST, then bisect. Bracketing straight on [0, keMax] does not work: klEx returns its
   // -1e9 sentinel at very low KE, where the two-body relation has no solution at this angle, so a
   // fixed lower bracket is invalid and every call bails out. Walk up the range instead and take
   // the first sign change between two VALID samples.
   const int nS = 400;
   double prevKE = -1, prevF = 0;
   for (int i = 1; i <= nS; ++i) {
      const double ke = keMax * i / nS;
      const double v = klEx(m1, m2, m3, m4, Ebeam, thLabDeg, ke);
      if (v < -1e8) { // unphysical here; reset the bracket
         prevKE = -1;
         continue;
      }
      const double fv = v - ex;
      if (prevKE > 0 && prevF * fv <= 0) {
         double lo = prevKE, hi = ke, flo = prevF;
         for (int it = 0; it < 80; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fm = klEx(m1, m2, m3, m4, Ebeam, thLabDeg, mid) - ex;
            if (flo * fm <= 0)
               hi = mid;
            else {
               lo = mid;
               flo = fm;
            }
         }
         return 0.5 * (lo + hi);
      }
      prevKE = ke;
      prevF = fv;
   }
   return -1;
}

bool klReadFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
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
double klInterp(const std::vector<double> &th, const std::vector<double> &xs, double x)
{
   if (th.empty() || x < th.front() || x > th.back())
      return 0;
   const size_t i = std::min((size_t)(x - th.front()), th.size() - 2);
   const double f = (x - th[i]) / (th[i + 1] - th[i]);
   return xs[i] * (1 - f) + xs[i + 1] * f;
}
} // namespace

void kine_lines_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString cfgTag = "pp_b285",
                    TString frescoDir = "./fresco/", Double_t chi2Cut = 5.0, TString outDir = "./plots/",
                    Double_t thZoomLo = 48.0, Double_t thZoomHi = 78.0, Double_t slLo = 58.0, Double_t slHi = 68.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};
   const char *stName[nL] = {"3/2^{+} g.s.", "1/2^{+} 217 keV", "5/2^{+} 332 keV"};
   const int col[nL] = {kBlack, kRed + 1, kBlue + 1};

   TString chan(cfgTag);
   chan.Remove(2);
   TString btag(cfgTag);
   btag.Remove(0, 3);

   const double uAmu = 931.49401;
   const double m1 = 17.0225787 * uAmu;                            // beam 17C
   const double mL = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu;  // target = ejectile
   const double m4 = m1;                                           // residual IS the beam
   const double E0 = 135.0;                                        // reference beam energy
   const double keMax = (chan == "dd") ? 60.0 : 34.0;

   printf("\n\033[1;33m===== reconstructed kinematic lines, %s =====\033[0m\n", cfgTag.Data());

   std::vector<double> fth[2], fxs[2];
   const bool haveFresco = klReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) &&
                           klReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1]);
   if (!haveFresco)
      printf("  [no FRESCO tables in %s -- events drawn with the flat generator's density]\n", frescoDir.Data());

   // ---- histograms ------------------------------------------------------------------------------
   const int nTh = 150, nKE = 400;
   TH2D *hRaw = new TH2D("hRaw", "", nTh, 15, 90, nKE, 0, keMax);
   TH2D *hCor = new TH2D("hCor", "", nTh, 15, 90, nKE, 0, keMax);
   // The residual against the GROUND-STATE kinematic line, dKE = KE - KE_gs(theta_lab). This is
   // the variable that actually shows the loci: inside any theta slice wide enough to hold
   // statistics the kinematics itself sweeps several MeV, which buries a 227 keV separation. The
   // subtraction removes exactly that sweep and leaves three horizontal bands at 0, -0.22, -0.45.
   TH2D *hRes = new TH2D("hRes", "", 100, thZoomLo, thZoomHi, 260, -0.95, 0.35);
   TH1D *hSl[nL];
   TH1D *hSlRaw = new TH1D("hSlRaw", "", 130, -0.95, 0.35);
   hSlRaw->SetDirectory(nullptr);
   hSlRaw->SetLineColor(kGray + 2);
   hSlRaw->SetLineStyle(2);
   hSlRaw->SetLineWidth(2);
   for (int l = 0; l < nL; ++l) {
      hSl[l] = new TH1D(Form("hSl%d", l), "", 130, -0.95, 0.35);
      hSl[l]->SetDirectory(nullptr);
      hSl[l]->SetLineColor(col[l]);
      hSl[l]->SetLineWidth(2);
   }

   // The analytic g.s. line at the reference energy, cached on a fine grid so the per-event
   // subtraction is a lookup rather than 400 klEx calls.
   const int nGrid = 1500;
   const double gLo = 15.0, gHi = 90.0;
   std::vector<double> keGS(nGrid + 1);
   for (int i = 0; i <= nGrid; ++i)
      keGS[i] = klKE(m1, mL, mL, m4, E0, gLo + (gHi - gLo) * i / nGrid, 0.0, keMax);
   auto keGSat = [&](double th) {
      if (th <= gLo || th >= gHi)
         return -1.0;
      const double x = (th - gLo) / (gHi - gLo) * nGrid;
      const int i = (int)x;
      if (keGS[i] < 0 || keGS[i + 1] < 0)
         return -1.0;
      return keGS[i] + (x - i) * (keGS[i + 1] - keGS[i]);
   };

   double ebz_a = 0, ebz_b = 0;
   bool haveEbz = false;
   long nUsed = 0;

   for (int l = 0; l < nL; ++l) {
      TString dirq = TString("\"") + root + "/" + cfgTag + "\"";
      TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
      TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-8s MISSING\033[0m\n", stTag[l]);
         continue;
      }
      TFile *f = TFile::Open(found);
      TTree *t = f ? (TTree *)f->Get("res") : nullptr;
      if (!t) {
         printf("\033[1;31m  no res tree in %s\033[0m\n", found.Data());
         continue;
      }
      double thT, thR, keT, keR, cmT, c2n, zT, zR;
      t->SetBranchAddress("thTrue", &thT);
      t->SetBranchAddress("thReco", &thR);
      t->SetBranchAddress("keTrue", &keT);
      t->SetBranchAddress("keReco", &keR);
      t->SetBranchAddress("cmTrue", &cmT);
      t->SetBranchAddress("chi2ndf", &c2n);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);

      // E_beam(z) from truth, on the first level present
      if (!haveEbz) {
         TGraph g;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = E0 - 40, hi = E0 + 40;
            if ((klEx(m1, mL, mL, m4, lo, thT, keT) - ExGen[l]) * (klEx(m1, mL, mL, m4, hi, thT, keT) - ExGen[l]) > 0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((klEx(m1, mL, mL, m4, lo, thT, keT) - ExGen[l]) * (klEx(m1, mL, mL, m4, mid, thT, keT) - ExGen[l]) <=
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
            printf("  E_beam(z) = %.4f %+.6f * z[mm], reference E0 = %.1f MeV\n", ebz_a, ebz_b, E0);
         }
      }

      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (c2n >= chi2Cut)
            continue;
         double w = 1.0;
         if (haveFresco && l > 0)
            w = klInterp(fth[l - 1], fxs[l - 1], cmT);
         else if (haveFresco && l == 0)
            w = klInterp(fth[0], fxs[0], cmT); // no elastic calculation; 1/2+ shape stands in
         if (w <= 0)
            continue;

         hRaw->Fill(thR, keR, w);

         // vertex correction: the excitation this event implies at ITS OWN beam energy, then the
         // KE that same excitation would give at the reference energy and the same lab angle.
         if (!haveEbz)
            continue;
         const double exHere = klEx(m1, mL, mL, m4, ebz_a + ebz_b * zR, thR, keR);
         if (exHere < -1e8)
            continue;
         const double keCorr = klKE(m1, mL, mL, m4, E0, thR, exHere, keMax);
         if (keCorr <= 0)
            continue;
         hCor->Fill(thR, keCorr, w);
         const double kgs = keGSat(thR);
         if (kgs > 0) {
            if (thR >= thZoomLo && thR <= thZoomHi)
               hRes->Fill(thR, keCorr - kgs, w);
            if (thR >= slLo && thR < slHi) {
               hSl[l]->Fill(keCorr - kgs, w);
               hSlRaw->Fill(keR - kgs, w); // the same residual WITHOUT the vertex correction
            }
         }
         ++nUsed;
      }
      f->Close();
   }
   printf("  %ld events on the corrected plot\n", nUsed);

   // ---- analytic lines ---------------------------------------------------------------------------
   TGraph *gLine[nL];
   for (int l = 0; l < nL; ++l) {
      gLine[l] = new TGraph();
      for (double th = 16; th <= 89; th += 0.25) {
         const double k = klKE(m1, mL, mL, m4, E0, th, ExGen[l], keMax);
         if (k > 0)
            gLine[l]->SetPoint(gLine[l]->GetN(), th, k);
      }
      gLine[l]->SetLineColor(col[l]);
      gLine[l]->SetLineWidth(2);
      gLine[l]->SetLineStyle(l == 0 ? 1 : (l == 1 ? 7 : 3));
   }

   // How far apart the lines are, in the measured variable
   printf("\n  %10s %12s %12s %12s   %s\n", "theta_lab", "KE(gs)", "KE(217)", "KE(332)", "gap 217-332 [keV]");
   for (double th : {30., 40., 50., 55., 60., 65., 70., 75., 80.}) {
      const double k0 = klKE(m1, mL, mL, m4, E0, th, ExGen[0], keMax);
      const double k1 = klKE(m1, mL, mL, m4, E0, th, ExGen[1], keMax);
      const double k2 = klKE(m1, mL, mL, m4, E0, th, ExGen[2], keMax);
      if (k1 < 0 || k2 < 0)
         continue;
      printf("  %10.0f %12.3f %12.3f %12.3f   %14.0f\n", th, k0, k1, k2, (k1 - k2) * 1000);
   }

   // ---- figure -------------------------------------------------------------------------------------
   if (!outDir.Length())
      return;
   gSystem->mkdir(outDir, kTRUE);
   TCanvas *cv = new TCanvas("cvkl", "kinematic lines", 1500, 1050);
   cv->Divide(2, 2);

   auto drawLines = [&](bool leg) {
      for (int l = 0; l < nL; ++l)
         gLine[l]->Draw("L SAME");
      if (leg) {
         TLegend *lg = new TLegend(0.58, 0.66, 0.90, 0.89);
         lg->SetBorderSize(0);
         lg->SetFillStyle(0);
         for (int l = 0; l < nL; ++l)
            lg->AddEntry(gLine[l], stName[l], "l");
         lg->Draw();
      }
   };

   cv->cd(1);
   gPad->SetLogz();
   hRaw->SetTitle(Form("%s  RAW reconstruction  (each event at its own vertex beam energy);"
                       "#theta_{lab} [deg];KE of the recoil [MeV]",
                       cfgTag.Data()));
   hRaw->Draw("COLZ");
   drawLines(true);

   cv->cd(2);
   gPad->SetLogz();
   hCor->SetTitle(Form("%s  VERTEX-CORRECTED to E_{beam} = %.0f MeV;#theta_{lab} [deg];KE of the recoil [MeV]",
                       cfgTag.Data(), E0));
   hCor->Draw("COLZ");
   drawLines(false);

   cv->cd(3);
   gPad->SetLogz();
   hRes->SetTitle("vertex-corrected, g.s. line subtracted;#theta_{lab} [deg];"
                  "KE #minus KE_{g.s.}(#theta_{lab}) [MeV]");
   hRes->Draw("COLZ");
   {
      // the three levels are now horizontal bands; draw where each should sit
      TGraph *gb[nL];
      TLegend *lg = new TLegend(0.15, 0.16, 0.52, 0.38);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      for (int l = 0; l < nL; ++l) {
         gb[l] = new TGraph();
         for (double th = thZoomLo; th <= thZoomHi; th += 0.5) {
            const double k = klKE(m1, mL, mL, m4, E0, th, ExGen[l], keMax);
            const double k0 = klKE(m1, mL, mL, m4, E0, th, 0.0, keMax);
            if (k > 0 && k0 > 0)
               gb[l]->SetPoint(gb[l]->GetN(), th, k - k0);
         }
         gb[l]->SetLineColor(col[l]);
         gb[l]->SetLineWidth(2);
         gb[l]->SetLineStyle(l == 0 ? 1 : (l == 1 ? 7 : 3));
         gb[l]->Draw("L SAME");
         lg->AddEntry(gb[l], stName[l], "l");
      }
      lg->Draw();
   }

   cv->cd(4);
   double mx = 0;
   for (int l = 0; l < nL; ++l)
      mx = std::max(mx, hSl[l]->GetMaximum());
   if (hSlRaw->Integral() > 0)
      hSlRaw->Scale(mx / std::max(1e-9, hSlRaw->GetMaximum()));
   for (int l = 0; l < nL; ++l) {
      hSl[l]->SetMaximum(1.35 * mx);
      hSl[l]->SetTitle(Form("%.0f < #theta_{lab} < %.0f deg;KE #minus KE_{g.s.}(#theta_{lab}) [MeV];yield [arb]",
                            slLo, slHi));
      hSl[l]->Draw(l ? "HIST SAME" : "HIST");
   }
   hSlRaw->Draw("HIST SAME");
   {
      TLegend *lg = new TLegend(0.14, 0.62, 0.52, 0.89);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      for (int l = 0; l < nL; ++l)
         lg->AddEntry(hSl[l], stName[l], "l");
      lg->AddEntry(hSlRaw, "sum, NO vertex correction", "l");
      lg->Draw();
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.032);
      const double th = 0.5 * (slLo + slHi);
      const double k1 = klKE(m1, mL, mL, m4, E0, th, ExGen[1], keMax);
      const double k2 = klKE(m1, mL, mL, m4, E0, th, ExGen[2], keMax);
      tx.DrawLatex(0.56, 0.52, Form("115 keV level gap ="));
      tx.DrawLatex(0.56, 0.47, Form("%.0f keV of recoil energy", (k1 - k2) * 1000));
   }

   cv->SaveAs(outDir + "kinelines_" + cfgTag + ".png");
   printf("\n  wrote %skinelines_%s.png\n\n", outDir.Data(), cfgTag.Data());
}
