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
/// @param R  N(elastic)/N(1/2+ 217), the elastic normalisation. THIS IS NOT MEASURED. No elastic
///           FRESCO calculation was supplied with the proposal, so the elastic events are given the
///           1/2+ angular SHAPE and their total is set by this number. An earlier version of this
///           macro left the elastic weighted by the 1/2+ cross section outright, i.e. R = 1, which
///           drew a comfortable three-peak spectrum and badly understated the problem: elastic
///           scattering is 10-100x the inelastic here, and the two states of interest are structure
///           on its low-energy flank. Vary R until the elastic curve exists.
/// @param slices  theta_lab window edges. The separation is strongly angle-dependent -- sigma(KE)
///           varies twentyfold across this range -- so one window cannot represent the measurement
///           and the panels are the point of the figure.
void kine_lines_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString cfgTag = "pp_b285",
                    TString frescoDir = "./fresco/", Double_t chi2Cut = 5.0, TString outDir = "./plots/",
                    Double_t R = 10.0, Double_t resLo = -1.5, Double_t resHi = 0.7)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};
   const char *stName[nL] = {"3/2^{+} g.s. (elastic)", "1/2^{+} 217 keV", "5/2^{+} 332 keV"};
   const int col[nL] = {kGray + 2, kRed + 1, kBlue + 1};

   // five theta_lab windows, forward to backward
   const int nW = 5;
   const double wLo[nW] = {30, 40, 50, 60, 70};
   const double wHi[nW] = {40, 50, 60, 70, 80};

   TString chan(cfgTag);
   chan.Remove(2);
   TString btag(cfgTag);
   btag.Remove(0, 3);

   const double uAmu = 931.49401;
   const double m1 = 17.0225787 * uAmu;                           // beam 17C
   const double mL = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu; // target = ejectile
   const double m4 = m1;                                          // residual IS the beam
   const double E0 = 135.0;                                       // reference beam energy
   const double keMax = (chan == "dd") ? 60.0 : 34.0;

   printf("\n\033[1;33m===== reconstructed kinematic lines, %s (N_elastic/N_217 = %.0f) =====\033[0m\n",
          cfgTag.Data(), R);

   std::vector<double> fth[2], fxs[2];
   const bool haveFresco = klReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) &&
                           klReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1]);
   if (!haveFresco) {
      printf("\033[1;31m  no FRESCO tables in %s -- cannot weight the events. Aborting.\033[0m\n", frescoDir.Data());
      return;
   }

   // ---- the analytic g.s. line, cached, so the per-event subtraction is a lookup ---------------
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

   // ---- histograms -----------------------------------------------------------------------------
   const int nB = (int)std::lround((resHi - resLo) / 0.025); // 25 keV bins
   TH2D *hCor = new TH2D("hCor", "", 150, 15, 90, 400, 0, keMax);
   TH2D *hRes = new TH2D("hRes", "", 130, 20, 85, 200, resLo, resHi);
   TH1D *hW[nW][nL];
   TH1D *hWraw[nW]; // the same window with NO vertex correction, for the contrast
   std::vector<double> resW[nW][nL];
   for (int w = 0; w < nW; ++w) {
      hWraw[w] = new TH1D(Form("hWraw%d", w), "", nB, resLo, resHi);
      hWraw[w]->SetDirectory(nullptr);
      hWraw[w]->SetLineColor(kGray + 1);
      hWraw[w]->SetLineStyle(2);
      hWraw[w]->SetLineWidth(2);
      for (int l = 0; l < nL; ++l) {
         hW[w][l] = new TH1D(Form("hW%d_%d", w, l), "", nB, resLo, resHi);
         hW[w][l]->SetDirectory(nullptr);
         hW[w][l]->SetLineColor(col[l]);
         hW[w][l]->SetLineWidth(2);
      }
   }

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
         return;
      }
      TFile *f = TFile::Open(found);
      TTree *t = f ? (TTree *)f->Get("res") : nullptr;
      if (!t) {
         printf("\033[1;31m  no res tree in %s\033[0m\n", found.Data());
         return;
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
         if (c2n >= chi2Cut || !haveEbz)
            continue;
         // ANGULAR WEIGHT. The 1/2+ and 5/2+ carry their own FRESCO dsigma/dOmega in absolute
         // mb/sr, so their relative areas ARE the cross-section ratio (2.38). The elastic has no
         // calculation: it borrows the 1/2+ shape and its total is fixed to R x N(217) below, which
         // is an assumption and is labelled as one on the figure.
         const double w = klInterp(fth[l == 0 ? 0 : l - 1], fxs[l == 0 ? 0 : l - 1], cmT);
         if (w <= 0)
            continue;

         const double exHere = klEx(m1, mL, mL, m4, ebz_a + ebz_b * zR, thR, keR);
         if (exHere < -1e8)
            continue;
         const double keCorr = klKE(m1, mL, mL, m4, E0, thR, exHere, keMax);
         if (keCorr <= 0)
            continue;
         hCor->Fill(thR, keCorr, w);

         const double kgs = keGSat(thR);
         if (kgs <= 0)
            continue;
         hRes->Fill(thR, keCorr - kgs, w);
         for (int wi = 0; wi < nW; ++wi)
            if (thR >= wLo[wi] && thR < wHi[wi]) {
               hW[wi][l]->Fill(keCorr - kgs, w);
               hWraw[wi]->Fill(keR - kgs, w);
               resW[wi][l].push_back(keCorr - kgs);
            }
         ++nUsed;
      }
      f->Close();
   }
   printf("  %ld events on the corrected plot\n", nUsed);

   // ---- put the elastic on its stated footing ---------------------------------------------------
   // Weights alone give 217:332 = the cross-section ratio, which is right. The elastic was filled
   // with the 1/2+ shape, so its normalisation carries no meaning until it is set here.
   for (int w = 0; w < nW; ++w) {
      const double i217 = hW[w][1]->Integral();
      const double igs = hW[w][0]->Integral();
      if (igs > 0 && i217 > 0)
         hW[w][0]->Scale(R * i217 / igs);
   }

   // ---- analytic lines, and the numbers behind the panels ----------------------------------------
   TGraph *gLine[nL], *gBand[nL];
   for (int l = 0; l < nL; ++l) {
      gLine[l] = new TGraph();
      gBand[l] = new TGraph();
      for (double th = 16; th <= 89; th += 0.25) {
         const double k = klKE(m1, mL, mL, m4, E0, th, ExGen[l], keMax);
         const double k0 = klKE(m1, mL, mL, m4, E0, th, 0.0, keMax);
         if (k > 0)
            gLine[l]->SetPoint(gLine[l]->GetN(), th, k);
         if (k > 0 && k0 > 0 && th >= 20 && th <= 85)
            gBand[l]->SetPoint(gBand[l]->GetN(), th, k - k0);
      }
      for (TGraph *g : {gLine[l], gBand[l]}) {
         g->SetLineColor(l == 0 ? kBlack : col[l]);
         g->SetLineWidth(2);
         g->SetLineStyle(l == 0 ? 1 : (l == 1 ? 7 : 3));
      }
   }

   auto iqrSigma = [](std::vector<double> v) {
      if (v.size() < 20)
         return -1.0;
      std::sort(v.begin(), v.end());
      return (v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())]) / 1.349;
   };

   printf("\n  %-12s %10s %10s %10s %12s %12s\n", "theta_lab", "KE(gs)", "gap 0-217", "gap 217-332", "sigma(KE)",
          "separation");
   printf("  %-12s %10s %10s %10s %12s %12s\n", "[deg]", "[MeV]", "[keV]", "[keV]", "[MeV]", "217/332");
   double sepW[nW];
   for (int w = 0; w < nW; ++w) {
      const double th = 0.5 * (wLo[w] + wHi[w]);
      const double k0 = klKE(m1, mL, mL, m4, E0, th, ExGen[0], keMax);
      const double k1 = klKE(m1, mL, mL, m4, E0, th, ExGen[1], keMax);
      const double k2 = klKE(m1, mL, mL, m4, E0, th, ExGen[2], keMax);
      const double s1 = iqrSigma(resW[w][1]), s2 = iqrSigma(resW[w][2]);
      sepW[w] = (s1 > 0 && s2 > 0) ? (k1 - k2) / (s1 + s2) : -1;
      printf("  %5.0f-%-6.0f %10.3f %10.0f %10.0f %12.3f %12.2f\n", wLo[w], wHi[w], k0, (k0 - k1) * 1000,
             (k1 - k2) * 1000, s1, sepW[w]);
   }

   // ---- figure ------------------------------------------------------------------------------------
   if (!outDir.Length())
      return;
   gSystem->mkdir(outDir, kTRUE);
   // 4x2: the two maps, the five angular windows, and a summary pad. On a 3x2 the fifth
   // window silently landed on a pad that does not exist and overwrote the fourth.
   TCanvas *cv = new TCanvas("cvkl", "kinematic lines", 2000, 1000);
   cv->Divide(4, 2);

   cv->cd(1);
   gPad->SetLogz();
   hCor->SetTitle(Form("%s  vertex-corrected to E_{beam} = %.0f MeV;#theta_{lab} [deg];KE of the recoil [MeV]",
                       cfgTag.Data(), E0));
   hCor->Draw("COLZ");
   for (int l = 0; l < nL; ++l)
      gLine[l]->Draw("L SAME");
   {
      TLegend *lg = new TLegend(0.52, 0.68, 0.90, 0.89);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      for (int l = 0; l < nL; ++l)
         lg->AddEntry(gLine[l], stName[l], "l");
      lg->Draw();
   }

   cv->cd(2);
   gPad->SetLogz();
   hRes->SetTitle("g.s. line subtracted -- the loci as flat bands;#theta_{lab} [deg];"
                  "KE #minus KE_{g.s.} [MeV]");
   hRes->GetYaxis()->SetTitleOffset(1.15);
   hRes->Draw("COLZ");
   for (int l = 0; l < nL; ++l)
      gBand[l]->Draw("L SAME");
   {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.030);
      tx.DrawLatex(0.14, 0.20, "the five windows below are slices of this panel");
   }

   // the five angular windows
   for (int w = 0; w < nW; ++w) {
      cv->cd(3 + w);
      // LOG y. At a realistic elastic ratio the components span two decades, and on a linear axis
      // the 1/2+ and 5/2+ are invisible under the elastic -- which reads as "there is nothing
      // there" when the fit in decompose_C17.C recovers them to 5.6 % and 9.8 %. The extraction
      // works off the flank, not off a visible bump, and only a log axis shows the flank.
      gPad->SetLogy();
      TH1D *sum = (TH1D *)hW[w][0]->Clone(Form("sum%d", w));
      sum->SetDirectory(nullptr);
      for (int l = 1; l < nL; ++l)
         sum->Add(hW[w][l]);
      sum->SetLineColor(kBlack);
      sum->SetLineWidth(3);
      sum->SetTitle(Form("%.0f < #theta_{lab} < %.0f deg   (sep %.2f);"
                         "KE #minus KE_{g.s.}(#theta_{lab}) [MeV];counts / 25 keV",
                         wLo[w], wHi[w], sepW[w]));
      sum->SetMaximum(4.0 * sum->GetMaximum());
      sum->SetMinimum(0.5);
      sum->Draw("HIST");
      for (int l = 0; l < nL; ++l)
         hW[w][l]->Draw("HIST SAME");
      if (hWraw[w]->Integral() > 0) {
         // Same INTEGRAL as the sum, not the same peak. Peak-matching on a log axis puts this
         // curve above everything and it reads as a dominant component; it is only a shape, drawn
         // to show how much wider the spectrum is without the vertex correction.
         hWraw[w]->Scale(sum->Integral() / hWraw[w]->Integral());
         hWraw[w]->Draw("HIST SAME");
      }
      if (w == 0) {
         TLegend *lg = new TLegend(0.13, 0.60, 0.55, 0.89);
         lg->SetBorderSize(0);
         lg->SetFillStyle(0);
         lg->AddEntry(sum, "sum = what is measured", "l");
         for (int l = 0; l < nL; ++l)
            lg->AddEntry(hW[w][l], stName[l], "l");
         lg->AddEntry(hWraw[w], "sum, NO vertex corr. (shape)", "l");
         lg->Draw();
      }
   }

   cv->cd(8);
   {
      TLatex tx;
      tx.SetTextSize(0.055);
      tx.SetTextFont(102);
      tx.DrawLatexNDC(0.04, 0.92, Form("%s", cfgTag.Data()));
      tx.SetTextSize(0.048);
      tx.DrawLatexNDC(0.04, 0.83, "#theta_{lab}    #sigma(KE)   sep(217/332)");
      for (int w = 0; w < nW; ++w)
         tx.DrawLatexNDC(0.04, 0.74 - 0.09 * w,
                         Form("%2.0f-%2.0f    %5.3f       %4.2f", wLo[w], wHi[w], iqrSigma(resW[w][1]), sepW[w]));
      tx.SetTextFont(42);
      tx.SetTextSize(0.042);
      tx.DrawLatexNDC(0.04, 0.22, "separation = #Delta / (#sigma_{1}+#sigma_{2})");
      tx.DrawLatexNDC(0.04, 0.14, "best at #theta_{lab} 50-70 deg,");
      tx.DrawLatexNDC(0.04, 0.07, "which is where the yield sits");
   }

   cv->cd(0);
   {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.017);
      tx.DrawLatex(0.005, 0.005,
                   Form("elastic normalisation ASSUMED: N(elastic)/N(1/2+) = %.0f, with the 1/2+ angular shape "
                        "-- no elastic FRESCO calculation exists. 1/2+ : 5/2+ = 2.38 is the real cross-section "
                        "ratio.",
                        R));
   }

   cv->SaveAs(outDir + "kinelines_" + cfgTag + ".png");
   printf("\n  wrote %skinelines_%s.png\n\n", outDir.Data(), cfgTag.Data());
}
