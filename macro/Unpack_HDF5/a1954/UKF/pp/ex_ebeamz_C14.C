/// @file ex_ebeamz_C14.C
/// @brief Re-reconstruct the a1954 14C(p,p') excitation energy with the beam energy taken at the
/// RECONSTRUCTED VERTEX instead of one constant for the whole chamber, and plot the difference.
///
/// WHY. `ex_C14.C` uses a single beam energy for every event. The 14C beam loses ~11 MeV crossing
/// the metre of H2 at 300 torr and the vertex is spread over the whole drift, so that one number
/// is wrong by several MeV in a way perfectly correlated with a quantity the fit already
/// measures. On simulation this is the DOMINANT term in the Ex width: rebuilding Ex from perfect
/// truth kinematics under the constant-beam assumption still gives sigma = 0.25 MeV, the same as
/// the fully reconstructed value, and no pad pitch or magnetic field changes it. See
/// macro/Simulation/ATTPC/14C_pp/highfield/RESULTS.md.
///
/// THE PROFILE COMES FROM CATIMA, NOT FROM A FIT. On the simulation the profile could be
/// extracted from MC truth; on data it cannot, so E(z) is integrated here from the CATIMA
/// stopping power of 14C in H2 at the a1954 density. That is the one step this analysis needed
/// that the simulation study did not.
///
/// THE ABSOLUTE SCALE IS PRESERVED ON PURPOSE. E0 is set so that E(z) equals the adopted constant
/// beam energy at the MEAN VERTEX of the elastic events. The correction is then a pure
/// de-trending: it removes the z-dependence and leaves the calibration -- which was anchored on
/// known levels -- where it was. Any residual shift of the peaks after that is a result, not an
/// input.
///
/// Z HANDEDNESS IS CHECKED, NOT ASSUMED. The fitters work in a lab frame with z = ZPadPlane -
/// z_digi, i.e. z = 0 at the entrance window and z = 1000 mm at the pad plane, so the beam energy
/// FALLS with increasing z. The macro applies both signs and reports which one flattens Ex vs z,
/// so if the data convention were the other way round the printout says so instead of silently
/// doubling the drift.
///
///   root -b -q 'ex_ebeamz_C14.C("plots/proton_kin_cat5_s013.root",159.75,"s013")'

#include <tuple>
#include <vector>

static double ez_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics, verbatim from pp/ex_C14.C: returns {Ex, theta_cm [deg]}
static std::tuple<double, double> ez_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * ez_om2(s, m1 * m1, m2 * m2) * ez_om2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (ez_om2(s, m1 * m1, m2 * m2) * ez_om2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

/// the recoil-proton kinetic energy that puts a given excitation at a given lab angle -- found by
/// bisection on the SAME expression the analysis inverts, so the locus and the data cannot
/// disagree through a different kinematics convention
static double ez_locus(double m1, double m2, double m3, double m4gs, double Ebeam, double thRad, double exWanted)
{
   auto F = [&](double ke) { return std::get<0>(ez_kine(m1, m2, m3, m4gs, Ebeam, thRad, ke)) - exWanted; };
   // Walk the upper bracket down instead of assuming one: above the kinematic limit the
   // expression under the square root goes negative and F is NaN, and a fixed hi silently
   // returned "no solution" for every angle the first time this ran.
   double lo = 0.02, hi = 60.0;
   double flo = F(lo);
   if (std::isnan(flo))
      return -1;
   double fhi = F(hi);
   while ((std::isnan(fhi) || fhi * flo > 0) && hi > lo + 0.05) {
      hi *= 0.9;
      fhi = F(hi);
   }
   if (std::isnan(fhi) || fhi * flo > 0)
      return -1;
   for (int i = 0; i < 80; ++i) {
      double mid = 0.5 * (lo + hi), fm = F(mid);
      if (std::isnan(fm))
         return -1;
      if (fm * flo <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
   }
   return 0.5 * (lo + hi);
}

/// @param cache      a proton_kin_*.root written by ex_C14.C / apply_theta_corr_C14.C
/// @param Ebeam      the constant beam energy that cache was built with (159.75 for cat5_*)
/// @param driftLen   mm; z = 0 at the entrance window, z = driftLen at the pad plane
void ex_ebeamz_C14(TString cache = "plots/proton_kin_cat5_s013.root", Double_t Ebeam = 159.75, TString tag = "s013",
                   Double_t density = 3.308e-5, Double_t driftLen = 1000.0, Double_t zLo = -50., Double_t zHi = 1050.)
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(__FILE__);
   if (!cache.BeginsWith("/"))
      cache = here + "/" + cache;
   TString plotDir = here + "/plots/";

   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;

   // ---- 1. the beam energy profile, from CATIMA ------------------------------------------------
   auto *eloss = new AtTools::AtELossCATIMA(density); // g/cm3 -- NOT the mg/cm3 of AtDigiPar
   eloss->SetProjectile(14, 6, 14.003242);
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen
   eloss->SetMaterial(mat);

   // integrate dE/dx forward from an arbitrary entrance energy; the absolute level is fixed
   // afterwards by the anchor, only the SHAPE matters here
   const int nStep = 2000;
   const double dz = driftLen / nStep; // mm
   std::vector<double> zg(nStep + 1), dEg(nStep + 1);
   double E = Ebeam + 20.0; // start high; the loss curve is flat enough that the shape is stable
   double acc = 0;
   for (int i = 0; i <= nStep; ++i) {
      zg[i] = i * dz;
      dEg[i] = acc;
      double dedx = eloss->GetdEdx(E); // MeV/mm at this density
      E -= dedx * dz;
      acc += dedx * dz;
   }
   TGraph gLoss(nStep + 1, zg.data(), dEg.data());
   auto lossAt = [&](double z) {
      double zz = std::max(0.0, std::min(driftLen, z));
      return gLoss.Eval(zz);
   };
   printf("\nCATIMA 14C in H2 at %.4e g/cm3: loss over %.0f mm = %.2f MeV "
          "(%.2f at 250 mm, %.2f at 500, %.2f at 750)\n",
          density, driftLen, lossAt(driftLen), lossAt(250), lossAt(500), lossAt(750));

   // ---- 2. the cache -----------------------------------------------------------------------------
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   auto *t = (TNtuple *)f->Get("pk");
   if (!t) {
      printf("\033[1;31mno pk ntuple in %s\033[0m\n", cache.Data());
      return;
   }
   float ke, th, vz, thcm, ex, c2n;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("thcm", &thcm);
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("chi2ndf", &c2n);
   const Long64_t N = t->GetEntries();

   // anchor: mean vertex of the elastic events, so the correction is a pure de-trending
   double zsum = 0;
   long nEl = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (std::fabs(ex) < 1.0 && vz > 0 && vz < driftLen) { zsum += vz; ++nEl; }
   }
   const double zRef = nEl ? zsum / nEl : 0.5 * driftLen;
   printf("anchor: %ld elastic events, mean vertex z = %.1f mm -> E(z_ref) held at %.2f MeV\n", nEl, zRef, Ebeam);

   // ---- 3. both handednesses, to settle the frame rather than assume it -------------------------
   // sign +1: beam enters at z = 0, energy FALLS with z (the ZPadPlane - z_digi lab frame)
   // sign -1: the mirror image
   auto Ebeam_at = [&](double z, int sign) {
      double l = (sign > 0) ? lossAt(z) : lossAt(driftLen - z);
      double lref = (sign > 0) ? lossAt(zRef) : lossAt(driftLen - zRef);
      return Ebeam + (lref - l);
   };

   double slope[3] = {0, 0, 0}; // [0] uncorrected, [1] sign +1, [2] sign -1
   for (int mode = 0; mode < 3; ++mode) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      long n = 0;
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (vz < 0 || vz > driftLen)
            continue;
         double e;
         if (mode == 0)
            e = ex;
         else {
            auto [ee, cc] = ez_kine(m_C14, m_p, m_p, m_C14, Ebeam_at(vz, mode == 1 ? +1 : -1), th * TMath::DegToRad(), ke);
            if (std::isnan(ee)) continue;
            e = ee;
         }
         if (std::fabs(e) > 1.5)
            continue; // elastic only: a discrete state must be flat in z
         sx += vz; sy += e; sxx += vz * vz; sxy += vz * e; ++n;
      }
      slope[mode] = (n * sxy - sx * sy) / (n * sxx - sx * sx);
      printf("  %-28s elastic Ex-vs-z slope = %+.5f MeV/mm  (%+.3f MeV across the drift, n=%ld)\n",
             mode == 0 ? "uncorrected" : (mode == 1 ? "E(z), beam enters at z=0" : "E(z), mirrored"), slope[mode],
             slope[mode] * driftLen, n);
   }
   const int useSign = (std::fabs(slope[1]) <= std::fabs(slope[2])) ? +1 : -1;
   printf("  -> using sign %+d (%s)\n\n", useSign,
          useSign > 0 ? "beam enters at z = 0, as the fitter lab frame says" : "MIRRORED -- check the z convention");

   // ---- 4. rebuild, plot, and write a new cache --------------------------------------------------
   TFile fo((plotDir + "proton_kin_" + tag + "_ez.root").Data(), "RECREATE");
   TNtuple *out = new TNtuple("pk", "proton kinematics (E_beam at the vertex)",
                              "ke:theta:vertexz:thcm:ex:chi2ndf:exold:ebeam");

   const double exLo = -3, exHi = 14;
   auto *hOld = new TH1F("hOld", ";E_{x} [MeV];counts / 50 keV", (int)((exHi - exLo) / 0.05), exLo, exHi);
   auto *hNew = new TH1F("hNew", ";E_{x} [MeV];counts / 50 keV", (int)((exHi - exLo) / 0.05), exLo, exHi);
   auto *hVzOld = new TH2F("hVzOld", "as analysed;vertex z [mm];E_{x} [MeV]", 55, zLo, zHi, 170, -3, 14);
   auto *hVzNew = new TH2F("hVzNew", "E_{beam}(z_{vertex});vertex z [mm];E_{x} [MeV]", 55, zLo, zHi, 170, -3, 14);
   auto *hCmOld = new TH2F("hCmOld", "as analysed;#theta_{cm} [deg];E_{x} [MeV]", 90, 0, 180, 170, -3, 14);
   auto *hCmNew = new TH2F("hCmNew", "E_{beam}(z_{vertex});#theta_{cm} [deg];E_{x} [MeV]", 90, 0, 180, 170, -3, 14);
   auto *hKeTh = new TH2F("hKeTh", ";#theta_{lab} [deg];proton KE [MeV]", 170, 10, 180, 200, 0, 30);
   // The inelastic group RESTRICTED to theta_cm > cmSens, where dEx/dE_beam is an order of
   // magnitude larger than it is under the elastic peak. Integrated over all angles the effect is
   // diluted by the forward-CM events that carry most of the yield and none of the sensitivity,
   // so the integrated spectrum is the wrong place to look for it.
   const double cmSens = 50.0;
   auto *hOldI = new TH1F("hOldI", ";E_{x} [MeV];counts / 50 keV", 120, 4.5, 10.5);
   auto *hNewI = new TH1F("hNewI", ";E_{x} [MeV];counts / 50 keV", 120, 4.5, 10.5);

   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      hOld->Fill(ex);
      hVzOld->Fill(vz, ex);
      hCmOld->Fill(thcm, ex);
      hKeTh->Fill(th, ke);
      if (vz < 0 || vz > driftLen)
         continue;
      double eb = Ebeam_at(vz, useSign);
      auto [exN, cmN] = ez_kine(m_C14, m_p, m_p, m_C14, eb, th * TMath::DegToRad(), ke);
      if (std::isnan(exN))
         continue;
      hNew->Fill(exN);
      if (thcm > cmSens) {
         hOldI->Fill(ex);
         hNewI->Fill(exN);
      }
      hVzNew->Fill(vz, exN);
      hCmNew->Fill(cmN, exN);
      out->Fill(ke, th, vz, cmN, exN, c2n, ex, eb);
   }

   // ---- 5. what it did to the peaks ---------------------------------------------------------------
   auto peak = [&](TH1F *h, double lo, double hi, const char *what) {
      int b1 = h->FindBin(lo), b2 = h->FindBin(hi);
      double best = 0, bc = 0;
      for (int b = b1; b <= b2; ++b)
         if (h->GetBinContent(b) > best) { best = h->GetBinContent(b); bc = h->GetBinCenter(b); }
      TF1 g("g", "gaus", bc - 0.6, bc + 0.6);
      h->Fit(&g, "QNR");
      printf("  %-22s %-8s centroid %+7.3f   sigma %6.3f   FWHM %6.3f MeV\n", h->GetName(), what, g.GetParameter(1),
             std::fabs(g.GetParameter(2)), 2.355 * std::fabs(g.GetParameter(2)));
      return std::make_pair(g.GetParameter(1), std::fabs(g.GetParameter(2)));
   };
   printf("gaussian fits to the dominant peaks:\n");
   auto pOld = peak(hOld, -1.5, 1.5, "g.s.");
   auto pNew = peak(hNew, -1.5, 1.5, "g.s.");
   peak(hOld, 5.6, 6.6, "6.094");
   peak(hNew, 5.6, 6.6, "6.094");
   printf("\n  g.s. width %.3f -> %.3f MeV sigma  (%.0f %% narrower)\n", pOld.second, pNew.second,
          100 * (1 - pNew.second / pOld.second));

   // ---- 5b. WHERE the correction acts ------------------------------------------------------------
   // dEx/dE_beam is not a constant: it runs from ~0.12 at theta_cm 140 deg to ~0.001 at 12 deg,
   // so the same 10 MeV of beam-energy swing is worth 1.2 MeV at the backward-CM end and 15 keV
   // at the forward one. The a1954 acceptance sits almost entirely at the insensitive end, which
   // is why the integrated spectrum barely moves -- and why the angular distributions, which run
   // out to theta_cm ~115 deg, are a different matter.
   {
      const int nS = 6;
      const double cmLo[nS] = {10, 30, 50, 70, 90, 110};
      const double cmHi[nS] = {30, 50, 70, 90, 110, 140};
      std::vector<double> dOld[nS], dNew[nS], sens[nS];
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (vz < 0 || vz > driftLen)
            continue;
         double eb = Ebeam_at(vz, useSign);
         auto [exN, cmN] = ez_kine(m_C14, m_p, m_p, m_C14, eb, th * TMath::DegToRad(), ke);
         if (std::isnan(exN))
            continue;
         double a = std::get<0>(ez_kine(m_C14, m_p, m_p, m_C14, Ebeam + 0.5, th * TMath::DegToRad(), ke));
         double b = std::get<0>(ez_kine(m_C14, m_p, m_p, m_C14, Ebeam - 0.5, th * TMath::DegToRad(), ke));
         for (int s2 = 0; s2 < nS; ++s2)
            if (thcm >= cmLo[s2] && thcm < cmHi[s2]) {
               if (std::fabs(ex) < 1.5) dOld[s2].push_back(ex);
               if (std::fabs(exN) < 1.5) dNew[s2].push_back(exN);
               if (!std::isnan(a) && !std::isnan(b)) sens[s2].push_back(a - b);
            }
      }
      auto med = [](std::vector<double> v, double q) {
         if (v.size() < 20) return std::numeric_limits<double>::quiet_NaN();
         size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, q * (v.size() - 1)));
         std::nth_element(v.begin(), v.begin() + k, v.end());
         return v[k];
      };
      printf("\nelastic peak per theta_cm slice (median and IQR/1.349, |Ex| < 1.5 MeV):\n");
      printf("  %-10s %8s %9s %9s %9s %9s %12s\n", "theta_cm", "n", "med old", "sig old", "med new", "sig new",
             "dEx/dEbeam");
      for (int s2 = 0; s2 < nS; ++s2) {
         if (dOld[s2].size() < 20) {
            printf("  %3.0f-%-6.0f %8zu %9s %9s %9s %9s %12s\n", cmLo[s2], cmHi[s2], dOld[s2].size(), "-", "-", "-",
                   "-", "-");
            continue;
         }
         double so = (med(dOld[s2], .75) - med(dOld[s2], .25)) / 1.349;
         double sn = (med(dNew[s2], .75) - med(dNew[s2], .25)) / 1.349;
         printf("  %3.0f-%-6.0f %8zu %+9.3f %9.3f %+9.3f %9.3f %12.4f\n", cmLo[s2], cmHi[s2], dOld[s2].size(),
                med(dOld[s2], .5), so, med(dNew[s2], .5), sn, med(sens[s2], .5));
      }
   }

   // ---- 6. figures ---------------------------------------------------------------------------------
   // (a) kinematics: the measured plane with the two-body loci of the levels that matter
   const double exLoci[5] = {0.0, 6.094, 6.728, 7.012, 8.317};
   const int lc[5] = {kBlack, kAzure + 2, kGreen + 2, kOrange + 7, kMagenta + 1};
   auto *cK = new TCanvas("cK", "kinematics", 1500, 560);
   cK->Divide(3, 1);
   cK->cd(1);
   gPad->SetLogz();
   hKeTh->SetTitle("14C(p,p') recoil protons, E_{beam} = 159.75 MeV");
   hKeTh->Draw("colz");
   auto *leg = new TLegend(0.55, 0.62, 0.89, 0.89);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   std::vector<TGraph *> loci;
   for (int k = 0; k < 5; ++k) {
      auto *g = new TGraph();
      for (double a = 12; a <= 178; a += 0.5) {
         double v = ez_locus(m_C14, m_p, m_p, m_C14, Ebeam, a * TMath::DegToRad(), exLoci[k]);
         if (v > 0 && v < 40) g->SetPoint(g->GetN(), a, v);
      }
      g->SetLineColor(lc[k]);
      g->SetLineWidth(2);
      g->Draw("l same");
      leg->AddEntry(g, TString::Format("E_{x} = %.3f MeV", exLoci[k]), "l");
      loci.push_back(g);
   }
   leg->Draw();
   cK->cd(2);
   gPad->SetLogz();
   hCmOld->Draw("colz");
   cK->cd(3);
   gPad->SetLogz();
   hCmNew->Draw("colz");
   cK->SaveAs(plotDir + "ez_kinematics_" + tag + ".png");

   // (b) the excitation energy itself
   auto *cE = new TCanvas("cE", "excitation", 1500, 900);
   cE->Divide(2, 2);
   cE->cd(1);
   hOld->SetLineColor(kGray + 2);
   hOld->SetLineWidth(2);
   hNew->SetLineColor(kAzure + 2);
   hNew->SetLineWidth(2);
   hOld->SetTitle("14C(p,p') excitation energy");
   hOld->GetYaxis()->SetRangeUser(0, 1.15 * std::max(hOld->GetMaximum(), hNew->GetMaximum()));
   hOld->Draw("hist");
   hNew->Draw("hist same");
   auto *lg2 = new TLegend(0.52, 0.68, 0.89, 0.88);
   lg2->SetBorderSize(0);
   lg2->SetFillStyle(0);
   lg2->AddEntry(hOld, "constant E_{beam} (as analysed)", "l");
   lg2->AddEntry(hNew, "E_{beam}(z_{vertex})", "l");
   lg2->Draw();
   cE->cd(2);
   auto *hOldZ = hOldI;
   auto *hNewZ = hNewI;
   hOldZ->SetTitle(TString::Format("the 6-8 MeV group, #theta_{cm} > %.0f#circ", cmSens));
   hOldZ->SetLineColor(kGray + 2);
   hOldZ->SetLineWidth(2);
   hNewZ->SetLineColor(kAzure + 2);
   hNewZ->SetLineWidth(2);
   double ym = 0;
   for (int b = 1; b <= hOldZ->GetNbinsX(); ++b)
      ym = std::max(ym, std::max(hOldZ->GetBinContent(b), hNewZ->GetBinContent(b)));
   hOldZ->GetYaxis()->SetRangeUser(0, 1.25 * ym);
   hOldZ->Draw("hist");
   hNewZ->Draw("hist same");
   auto *lgI = new TLegend(0.52, 0.72, 0.89, 0.88);
   lgI->SetBorderSize(0);
   lgI->SetFillStyle(0);
   lgI->AddEntry(hOldZ, "constant E_{beam}", "l");
   lgI->AddEntry(hNewZ, "E_{beam}(z_{vertex})", "l");
   lgI->Draw();
   for (int k = 1; k < 5; ++k) {
      auto *l = new TLine(exLoci[k], 0, exLoci[k], 1.2 * ym);
      l->SetLineStyle(2);
      l->SetLineColor(kGray + 1);
      l->Draw();
   }
   cE->cd(3);
   gPad->SetLogz();
   hVzOld->Draw("colz");
   cE->cd(4);
   gPad->SetLogz();
   hVzNew->Draw("colz");
   cE->SaveAs(plotDir + "ez_excitation_" + tag + ".png");

   // (c) WHY the integrated spectrum barely moves: the sensitivity and the yield sit in different
   //     places. This panel decides whether the correction was worth applying, and it also shows
   //     that the correction is NOT the source of the known Ex-vs-angle drift -- it is roughly
   //     thirty times too small for that.
   {
      const int nS = 14;
      auto *gMedOld = new TGraphErrors();
      auto *gMedNew = new TGraphErrors();
      auto *gSigOld = new TGraph();
      auto *gSigNew = new TGraph();
      auto *gSens = new TGraph();
      auto *hYield = new TH1F("hYield", ";#theta_{cm} [deg];elastic yield", nS, 0, 140);
      std::vector<double> vOld[nS], vNew[nS], vSen[nS];
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (vz < 0 || vz > driftLen || thcm < 0 || thcm >= 140)
            continue;
         int b = (int)(thcm / (140.0 / nS));
         if (b < 0 || b >= nS)
            continue;
         double eb = Ebeam_at(vz, useSign);
         auto [exN, cmN] = ez_kine(m_C14, m_p, m_p, m_C14, eb, th * TMath::DegToRad(), ke);
         double a = std::get<0>(ez_kine(m_C14, m_p, m_p, m_C14, Ebeam + 0.5, th * TMath::DegToRad(), ke));
         double c = std::get<0>(ez_kine(m_C14, m_p, m_p, m_C14, Ebeam - 0.5, th * TMath::DegToRad(), ke));
         if (std::fabs(ex) < 1.5) {
            vOld[b].push_back(ex);
            hYield->Fill(thcm);
         }
         if (!std::isnan(exN) && std::fabs(exN) < 1.5)
            vNew[b].push_back(exN);
         if (!std::isnan(a) && !std::isnan(c))
            vSen[b].push_back(a - c);
      }
      auto qf = [](std::vector<double> v, double q) {
         size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, q * (v.size() - 1)));
         std::nth_element(v.begin(), v.begin() + k, v.end());
         return v[k];
      };
      for (int b = 0; b < nS; ++b) {
         double xc = (b + 0.5) * (140.0 / nS);
         if (vSen[b].size() > 20)
            gSens->SetPoint(gSens->GetN(), xc, qf(vSen[b], .5));
         if (vOld[b].size() < 40 || vNew[b].size() < 40)
            continue;
         double so = (qf(vOld[b], .75) - qf(vOld[b], .25)) / 1.349;
         double sn = (qf(vNew[b], .75) - qf(vNew[b], .25)) / 1.349;
         int i0 = gMedOld->GetN();
         gMedOld->SetPoint(i0, xc, qf(vOld[b], .5));
         gMedOld->SetPointError(i0, 0, 1.253 * so / std::sqrt((double)vOld[b].size()));
         gMedNew->SetPoint(i0, xc + 1.5, qf(vNew[b], .5));
         gMedNew->SetPointError(i0, 0, 1.253 * sn / std::sqrt((double)vNew[b].size()));
         gSigOld->SetPoint(gSigOld->GetN(), xc, so);
         gSigNew->SetPoint(gSigNew->GetN(), xc, sn);
      }
      auto *cW = new TCanvas("cW", "where", 1500, 480);
      cW->Divide(3, 1);
      cW->cd(1);
      gMedOld->SetTitle("elastic peak position;#theta_{cm} [deg];median E_{x} [MeV]");
      gMedOld->SetMarkerStyle(20);
      gMedOld->SetMarkerColor(kGray + 2);
      gMedOld->SetLineColor(kGray + 2);
      gMedOld->GetYaxis()->SetRangeUser(-1.2, 0.6);
      gMedOld->Draw("ap");
      gMedNew->SetMarkerStyle(24);
      gMedNew->SetMarkerColor(kAzure + 2);
      gMedNew->SetLineColor(kAzure + 2);
      gMedNew->Draw("p same");
      auto *l0 = new TLine(0, 0, 140, 0);
      l0->SetLineStyle(2);
      l0->SetLineColor(kGray + 1);
      l0->Draw();
      auto *lg3 = new TLegend(0.45, 0.18, 0.89, 0.34);
      lg3->SetBorderSize(0);
      lg3->SetFillStyle(0);
      lg3->AddEntry(gMedOld, "constant E_{beam}", "p");
      lg3->AddEntry(gMedNew, "E_{beam}(z_{vertex})", "p");
      lg3->Draw();
      cW->cd(2);
      gSigOld->SetTitle("elastic peak width;#theta_{cm} [deg];IQR/1.349 [MeV]");
      gSigOld->SetMarkerStyle(20);
      gSigOld->SetMarkerColor(kGray + 2);
      gSigOld->GetYaxis()->SetRangeUser(0, 1.0);
      gSigOld->Draw("apl");
      gSigNew->SetMarkerStyle(24);
      gSigNew->SetMarkerColor(kAzure + 2);
      gSigNew->SetLineColor(kAzure + 2);
      gSigNew->Draw("pl same");
      cW->cd(3);
      gPad->SetRightMargin(0.13);
      gSens->SetTitle("sensitivity vs where the yield is;#theta_{cm} [deg];dE_{x}/dE_{beam}");
      gSens->SetMarkerStyle(21);
      gSens->SetMarkerColor(kRed + 1);
      gSens->SetLineColor(kRed + 1);
      gSens->SetLineWidth(2);
      gSens->GetYaxis()->SetRangeUser(0, 0.13);
      gSens->Draw("apl");
      double smax = 0;
      for (int b = 1; b <= hYield->GetNbinsX(); ++b)
         smax = std::max(smax, hYield->GetBinContent(b));
      hYield->Scale(0.12 / (smax > 0 ? smax : 1));
      hYield->SetFillColorAlpha(kAzure + 1, 0.30);
      hYield->SetLineColor(kAzure + 2);
      hYield->Draw("hist same");
      gSens->Draw("pl same");
      auto *lg4 = new TLegend(0.40, 0.70, 0.86, 0.88);
      lg4->SetBorderSize(0);
      lg4->SetFillStyle(0);
      lg4->AddEntry(gSens, "dE_{x}/dE_{beam}", "pl");
      lg4->AddEntry(hYield, "elastic yield (scaled)", "f");
      lg4->Draw();
      cW->SaveAs(plotDir + "ez_where_" + tag + ".png");
      printf("wrote %sez_where_%s.png\n", plotDir.Data(), tag.Data());
   }

   fo.cd();
   out->Write();
   hOld->Write(); hNew->Write(); hOldI->Write(); hNewI->Write(); hVzOld->Write(); hVzNew->Write(); hCmOld->Write(); hCmNew->Write();
   TString recTitle = TString::Format("CATIMA E(z), rho=%.4e g/cm3, anchor E(%.1f mm)=%.2f MeV, sign %+d", density,
                                      zRef, Ebeam, useSign);
   TNamed rec(TString("ebeamz"), recTitle);
   rec.Write();
   fo.Close();
   printf("\nwrote %sproton_kin_%s_ez.root, %sez_kinematics_%s.png, %sez_excitation_%s.png\n",
          plotDir.Data(), tag.Data(), plotDir.Data(), tag.Data(), plotDir.Data(), tag.Data());
   printf("ebeamz data done\n\n");
}
