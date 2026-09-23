/// @file exvz_dp.C
/// @brief Recompute the 16C(d,p) excitation energy with a VERTEX-DEPENDENT BEAM ENERGY, and show
/// the bound region before and after.
///
/// WHY. The analysis uses ONE beam energy for every vertex. The beam loses energy traversing the
/// gas, so a reaction 800 mm in happens at a lower energy than one 100 mm in, and the inversion
/// then returns an Ex that is too high -- by more at the far end. Measured on the g.s. simulation
/// with TRUTH ONLY (no reconstruction anywhere in it), generated Ex = 0 came back as:
///        <vz> mm      71    200    297    398    498    600    698    800    877
///        Ex (MeV)  +0.068 +0.134 +0.189 +0.243 +0.296 +0.356 +0.406 +0.461 +0.514
/// i.e. a **+0.45 MeV drift across the chamber that is not reconstruction at all**. The data's
/// bound peak walks +0.30 MeV over 150->850 where the simulation's truth walks +0.33
/// ([[project_a1975_dp_vz_drift]]).
///
/// THE CORRECTION. Scanning the beam energy that zeroes truth's Ex in each slice gives a line to
/// ~3 %:   Ebeam(vz) = 183.613 - 0.011030 * vz      [MeV, vz in mm]
/// i.e. **11.0 MeV/m**, ~8.9 MeV lost across the chamber.
///
/// *** ONLY THE SLOPE IS TAKEN FROM THE SIMULATION. THE ANCHOR IS THE DATA'S. ***
/// The fitted intercept is 183.613, not the 184.25 the simulation was GENERATED at -- the
/// difference is energy lost before the first slice plus whatever the generator and the analysis
/// inversion do not share. Importing that absolute offset into the data would silently re-calibrate
/// the energy scale on a simulation artefact. So the correction is applied as a PIVOT about a
/// reference vertex:
///        Ebeam(vz) = Ebeam0 - slope * (vz - vzRef)
/// which removes the DRIFT and leaves the overall scale exactly where the data's own calibration
/// put it. vzRef defaults to 500 mm, mid-chamber.
///
/// SIGN CHECK, and it matters: applied backwards this DOUBLES the drift instead of removing it.
/// The data's peak rises with vz (-0.190 -> +0.110) and the simulation's truth rises with vz
/// (+0.134 -> +0.461). Same direction, so the same sign applies. That is an empirical check, not
/// an assumption about which end the beam enters.
///
/// STILL NOT CORRECTED BY THIS: the sigma 0.927 resolution collapse at vz 650-750. The simulation
/// gives 0.564 there against 0.372-0.614 elsewhere -- no collapse -- so that one is in the data and
/// is not beam energy loss.
///
///   root -l 'exvz_dp.C()'

#include "exsel_dp.C"

namespace {
double vmed(std::vector<double> v)
{
   if (v.empty()) return -999;
   std::sort(v.begin(), v.end());
   const size_t n = v.size();
   return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
} // namespace

void exvz_dp(const char *fname = "/mnt/f/a1975/caches/.explorer/exp_dp_find805.root",
             double slope = 0.011030, double vzRef = 500.0, const char *png = "")
{
   using namespace exsel;
   Cuts s; // the adopted GOOD_ selection, keOff +0.25 and the theta correction included

   TFile *f = TFile::Open(fname);
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("ERROR: no tree pk in %s\n", fname); return; }
   float ke = 0, th = 0, vz = 0, c2 = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);

   const int NS = 8;
   const double zlo = 150, zw = 100;
   TH1D *hBefore = new TH1D("hBefore", "", 60, -3, 3);
   TH1D *hAfter = new TH1D("hAfter", "", 60, -3, 3);
   TH1D *hb[NS], *ha[NS];
   for (int i = 0; i < NS; ++i) {
      hb[i] = new TH1D(Form("hb%d", i), "", 60, -3, 3);
      ha[i] = new TH1D(Form("ha%d", i), "", 60, -3, 3);
   }
   long n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      const double k = ke + s.keOff, thc = thCorr(s, th, k);
      if (!(c2 <= s.chi2)) continue;
      if (!(thc >= s.thLo && thc <= s.thHi)) continue;
      if (!(k >= s.keLo && k <= s.keHi)) continue;
      if (!(vz >= s.vzLo && vz <= s.vzHi)) continue;
      double exB = 0, tcB = 0, exA = 0, tcA = 0;
      const double ebCorr = s.ebeam - slope * (vz - vzRef);
      if (!kine2b(s.ebeam, thc * TMath::Pi() / 180.0, k, exB, tcB)) continue;
      if (!kine2b(ebCorr, thc * TMath::Pi() / 180.0, k, exA, tcA)) continue;
      ++n;
      hBefore->Fill(exB);
      hAfter->Fill(exA);
      const int b = (int)((vz - zlo) / zw);
      if (b >= 0 && b < NS) { hb[b]->Fill(exB); ha[b]->Fill(exA); }
   }

   printf("\n=== exvz_dp : vertex-dependent beam energy on the DATA ===\n");
   printf("  Ebeam(vz) = %.2f - %.6f * (vz - %.0f)   [%.1f MeV/m]\n", s.ebeam, slope, vzRef, slope * 1000);
   printf("  at vz 150 -> %.2f MeV,  at vz 940 -> %.2f MeV   (span %.2f MeV)\n", s.ebeam - slope * (150 - vzRef),
          s.ebeam - slope * (940 - vzRef), slope * 790);
   printf("  tracks %ld\n\n", n);

   printf("  BOUND PEAK per vz slice, gaussian over -1.5..1.5\n");
   printf("    vz slice      n     BEFORE peak  sigma      AFTER peak  sigma\n");
   std::vector<double> pB, pA;
   for (int b = 0; b < NS; ++b) {
      const double nb = hb[b]->Integral();
      if (nb < 25) { printf("    %4.0f-%4.0f   %5.0f    (too few)\n", zlo + zw * b, zlo + zw * (b + 1), nb); continue; }
      TF1 g1("g1", "gaus", -1.5, 1.5), g2("g2", "gaus", -1.5, 1.5);
      g1.SetParameters(hb[b]->GetMaximum(), 0.0, 0.4);
      g2.SetParameters(ha[b]->GetMaximum(), 0.0, 0.4);
      hb[b]->Fit(&g1, "QNR");
      ha[b]->Fit(&g2, "QNR");
      printf("    %4.0f-%4.0f   %5.0f      %+6.3f  %6.3f      %+6.3f  %6.3f\n", zlo + zw * b, zlo + zw * (b + 1), nb,
             g1.GetParameter(1), TMath::Abs(g1.GetParameter(2)), g2.GetParameter(1), TMath::Abs(g2.GetParameter(2)));
      pB.push_back(g1.GetParameter(1));
      pA.push_back(g2.GetParameter(1));
   }
   if (pB.size() > 2) {
      auto spread = [](std::vector<double> &v) {
         return *std::max_element(v.begin(), v.end()) - *std::min_element(v.begin(), v.end());
      };
      printf("\n  PEAK SPREAD ACROSS THE CHAMBER   before %.3f MeV   ->   after %.3f MeV\n", spread(pB), spread(pA));
      printf("  (the 17C bound states span 0.331 MeV, g.s. to 5/2+ -- that is the scale to beat)\n");
   }

   TF1 gb("gb", "gaus", -1.5, 1.5), ga("ga", "gaus", -1.5, 1.5);
   gb.SetParameters(hBefore->GetMaximum(), 0, 0.4);
   ga.SetParameters(hAfter->GetMaximum(), 0, 0.4);
   hBefore->Fit(&gb, "QNR");
   hAfter->Fit(&ga, "QNR");
   printf("\n  WHOLE SAMPLE SUMMED\n");
   printf("    before   peak %+.3f   sigma %.3f\n", gb.GetParameter(1), TMath::Abs(gb.GetParameter(2)));
   printf("    after    peak %+.3f   sigma %.3f   <- narrower means the states were being smeared\n",
          ga.GetParameter(1), TMath::Abs(ga.GetParameter(2)));
   printf("    simulation resolution FLOOR: sigma(Ex) = 0.475 MeV\n\n");

   if (png && png[0]) {
      gStyle->SetOptStat(0);
      TCanvas *c = new TCanvas("c", "", 1200, 700);
      hBefore->SetTitle("16C(d,p)17C bound region: constant vs vertex-dependent E_{beam};E_{x} [MeV];counts");
      hBefore->SetLineColor(kGray + 2);
      hBefore->SetLineWidth(2);
      hBefore->Draw("hist");
      hAfter->SetLineColor(kRed + 1);
      hAfter->SetLineWidth(2);
      hAfter->Draw("hist same");
      TLegend *l = new TLegend(0.62, 0.72, 0.89, 0.88);
      l->AddEntry(hBefore, "E_{beam} = 184.25 (constant)", "l");
      l->AddEntry(hAfter, "E_{beam}(v_{z}) corrected", "l");
      l->Draw();
      c->SaveAs(png);
      printf("  wrote %s\n\n", png);
   }
}
