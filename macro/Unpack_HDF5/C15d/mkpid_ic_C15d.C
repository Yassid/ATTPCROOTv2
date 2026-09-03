/// @file mkpid_ic_C15d.C
/// @brief PID plane built from the JOINED points file, so an IC beam window can be applied.
///
///   root -b -q 'mkpid_ic_C15d.C()'                        // carbon window, single pulse
///   root -b -q 'mkpid_ic_C15d.C(-1,-1)'                   // no IC cut, for comparison
///
/// WHY NOT mkpid_C15d.C. That one reads the per-run <run>_pid.root ntuples, which carry no ion
/// chamber -- the IC is joined later, in pid/make_points_C15d.C. So a plane from mkpid is always
/// the full beam cocktail. On these a2091 D2 runs the cocktail is Li/B/C/N/O with CARBON ONLY
/// ~40 %, so that plane mixes five beams.
///
/// This macro reads pid/points_C15d.root, which has run, sqrtdedx, brho, ic and npulse together,
/// and is the SAME file pid/gate_draw_C15d.C reads. That matters: the plane a gate is drawn on and
/// the plane it is judged on should be the same plane, and they have drifted apart before.
///
/// ⚠ SINGLE PULSE RIDES WITH THE WINDOW. apply_gate_C15d.C enforces npulse == 1 whenever an IC
/// window is active, so this does too. Applying the window without it gives a plane ~8 % larger
/// than the one the gate will actually be applied to.

#include "gain_C15d.h"

void mkpid_ic_C15d(Double_t icLo = 1045, Double_t icHi = 1225, TString pointsFile = "pid/points_C15d.root",
                   TString outDir = "plots/", TString outName = "pid_ic_C15d",
                   Int_t nbx = 340, Double_t xhi = 85.0, Int_t nby = 300, Double_t yhi = 2.5,
                   Bool_t requireSinglePulse = kTRUE, Int_t runMin = 13, Int_t runMax = 133)
{
   gSystem->mkdir(outDir, kTRUE);
   gStyle->SetOptStat(0);

   if (gSystem->AccessPathName(pointsFile)) {
      std::cout << "\033[1;31mERROR: " << pointsFile << " not found (run pid/make_points_C15d.C)\033[0m\n";
      return;
   }
   TFile *f = TFile::Open(pointsFile);
   TTree *t = f ? (TTree *)f->Get("pts") : nullptr;
   if (!t) {
      std::cout << "\033[1;31mERROR: no tree 'pts' in " << pointsFile << "\033[0m\n";
      return;
   }

   Int_t run = 0, npulse = 0;
   Float_t sq = 0, br = 0, ic = -1;
   t->SetBranchAddress("run", &run);
   t->SetBranchAddress("npulse", &npulse);
   t->SetBranchAddress("sqrtdedx", &sq);
   t->SetBranchAddress("brho", &br);
   t->SetBranchAddress("ic", &ic);

   const bool icOn = (icLo >= 0 && icHi > icLo);
   auto *h = new TH2F("hpid",
                      Form("C15d PID, gain matched, %s;#sqrt{dE/dx};B#rho [T#upointm]",
                           icOn ? Form("IC [%.0f, %.0f], single pulse", icLo, icHi) : "IC gate OFF"),
                      nbx, 0, xhi, nby, 0, yhi);

   Long64_t nAll = 0, nRange = 0, nNoIC = 0, nIn = 0;
   std::set<int> runs;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      ++nAll;
      if (run < runMin || run > runMax) continue;
      ++nRange;
      if (icOn) {
         // a track whose run has no IC carries ic = -1 and CANNOT satisfy a window. Dropping it is
         // right, but it has to be counted rather than vanish -- that is how a missing IC pass
         // turns into a mysteriously small sample.
         if (ic < 0) { ++nNoIC; continue; }
         if (ic < icLo || ic > icHi) continue;
         if (requireSinglePulse && npulse != 1) continue;
      }
      ++nIn;
      runs.insert(run);
      h->Fill(sq, br);
   }

   auto *c = new TCanvas("cpid", "pid", 1100, 850);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");
   c->SaveAs(outDir + outName + ".png");
   TFile fo(outDir + outName + ".root", "RECREATE");
   h->Write();
   fo.Close();
   f->Close();

   printf("\033[1;33m=== mkpid_ic_C15d ===\033[0m\n");
   printf("  points     : %s\n", pointsFile.Data());
   printf("  IC window  : %s%s\n", icOn ? Form("[%.0f, %.0f]", icLo, icHi) : "OFF",
          (icOn && requireSinglePulse) ? " with single pulse" : "");
   printf("  tracks     : %lld total, %lld in runs %d-%d\n", nAll, nRange, runMin, runMax);
   if (icOn) printf("  no usable IC: %lld dropped\n", nNoIC);
   printf("  \033[1;32mIN PLANE   : %lld  (%.1f %% of the run range), %zu runs\033[0m\n", nIn,
          nRange ? 100.0 * nIn / nRange : 0.0, runs.size());
   printf("  wrote %s%s.{png,root}\n", outDir.Data(), outName.Data());
}
