/// @file ic_spec_C15d.C
/// @brief Dump the ion-chamber spectrum that dd/inject_ic_panel.py injects into the explorer.
///
///   root -b -q 'dd/ic_spec_C15d.C()'
///   root -b -q 'dd/ic_spec_C15d.C("pid/points_C15d.root","dd/plots/ic_C15d.json")'
///
/// ★ WHY THIS EXISTS. The IC panel was injected BY HAND on 2026-09-01 and nothing in the pipeline
/// wrote its input, so every regeneration of the page since then silently dropped the panel --
/// the viewer kept only the numeric icLo/icHi boxes. `inject_ic_panel.py` names an "IC dump step"
/// in its docstring that did not exist in the tree. This is that step.
///
/// ★ THE SPECTRUM IS THE FULL COCKTAIL, THE FITTED SAMPLE IS NOT. This reads
/// pid/points_C15d.root, i.e. every reconstructed track's IC value (2.02 M), because seeing the
/// whole beam structure is the point of the panel. But the fits it filters were GATED BEFORE
/// FITTING on a window, so only that slice of the spectrum has tracks behind it. fitLo/fitHi
/// record the fitted range so the panel can grey out where dragging would select nothing --
/// without it the obvious experiment (drag onto the 2058 component) returns an empty page and
/// reads as a bug.
///
/// ★ npulse == 0 IS NOT MULTIPLICITY ZERO, IT IS NO PULSE FOUND, so its IC value means nothing and
/// it is counted with the no-usable-IC tracks rather than being dropped silently. On a2091 D2 that
/// is 38,511 tracks -- enough that losing them would move a percentage.
///
/// The JSON contract is set by the injector's JS, not by its docstring:
///   lo, hi, nb        histogram range and bin count
///   byNp[4][nb]       counts for multiplicity 1, 2, 3, 4+
///   cnt[4]            totals per multiplicity
///   nTot, nNoIC       all tracks, and those with no usable IC
///   gateLo, gateHi    the window the batch macros use, drawn dashed and used by the snap button
///   fitLo, fitHi      the range that actually has fitted tracks behind it (-1 = unrestricted)

void ic_spec_C15d(TString pointsFile = "pid/points_C15d.root", TString outJson = "dd/plots/ic_C15d.json",
                  Double_t lo = 0, Double_t hi = 3000, Int_t nb = 300,
                  // The a2091 CARBON peak, as gate_events_C15d.C applies it before fitting.
                  Double_t gateLo = 1045, Double_t gateHi = 1225,
                  // Range with fitted tracks behind it. Defaults to the gate because the fits were
                  // gated on it; pass -1,-1 once a pass has been fitted with the IC condition off.
                  Double_t fitLo = 1045, Double_t fitHi = 1225, Int_t runMin = 13, Int_t runMax = 133)
{
   if (gSystem->AccessPathName(pointsFile)) {
      printf("\033[1;31mERROR: %s not found (run pid/make_points_C15d.C)\033[0m\n", pointsFile.Data());
      return;
   }
   TFile *fp = TFile::Open(pointsFile);
   TTree *t = fp ? (TTree *)fp->Get("pts") : nullptr;
   if (!t) {
      printf("\033[1;31mERROR: no 'pts' tree in %s\033[0m\n", pointsFile.Data());
      return;
   }

   Int_t run = 0, np = 0;
   Float_t ic = -1;
   t->SetBranchAddress("run", &run);
   t->SetBranchAddress("npulse", &np);
   t->SetBranchAddress("ic", &ic);

   std::vector<std::vector<double>> byNp(4, std::vector<double>(nb, 0.0));
   std::vector<double> cnt(4, 0.0);
   Long64_t nTot = 0, nNoIC = 0, nOut = 0;

   const Long64_t n = t->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      if (run < runMin || run > runMax)
         continue;
      ++nTot;
      // A zero-pulse reading carries no amplitude; counting it as multiplicity would put a
      // meaningless ADC value into the histogram.
      if (!(ic >= 0) || np < 1) {
         ++nNoIC;
         continue;
      }
      const Int_t m = (np >= 4) ? 3 : np - 1;   // index 3 is "4 or more"
      cnt[m] += 1;
      const Int_t b = (Int_t)((ic - lo) / (hi - lo) * nb);
      if (b < 0 || b >= nb) {
         ++nOut;
         continue;
      }
      byNp[m][b] += 1;
   }
   fp->Close();

   gSystem->mkdir(gSystem->DirName(outJson), kTRUE);
   std::ofstream o(outJson.Data());
   o << "{\"lo\":" << lo << ",\"hi\":" << hi << ",\"nb\":" << nb << ",\"gateLo\":" << gateLo
     << ",\"gateHi\":" << gateHi << ",\"fitLo\":" << fitLo << ",\"fitHi\":" << fitHi
     << ",\"nTot\":" << nTot << ",\"nNoIC\":" << nNoIC << ",\"cnt\":[";
   for (int m = 0; m < 4; ++m)
      o << (m ? "," : "") << (Long64_t)cnt[m];
   o << "],\"byNp\":[";
   for (int m = 0; m < 4; ++m) {
      o << (m ? ",[" : "[");
      for (int b = 0; b < nb; ++b)
         o << (b ? "," : "") << (Long64_t)byNp[m][b];
      o << "]";
   }
   o << "]}\n";
   o.close();

   printf("\033[1;33m=== ion-chamber spectrum ===\033[0m\n");
   printf("  tracks    : %lld  (%lld with no usable IC = %.1f %%)\n", nTot, nNoIC, 100.0 * nNoIC / nTot);
   printf("  by np     : 1:%lld  2:%lld  3:%lld  4+:%lld\n", (Long64_t)cnt[0], (Long64_t)cnt[1],
          (Long64_t)cnt[2], (Long64_t)cnt[3]);
   if (nOut)
      printf("  \033[1;33moutside [%.0f, %.0f]: %lld\033[0m\n", lo, hi, nOut);
   printf("  gate      : [%.0f, %.0f]\n", gateLo, gateHi);
   if (fitLo >= 0)
      printf("  fitted    : [%.0f, %.0f]  -- outside this the panel greys out, no fits exist there\n",
             fitLo, fitHi);
   printf("  \033[1;32mwrote\033[0m %s\n", outJson.Data());
}
