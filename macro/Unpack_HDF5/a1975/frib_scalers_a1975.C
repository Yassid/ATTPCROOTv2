/// @file frib_scalers_a1975.C
/// @brief Turn the FRIBDAQ scalers into a LUMINOSITY. Reading is done by AtFribScalers (AtUnpack).
///
/// WHY THIS EXISTS. Both a1975 transfer channels are normalised by a luminosity, and each rested
/// on something indirect: (p,d) on the 16C(p,p) elastic yield divided by a DWBA, so it inherits
/// the optical-potential choice; (d,t) on a beam count taken from the 16C(d,p) analysis with no
/// independent check at all. The scalers give a third route that shares nothing with either -- a
/// hardware count of beam through the ion chamber, times the gas areal density. No elastic, no
/// optical model, no acceptance.
///
/// The reading lives in AtFribScalers, not here: it belongs to the framework, it is validated
/// bit-identical against Spyral on run 0016, and it handles both HDF5 layouts. See that header for
/// the counter indices and for the warning that THE COUNTERS ARE INCREMENTAL.
///
/// *** A RUN RANGE IS NOT A RUN LIST. *** The a1975 h5 directory holds 73 runs numbered 16-103
/// while the (d,t) analysis uses 47 of them; summing the range instead of the list gave 483.7
/// against 330.0 and looked perfectly plausible. Pass runListFile -- one run number per line,
/// generated from the reconstructed files the analysis actually consumed -- and the range is
/// ignored. Handed a bare range this macro says so in yellow.
///
/// THE LUMINOSITY, with every assumption named:
///     L [mb^-1] = ic_sca * f_beam * livetime * n_target / 1e27
///   f_beam    fraction of IC counts that are the species of interest. MEASURED per EVENT on the
///             H2 block at 0.613 carbon; measured per TRACK on D2 at 0.659 (see ic_spectrum_d2.C
///             in the analysis repo). The rest is mostly 21O -- NOT 16O: a separator at fixed Brho
///             selects one A/Q and 16O is 25% off it, which is why Be is absent from both blocks.
///             This multiplies L directly and is its softest input.
///   livetime  clock or trigger; they differ by ~2% and which belongs here depends on what the
///             yield was normalised to, so both are printed rather than one being chosen for you.
///   n_target  from P, T and the drift length. T is assumed, not logged.
///
///   root -b -q 'frib_scalers_a1975.C(0,0,"/mnt/f/a1975/h5/",2,300,293,97.1731,0.613,"d2_runs.txt")'
///   root -b -q 'frib_scalers_a1975.C(106,189)'      // H2 block, by range -- warns
void frib_scalers_a1975(int runLo = 16, int runHi = 103, TString h5dir = "/mnt/f/a1975/h5/",
                        // 2 for D2 and for H2 (two atoms per molecule), 1 for a monatomic gas
                        double atomsPerMolecule = 2.0, double pTorr = 300.0, double tempK = 293.0,
                        double lengthCm = 97.1731, double fBeam = 0.613, TString runListFile = "")
{
   gSystem->Load("libAtUnpack.so");

   std::vector<int> runs;
   if (runListFile != "") {
      std::ifstream in(runListFile.Data());
      if (!in) { printf("\n  cannot open run list %s\n\n", runListFile.Data()); return; }
      int r;
      while (in >> r) runs.push_back(r);
      printf("\n  run LIST %s: %zu runs\n", runListFile.Data(), runs.size());
   } else {
      for (int r = runLo; r <= runHi; ++r) runs.push_back(r);
      printf("\n  \033[1;33mrun RANGE %d-%d (%zu numbers). A range is not a run list -- if the analysis\n"
             "  used a subset, pass runListFile or this over-counts the beam.\033[0m\n",
             runLo, runHi, runs.size());
   }

   AtFribScalers total;
   int nMiss = 0;
   printf("\n  %-6s %10s %14s %14s %10s\n", "run", "scal ev", "ic_sca", "trigger_live", "live");
   for (int r : runs) {
      TString f = Form("%srun_%04d.h5", h5dir.Data(), r);
      AtFribScalers one;
      if (!one.ReadRun(f.Data())) { ++nMiss; continue; }
      printf("  %-6d %10ld %14lld %14lld %10.4f\n", r, one.GetScalerEvents(),
             one.Get(AtFribScalers::kIcSca), one.Get(AtFribScalers::kTriggerLive),
             one.GetLiveFractionClock());
      total.Add(one);
   }
   if (total.GetRuns() == 0) { printf("\n  no runs with scaler data under %s\n\n", h5dir.Data()); return; }

   printf("\n  ===== TOTALS (%d missing or without scalers) =====\n", nMiss);
   total.Print();

   const double kB = 1.380649e-23;
   const double nMol = pTorr * 133.322 / (kB * tempK) / 1e6; // molecules per cm3
   const double nTgt = atomsPerMolecule * nMol * lengthCm;   // atoms per cm2
   printf("\n    target %.0f torr, %.0f K, %.4f cm, %.1f atoms/molecule -> %.4e atoms/cm2\n",
          pTorr, tempK, lengthCm, atomsPerMolecule, nTgt);

   const double ic = total.Get(AtFribScalers::kIcSca);
   printf("\n  ===== LUMINOSITY, beam fraction %.3f =====\n", fBeam);
   printf("    %-34s %12s\n", "livetime treatment", "L [mb^-1]");
   printf("    %-34s %12.1f\n", "clock_live/clock_free", ic * fBeam * total.GetLiveFractionClock() * nTgt / 1e27);
   printf("    %-34s %12.1f\n", "trigger_live/trigger_free", ic * fBeam * total.GetLiveFractionTrigger() * nTgt / 1e27);
   printf("    %-34s %12.1f\n", "none", ic * fBeam * nTgt / 1e27);
   if (nMiss > 0)
      printf("\n  \033[1;33m%d of %zu runs had no scaler data. If they carried beam this is a LOWER BOUND;\n"
             "  scale by the run count only if the exposure per run was equal.\033[0m\n", nMiss, runs.size());
   printf("\n  f_beam is an ASSUMPTION unless measured on THIS run block.\n\n");
}
