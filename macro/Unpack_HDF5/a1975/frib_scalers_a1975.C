/// @file frib_scalers_a1975.C
/// @brief Read the FRIBDAQ SCALERS from the raw a1975 HDF5, and turn them into a LUMINOSITY.
///
/// WHY THIS EXISTS. Both a1975 transfer channels are normalised by a luminosity, and until now
/// each rested on something indirect: (p,d) on the 16C(p,p) elastic yield divided by a DWBA
/// (so it inherits the optical-potential choice), and (d,t) on a beam count taken from the
/// 16C(d,p) analysis with no independent check at all. The scalers give a THIRD route that shares
/// nothing with either: a hardware count of beam particles through the ion chamber, times the gas
/// areal density. No elastic, no optical model, no acceptance.
///
/// ATTPCROOT DOES NOT READ THESE. AtFRIBHDFUnpacker walks `frib/evt` and stops; AtMergerHDFUnpacker
/// documents `/scalers` as "(not used here)". The counters were simply never needed. Spyral does
/// read them (spyral/trace/frib_scalers.py -> {run}_scaler.parquet), and this macro is that reader
/// adapted to the framework -- VALIDATED BIT-IDENTICAL against Spyral's parquet on run 0016, all
/// eleven counters.
///
/// WHERE THEY LIVE, in the legacy remerged layout a1975 uses:
///     frib/scaler/scaler<N>_data      Dataset {32} uint32   <- the counters
///     frib/scaler/scaler<N>_header    Dataset {5}  uint32
/// Only the first eleven of the 32 words are named by FRIBDAQ; the rest are unassigned.
///
/// *** THE SCALERS ARE INCREMENTAL. *** Each scaler event holds the counts SINCE THE LAST ONE, so
/// the run total is the SUM over all scaler events -- not the last value, which is the mistake this
/// interface invites. A run has a few thousand of them.
///
///     index  name           meaning
///       0    clock_free     time elapsed while the DAQ was running
///       1    clock_live     time for which the DAQ could accept triggers
///       2    trigger_free   triggers received
///       3    trigger_live   triggers that actually made events
///       4    ic_sca         ION CHAMBER counts -- the beam counter
///       5    mesh_sca       mesh signals
///       6    si1_cfd        Si detector 1
///       7    si2_cfd        Si detector 2
///       8    sipm           unclear (FRIBDAQ's own docstring says so)
///       9    ic_ds          downscaled ion chamber
///      10    ic_cfd         unclear; equals ic_sca in every a1975 run checked
///
/// THE LUMINOSITY, and every assumption in it stated:
///     L [mb^-1] = ic_sca * f_beam * livetime * n_target / 1e27
///   f_beam    the fraction of IC counts that are the species of interest. MEASURED for the H2
///             block at 61.3% carbon (the rest is mostly 21O -- NOT 16O, a separator at fixed Brho
///             selects one A/Q and 16O is 25% off it). NOT measured for D2; passing the H2 value is
///             an assumption and the macro says so.
///   livetime  clock_live/clock_free. trigger_live/trigger_free differs by ~1.8% and which one the
///             original beam count used is not recorded, so both are printed.
///   n_target  from P, T and the drift length. T is assumed, not logged.
///
/// HDF5 is not on the interpreter's default include path, so the macro adds it itself rather than
/// making the caller export ROOT_INCLUDE_PATH -- one less thing to get wrong at 2 a.m.
///
///   root -b -q 'frib_scalers_a1975.C(0,0,"/mnt/f/a1975/h5/",2,300,293,97.1731,0.613,"d2runs.txt")'
///   root -b -q 'frib_scalers_a1975.C(106,189,"/mnt/f/a1975/h5/",1.0,300,293,97.1731,0.613)'
// R__LOAD_LIBRARY at FILE scope, not gSystem->Load inside the function: cling JITs the whole
// function before it runs, so a load in the body is already too late and every H5:: symbol comes
// back unresolved. AtUnpack links these at build time; an interpreted macro must ask up front.
R__ADD_INCLUDE_PATH(/home/yassid/fair_install/hdf5/include)
R__LOAD_LIBRARY(/home/yassid/fair_install/hdf5/lib/libhdf5.so)
R__LOAD_LIBRARY(/home/yassid/fair_install/hdf5/lib/libhdf5_cpp.so)
#include <H5Cpp.h>

/// *** A RUN RANGE IS NOT A RUN LIST. *** The h5 directory holds 73 runs numbered 16-103 while the
/// (d,t) analysis uses 47 of them; summing the range instead of the list inflates the luminosity by
/// 56%. Pass runListFile -- one run number per line, e.g. generated from the reconstructed files
/// the analysis actually consumed -- and the range arguments are ignored.
void frib_scalers_a1975(int runLo = 16, int runHi = 103, TString h5dir = "/mnt/f/a1975/h5/",
                        // target: 2 for D2 (two deuterons per molecule), 2 for H2, 1 for a monatomic gas
                        double atomsPerMolecule = 2.0, double pTorr = 300.0, double tempK = 293.0,
                        double lengthCm = 97.1731, double fBeam = 0.613, TString runListFile = "")
{
   // The scaler count is not stored anywhere, so the only way to find it is to scan until the
   // dataset is missing -- which HDF5 reports as a 20-line diagnostic per run. That failure is the
   // NORMAL termination, not an error, so silence the printer and keep the exceptions.
   H5::Exception::dontPrint();

   const char *NAME[11] = {"clock_free", "clock_live", "trigger_free", "trigger_live", "ic_sca",
                           "mesh_sca",   "si1_cfd",    "si2_cfd",      "sipm",         "ic_ds",
                           "ic_cfd"};
   // the run list wins over the range, always
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
             "  used a subset, pass runListFile instead or this over-counts the beam.\033[0m\n",
             runLo, runHi, runs.size());
   }

   double tot[11] = {0};
   int nRun = 0, nMiss = 0;
   long nScalerEv = 0;
   printf("\n  %-6s %10s %14s %14s %10s\n", "run", "scal ev", "ic_sca", "trigger_live", "live");
   for (int r : runs) {
      TString f = Form("%srun_%04d.h5", h5dir.Data(), r);
      if (gSystem->AccessPathName(f)) { ++nMiss; continue; }
      double sub[11] = {0};
      long nev = 0;
      try {
         H5::H5File file(f.Data(), H5F_ACC_RDONLY);
         H5::Group g = file.openGroup("frib/scaler");
         while (true) {
            uint32_t buf[32];
            try {
               H5::DataSet ds = g.openDataSet(Form("scaler%ld_data", nev));
               ds.read(buf, H5::PredType::NATIVE_UINT32);
            } catch (...) { break; }
            for (int i = 0; i < 11; ++i) sub[i] += buf[i];
            ++nev;
         }
      } catch (...) { ++nMiss; continue; }
      if (nev == 0) { ++nMiss; continue; }
      ++nRun; nScalerEv += nev;
      for (int i = 0; i < 11; ++i) tot[i] += sub[i];
      const double lv = sub[0] > 0 ? sub[1]/sub[0] : 0;
      printf("  %-6d %10ld %14.0f %14.0f %10.4f\n", r, nev, sub[4], sub[3], lv);
   }
   if (nRun == 0) { printf("\n  no runs with scaler data found in %s\n\n", h5dir.Data()); return; }

   printf("\n  ===== TOTALS over %d runs (%d missing or without scalers), %ld scaler events =====\n",
          nRun, nMiss, nScalerEv);
   for (int i = 0; i < 11; ++i) printf("    %-14s %20.0f\n", NAME[i], tot[i]);

   const double liveClock = tot[0] > 0 ? tot[1]/tot[0] : 1.0;
   const double liveTrig = tot[2] > 0 ? tot[3]/tot[2] : 1.0;
   printf("\n    livetime  clock_live/clock_free     %.4f\n", liveClock);
   printf("    livetime  trigger_live/trigger_free %.4f\n", liveTrig);

   // target areal density
   const double kB = 1.380649e-23;
   const double nMol = pTorr*133.322/(kB*tempK)/1e6;      // molecules per cm3
   const double nTgt = atomsPerMolecule*nMol*lengthCm;    // atoms per cm2
   printf("\n    target %.0f torr, %.0f K, %.4f cm, %.1f atoms/molecule\n",
          pTorr, tempK, lengthCm, atomsPerMolecule);
   printf("    areal density %.4e atoms/cm2\n", nTgt);

   printf("\n  ===== LUMINOSITY, beam fraction %.3f =====\n", fBeam);
   printf("    %-34s %12s\n", "livetime treatment", "L [mb^-1]");
   printf("    %-34s %12.1f\n", "clock_live/clock_free", tot[4]*fBeam*liveClock*nTgt/1e27);
   printf("    %-34s %12.1f\n", "trigger_live/trigger_free", tot[4]*fBeam*liveTrig*nTgt/1e27);
   printf("    %-34s %12.1f\n", "none", tot[4]*fBeam*nTgt/1e27);
   if (nMiss > 0)
      printf("\n  \033[1;33m%d runs in [%d,%d] had no scaler data. If they carried beam, this is a\n"
             "  LOWER BOUND -- scale by the run count only if the exposure per run was equal.\033[0m\n",
             nMiss, runs.front(), runs.back());
   printf("\n  f_beam is an ASSUMPTION unless measured on this run block: 0.613 was measured on H2.\n\n");
}
