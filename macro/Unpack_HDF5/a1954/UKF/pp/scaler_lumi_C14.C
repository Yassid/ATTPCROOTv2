/// @file scaler_lumi_C14.C
/// @brief The a1954 luminosity from the FRIBDAQ scalers -- no optical model, no elastic, no acceptance.
///
/// The guide (13_absnorm.tex) opens "The measurement has no independent luminosity: there is no
/// beam counter reading to divide by", and the whole normalisation chain is built on that. It is
/// not true. Every remerged a1954 run carries /frib/scaler with 32 words per record, and
/// AtUnpack/AtFribScalers already reads it (validated bit-identical against Spyral on a1975).
/// No a1954 macro had ever referenced it.
///
///     L [mb^-1] = ic_sca * f_beam * livetime * n_target / 1e27
///
/// ic_sca      hardware ion-chamber count. THE COUNTERS ARE INCREMENTAL -- a run total is the sum
///             over every scaler record, not the last one. AtFribScalers handles that.
/// f_beam      fraction of IC counts that are 14C. MEASURE it with ic_fraction_C14.C; do not
///             inherit a1975's 0.613. It multiplies L directly and is the softest input.
/// livetime    clock and trigger differ by ~8% here; both are written out rather than one chosen.
/// n_target    from P, T and the SAME z window the yields were counted over (10-490 mm = 48 cm).
///
/// *** A RUN RANGE IS NOT A RUN LIST. *** The a1954 h5 directory holds 96 runs; the C14 analysis
/// uses 14 of them, and 0067 is absent from the middle of the block. Pass the list.
///
/// MEASURED, 14 runs, 2026-08-27 (plots/ is gitignored, so recorded here):
///     scaler records 40318, ic_sca 150422640, live fraction clock 0.9535 / trigger 0.8883
///     n_target 9.4917e20 /cm2 (300 torr, 293 K, 48 cm)
///     with f_beam 0.7112 ->  L = 96.8 (clock) / 90.2 (trigger) / 101.5 (none) mb^-1
///
/// WHAT THIS IS FOR, AND WHAT IT IS NOT FOR. It is NOT the normalisation. Cross sections must use
/// the ELASTIC luminosity, because that is a ratio measurement: any efficiency common to the
/// elastic and inelastic channels cancels between them, whereas the scalers count beam and know
/// nothing about detection, so normalising to them leaves every cross section low by that
/// efficiency. What the scalers measure that nothing else can is the efficiency itself,
/// eps = L_elastic / L_scaler: 0.59 (KD03), 0.74 (Perey), 0.80 (CH89, Menet), 0.86 (Becchetti-
/// Greenlees). eps <= 1 is a hard bound and all five respect it, but the acceptance already
/// carries the chi2 cut, so eps is trigger x PID-gate x residual DAQ loss -- 14-20% is ordinary,
/// 41% is a strain, which is a third independent reason to distrust KD03.
///
///   root -b -q 'scaler_lumi_C14.C("55,56,57,58,59,60,61,62,63,64,65,66,68,69")'
void scaler_lumi_C14(TString runsCSV = "55,56,57,58,59,60,61,62,63,64,65,66,68,69", double fBeam = -1,
                     TString h5dir = "/mnt/h/a1954_remerged/", double pTorr = 300.0, double tempK = 293.0,
                     double lengthCm = 48.0, double atomsPerMolecule = 2.0)
{
   gSystem->Load("libAtUnpack.so");
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   // f_beam: prefer the measured value on disk over anything passed in
   TString fbSrc = "argument";
   if (fBeam <= 0) {
      std::ifstream fi((here + "/plots/ic_fraction.txt").Data());
      if (fi) {   // skip ANY number of comment lines, not exactly one
         std::string line;
         while (std::getline(fi, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream is(line);
            double v;
            if (is >> v) { fBeam = v; fbSrc = "plots/ic_fraction.txt"; }
            break;
         }
      }
      if (fBeam <= 0) { printf("\033[1;31m  no f_beam: run ic_fraction_C14.C or pass one\033[0m\n"); return; }
   }

   std::vector<int> runs;
   { TObjArray *a = runsCSV.Tokenize(","); for (auto o : *a) runs.push_back(((TObjString *)o)->String().Atoi()); }
   AtFribScalers total; int nMiss = 0;
   for (int r : runs) {
      AtFribScalers one;
      if (!one.ReadRun(Form("%srun_%04d.h5", h5dir.Data(), r))) { ++nMiss; continue; }
      total.Add(one);
   }
   if (!total.GetRuns()) { printf("\033[1;31m  no scaler data under %s\033[0m\n", h5dir.Data()); return; }

   const double kB = 1.380649e-23;
   const double nMol = pTorr * 133.322 / (kB * tempK) / 1e6;
   const double nTgt = atomsPerMolecule * nMol * lengthCm;
   const double ic = total.Get(AtFribScalers::kIcSca);
   const double lc = total.GetLiveFractionClock(), lt = total.GetLiveFractionTrigger();
   const double Lc = ic * fBeam * lc * nTgt / 1e27, Lt = ic * fBeam * lt * nTgt / 1e27,
                Ln = ic * fBeam * nTgt / 1e27;

   printf("\n  ===== a1954 scaler luminosity =====\n");
   printf("    runs                %d of %zu (%d without scalers)\n", total.GetRuns(), runs.size(), nMiss);
   printf("    scaler records      %ld\n", total.GetScalerEvents());
   printf("    ic_sca              %lld\n", total.Get(AtFribScalers::kIcSca));
   printf("    live fraction       clock %.4f   trigger %.4f\n", lc, lt);
   printf("    f_beam              %.4f   (%s)\n", fBeam, fbSrc.Data());
   printf("    n_target            %.4e atoms/cm2  (%.0f torr, %.0f K, %.1f cm)\n", nTgt, pTorr, tempK, lengthCm);
   printf("\n    L clock             %8.1f mb^-1\n    L trigger           %8.1f mb^-1\n"
          "    L no livetime       %8.1f mb^-1\n", Lc, Lt, Ln);

   std::ofstream o((here + "/plots/scaler_luminosity.txt").Data());
   o << "# a1954 scaler luminosity, scaler_lumi_C14.C\n";
   o << "# runs=" << runsCSV << " h5dir=" << h5dir << " P=" << pTorr << " T=" << tempK
     << " len=" << lengthCm << " fBeam=" << fBeam << " (" << fbSrc << ")\n";
   o << "# ic_sca liveClock liveTrigger nTarget L_clock L_trigger L_nolive\n";
   o << (long long)ic << " " << lc << " " << lt << " " << nTgt << " " << Lc << " " << Lt << " " << Ln << "\n";
   printf("\n  wrote plots/scaler_luminosity.txt\n\n");
}
