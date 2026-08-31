/// @file fitGenfitter_Ar46.C
/// @brief Fit the 46Ar(3He,d)47K simulation with EventFit::AtGenfitter, DEUTERON hypothesis.
///
/// Reads <tag>_reco.root (AtPatternEvent) -- the files the accumulation already wrote. Nothing is
/// re-digitised or re-reconstructed. Writes <tag>_genfitter_d<suffix>.root with AtTrackingEvent +
/// AtPIDEvent.
///
/// ============================================================================================
/// REWRITTEN 2026-08-30. EVERY FIT MADE WITH THE PREVIOUS VERSION IS VOID IN THE BACKWARD REGION.
/// ============================================================================================
///
/// This channel is a BACKWARD-EJECTILE channel: theta_lab(d) runs 131.4 deg at theta_cm 15 down to
/// 59.8 deg at theta_cm 80, so more than half the proposal's own angular range is past 90 deg. It
/// is therefore governed by the defects found on 14C(d,p) in August 2026, not by the (p,p')
/// lineage these macros were copied from:
///
///   * AtGenfitter added its measurements to the genfit track in ascending z_lab WHATEVER the
///     direction of travel, while the seed sits at the VERTEX -- which for a backward track is the
///     HIGHEST z_lab. TrackPoint 0 was therefore the STOPPING end: the momentum read out was the
///     one the ejectile had after shedding most of its energy, and the back-extrapolation ran the
///     energy loss the wrong way along the track. Fixed in 54a3726c (2026-08-29) by adding the
///     measurements in vertex-first order when backwardSeed is set. On 14C(d,p) at 2.85 T that
///     took backward KE_reco/KE_true from 0.834 to 0.999 and the spread from 19.1 % to 0.7 %.
///     Forward tracks are byte-for-byte unchanged.
///   * measSigma 4.0 mm (inherited from (p,p')) sits ~7x above the real hit residuals on the
///     AT-TPC pad plane, which chi2/ndf implies are 0.59-0.64 mm. At 4.0 a good fit lands at
///     chi2/ndf ~0.02, so the production cut at 5 selects on TRACK LENGTH, not fit quality, and
///     throws away 29 % of the longest spirals for nothing. DEFAULT IS NOW 0.6.
///   * matEffects + CATIMA dE/dx are the adopted setting. genfit applies ZERO stopping power below
///     beta*gamma = 0.05 unless a dE/dx source is loaded, and a 5 MeV deuteron is below that.
///
/// PAD PLANE: the real AT-TPC plane WITH the beam hole, deliberately. The 2 mm study is parked --
/// the hole has to stay so a telescope can be put behind it for the heavy fragment.
///
/// B FIELD SIGN IS NEGATIVE, MEASURED 2026-08-30, AND THE OLD HEADER HERE SAID THE OPPOSITE.
/// The previous version asserted +2.85 "because the simulation reverses the drift-z handedness in
/// digitisation". That handedness claim is true and unrelated: it concerns the UKF's bFieldSign
/// argument and the vertex-z mirror (r = -1.000, z_true + z_reco = 101.2 cm), not the SIGNED field
/// handed to AtGenfitter. The sign was never testable before, because it acts only through the
/// backward-seed path and that path was broken until 54a3726c.
///
/// A/B on gs_s3001, 8000 entries, identical in every other argument, against MC truth:
///
///     theta_lab      +2.85 T  KE_reco/KE_true (IQR)   |   -2.85 T  KE_reco/KE_true (IQR)
///      50- 90          1.003-1.009  (0.09-0.14)        |    1.003-1.009  (0.09-0.14)   identical
///      90-100          0.997  (0.071)                  |    0.997  (0.074)
///     100-110          0.993  (0.057)                  |    0.995  (0.051)
///     110-125          0.010  (0.005)   <-- COLLAPSE   |    0.991  (0.031)
///     125-145          0.974  (0.978)   <-- bimodal    |    0.984  (0.022)
///     fits/truth        94.1 %                         |     98.2 %
///     wall time         104 s                          |      16 s   (failures burn iterations)
///
/// The collapse is coherent, not noise: 194 tracks at 1 % of the true energy with an IQR of 0.005,
/// and a median dTheta of -0.86 deg. It sits exactly where the proposal's forward theta_cm
/// (15-30 deg) lands, so a fit made at the wrong sign destroys the most important part of the
/// range while leaving the forward arm untouched -- it does not present as an error.
///
/// NO PID GATE BY DEFAULT, ON PURPOSE. AtGenfitter owns its own AtSpyralPID and computes the gate
/// observables at the class default fMinPoints = 30 -- there is no SetMinPoints passthrough. A
/// gate drawn on the mp15 plane would be applied to mp30 observables inside the fitter. Since the
/// PID is fit-independent, the honest order is: fit everything, gate afterwards.
///
///   root -b -q 'fitGenfitter_Ar46.C("gs_s3001", 500, "/mnt/f/ar46_3hed/", "_probe")'   // quick
///   root -b -q 'fitGenfitter_Ar46.C("gs_s3001", -1, "/mnt/f/ar46_3hed/")'              // full

void fitGenfitter_Ar46(TString fileName = "gs_s3001", Long64_t nEvents = -1,
                       TString ioDir = "/mnt/f/ar46_3hed/", TString outSuffix = "", TString outDir = "",
                       Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5, TString pidGate = "",
                       // 0.6 mm, not 4.0: see the header. This is the AT-TPC-plane value; a 2 mm
                       // plane would want ~0.35, but that configuration is parked.
                       Double_t measSigma = 0.6, Double_t thetaMinDeg = 10.0, Double_t thetaMaxDeg = 170.0,
                       Bool_t matEffects = kTRUE, Bool_t backwardSeedFix = kTRUE,
                       TString geoName = "ATTPC_He3CO2_300torr", TString parName = "ATTPC.46Ar_3Hed_sim.par",
                       // matFallback OFF gives a clean matFX-only sample: tracks whose
                       // material-effects fit threw are DROPPED rather than silently refitted
                       // without material effects and mixed into the same output.
                       Bool_t matFallback = kFALSE, Bool_t backExtrap = kTRUE,
                       // Manual eloss is the matFX-OFF substitute for material effects, NOT a
                       // companion to them: with matEffects ON, extrapolateToLine already
                       // integrates the vertex gap and the two together count it twice.
                       Double_t manualElossDensity = 0.0, Int_t matA = 3,
                       Bool_t catimaELoss = kTRUE, Bool_t catimaELossFull = kFALSE,
                       Bool_t catimaMSC = kFALSE, Bool_t catimaStraggling = kFALSE,
                       // DO NOT USE THE RANGE CONSTRAINT without re-validating it here. On
                       // 14C(d,p) its gates reached only 4 % of backward spirals while DESTROYING
                       // 9 % of forward tracks (KE_fit/KE_true 0.993 -> 0.017), and because those
                       // then fail a truth match the damage presents as a null result.
                       Bool_t rangeConstraint = kFALSE, Double_t rangeDensity = 0.0, Int_t rangeMatA = 3,
                       Bool_t seedFromSpyral = kFALSE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   // DEUTERON ejectile for 46Ar(3He,d)47K
   const Int_t pdg = 1000010020;
   const Double_t massAmu = 2.01355321275;
   const Int_t Z = 1;

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/" + geoName + "_geomanager.root";
   TString digiParFile = dir + "/parameters/" + parName;
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter_d" + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(geoManFile.Data())) {
      std::cout << "\033[1;31mERROR: " << geoManFile << " not found; run geometry/" << geoName << ".C.\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(digiParFile.Data())) {
      std::cout << "\033[1;31mERROR: par " << digiParFile << " not found.\033[0m\n";
      return;
   }

   // FULL PARAMETER DUMP. A summary line once hid a fit parameter sitting on its bound; the same
   // applies to a configuration. Everything that can change the answer is printed, so a log file
   // is a complete record of the run that produced the file next to it.
   std::cout << "\033[1;33m=== fitGenfitter_Ar46 (46Ar(3He,d)47K, DEUTERON hyp) ===\033[0m\n"
             << "  in        : " << inputFile << "\n"
             << "  out       : " << outputFile << "\n"
             << "  geometry  : " << geoManFile << "\n"
             << "  par       : " << digiParFile << "\n"
             << "  particle  : pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << " (deuteron)\n"
             << "  gas       : 3He + 5% CO2, 300 torr, 8.3128e-5 g/cm3\n"
             << "  bField    : " << bField << " T   (NEGATIVE is correct here -- measured, see header)\n"
             << "  iter      : " << minIter << "-" << maxIter << "\n"
             << "  measSigma : " << measSigma << " mm\n"
             << "  theta win : [" << thetaMinDeg << ", " << thetaMaxDeg << "] deg\n"
             << "  matEffects: " << (matEffects ? "ON" : "OFF") << "   fallback "
             << (matFallback ? "on (retries flagged)" : "OFF (failures dropped)") << "\n"
             << "  backSeed  : " << (backwardSeedFix ? "ON" : "OFF") << "\n"
             << "  backExtrap: " << (backExtrap ? "ON" : "OFF") << "\n"
             << "  nEvents   : " << nEvents << "\n";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   // genfit material/propagation needs gGeoManager loaded before the fit
   if (gROOT->FindObject("FAIRGeom") == nullptr) {
      TFile *gf = TFile::Open(geoManFile);
      gf->Get("FAIRGeom");
   }

   // empty eloss file -> genfit internal; noMatEffects = !matEffects
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, /*elossFile*/ "", !matEffects,
                                                         minIter, maxIter);
   auto *fitterObs = fitter.get(); // non-owning: `fitter` is moved into the task below

   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);

   fitter->SetBackwardSeedFix(backwardSeedFix);
   if (backwardSeedFix)
      std::cout << "  \033[1;36mBackward seed-fix: ON\033[0m  (over half these deuterons are past 90 deg)\n";

   // The gas between the reaction vertex and the first cluster is unmeasured, so its energy loss
   // is missing from |p| at TrackPoint 0. Result lands in GetKinematicsXtr(), NOT GetKinematics().
   fitter->SetBackExtrapToAxis(backExtrap);
   if (backExtrap)
      std::cout << "  Back-extrapolation to beam axis: ON (read GetKinematicsXtr)\n";

   if (manualElossDensity > 0 && matEffects) {
      std::cout << "\033[1;31mREFUSING manual eloss: matEffects is ON, so extrapolateToLine already "
                   "integrates the vertex gap and this would count it twice. Ignored.\033[0m\n";
   } else if (manualElossDensity > 0) {
      fitter->SetManualELoss(manualElossDensity, matA);
      std::cout << "  Manual CATIMA eloss over the vertex gap: rho=" << manualElossDensity
                << " g/cm3, target A=" << matA << "\n";
   }

   if (catimaMSC || catimaStraggling) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catima material flags set but matEffects is OFF -- inert.\033[0m\n";
      fitter->SetCatimaMaterial(catimaMSC, catimaStraggling);
      std::cout << "  \033[1;35mCATIMA material model: MSC " << (catimaMSC ? "ON" : "off") << ", straggling "
                << (catimaStraggling ? "ON" : "off") << "\033[0m\n";
   }
   fitter->SetCatimaELoss(catimaELoss, catimaELossFull);
   if (catimaELoss) {
      if (!matEffects)
         std::cout << "\033[1;31mWARNING: catimaELoss set but matEffects is OFF -- inert.\033[0m\n";
      else
         std::cout << "  \033[1;32mdE/dx from CATIMA"
                   << (catimaELossFull ? " over the FULL range (Bethe-Bloch replaced too)"
                                       : " below beta*gamma=0.05 (Bethe-Bloch kept above)")
                   << ", per-step material\033[0m\n";
   } else if (matEffects) {
      std::cout << "\033[1;31mWARNING: matEffects ON with NO dE/dx source -- genfit applies zero "
                   "stopping power below beta*gamma=0.05 (KE 3.5 MeV for a triton, ~2.4 for a "
                   "deuteron), which is most of this channel's backward arm.\033[0m\n";
   }

   if (seedFromSpyral) {
      fitter->SetSeedFromSpyral(kTRUE);
      std::cout << "  SEED FROM SPYRAL: circle on the first arc, polar from a rho-vs-z regression\n";
   }
   if (rangeConstraint) {
      if (!(rangeDensity > 0))
         std::cout << "\033[1;31mWARNING: rangeConstraint on with density " << rangeDensity
                   << " -- the range<->energy table needs a real gas density, constraint will be inert.\033[0m\n";
      fitter->SetRangeConstraint(kTRUE, rangeDensity, rangeMatA);
      std::cout << "  RANGE CONSTRAINT on: rho = " << rangeDensity << " g/cm3, matA = " << rangeMatA << "\n";
   }

   fitter->SetMatEffectsFallback(matFallback);

   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "  \033[1;31m(applied at fMinPoints=30, NOT mp15)\033[0m\n";
      }
   } else {
      std::cout << "  PID gate: NONE (fit all tracks; apply the hand-drawn gate afterwards)\n";
   }

   AtFitterTask *fitterTask = new AtFitterTask(std::move(fitter));
   fitterTask->SetInputBranch("AtPatternEvent");
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   run->AddTask(pidTask);
   run->AddTask(fitterTask);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   t.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << t.RealTime() << " s)\n";
   if (rangeConstraint && fitterObs)
      std::cout << "  range constraint: contained " << fitterObs->GetNRangeContained() << ", applied "
                << fitterObs->GetNRangeApplied() << ", failed the Bragg test " << fitterObs->GetNRangeBraggFail()
                << "\n";
}
