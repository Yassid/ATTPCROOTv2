/// @file fitGenfit_C15d.C
/// @brief 15C + d  --  stage 2: GENFIT fit of a C15d reco, with the CATIMA material model.
///
/// Reads <run>_reco.root (AtPatternEvent), writes <run>_genfit_<tag><suffix>.root carrying
/// AtTrackingEvent + AtPIDEvent.
///
///   root -b -q 'fitGenfit_C15d.C("run_0017", 500, "/home/yassid/C15d_reco/", "_test", "/tmp/")'
///
/// ---------------------------------------------------------------------------------------
/// WHY THE DEFAULTS ARE WHAT THEY ARE
///
/// matEffects = kTRUE and catimaMSC = catimaStraggling = catimaELoss = kTRUE. This is the
/// whole point of the workspace, so it is the default rather than an opt-in. Three things
/// follow from that and are easy to get wrong:
///
///   1. The CATIMA flags are INERT with matEffects = kFALSE -- noiseCoulomb and
///      noiseBetheBloch are never called. An A/B against the matEffects-off configuration
///      is therefore NOT a test of CATIMA; it reads as a perfect null and looks like the
///      backend is broken when it was simply never reached.
///   2. MSC and straggling must BOTH be on. Either one alone leaves genfit's own model in
///      the loop for the other term, and the fit-collapse rate goes back up.
///   3. They need a GenFit built with -DGENFIT_USE_CATIMA=ON. The calls compile and run
///      against a stock GenFit but do nothing. Check with
///        ldd $GENFIT/lib/libgenfit2.so | grep catima
///      This installation: ~/fair_install/GenFitInst, branch catima-scattering, linked
///      against ~/fair_install/catima-inst.
///
/// matFallback = kFALSE. With it on, a track whose material-effects fit throws is silently
/// refitted with setNoEffects(true) and KEPT, so the sample quietly mixes matFX and no-matFX
/// tracks. Off, a failed fit simply drops out and the production stays one thing.
///
/// Backward tracks are kept (theta window 5-178 deg, backwardSeedFix ON): in (d,p) inverse
/// kinematics the proton goes largely backward in the lab, and the default 10-170 window
/// plus forward-biased seeding throws most of them away.
///
/// bField SIGN IS UNVERIFIED for this run set. Experimental handedness elsewhere on this
/// detector is Bz = -2.85; it is an argument so both can be tried. Pick it on the fitted
/// vertex (should sit near the beam axis) and a physical ejectile KE, not on chi2 alone.
/// ---------------------------------------------------------------------------------------

void fitGenfit_C15d(TString fileName = "run_0017", Long64_t nEvents = -1,
                    TString ioDir = "/home/yassid/C15d_reco/", TString outSuffix = "", TString outDir = "",
                    Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5, TString pidGate = "",
                    Double_t measSigma = 4.0, Double_t thetaMinDeg = 5.0, Double_t thetaMaxDeg = 178.0,
                    Bool_t matEffects = kTRUE, Bool_t backwardSeedFix = kTRUE,
                    // ejectile: proton = the (d,p) channel Spyral solved with p_gate.
                    //   (d,t) : pdg 1000010030, massAmu 3.01550072, Z 1, speciesTag "t"
                    //   (d,d) : pdg 1000010020, massAmu 2.01410178, Z 1, speciesTag "d"
                    Int_t pdg = 2212, Double_t massAmu = 1.00782503207, Int_t Z = 1, TString speciesTag = "p",
                    TString recoSuffix = "_reco", TString geoName = "ATTPC_D300torr_v2_geomanager.root",
                    TString parName = "ATTPC.C15d_D2_300torr.par",
                    // D2 target: A = 2, and rho matches media.geo TargetD2_300 exactly.
                    Double_t gasDensity = 6.5643e-5, Int_t matA = 2, Bool_t matFallback = kFALSE,
                    Bool_t catimaMSC = kTRUE, Bool_t catimaStraggling = kTRUE, Bool_t catimaELoss = kTRUE,
                    Bool_t catimaELossFull = kFALSE, Bool_t backExtrap = kFALSE, Bool_t rangeConstraint = kFALSE,
                    Bool_t seedFromSpyral = kFALSE,
                    // Per-run dE/dx gain matching of the persisted AtPIDEvent. Opt-in task; see the
                    // caveat printed below about SetPIDGate.
                    Bool_t gainMatch = kTRUE, TString gainTable = "")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   if (dir.Length() == 0) {
      std::cout << "\033[1;31mERROR: VMCWORKDIR unset -- source build/config.sh first.\033[0m\n";
      return;
   }
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/" + geoName;
   // The par MUST be the one the reco ran with: nothing cross-checks the two, so a reco made
   // with a variant drift velocity fitted against this one fails silently and plausibly.
   TString digiParFile = dir + "/parameters/" + parName;
   TString inputFile = ioDir + fileName + recoSuffix + ".root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfit_" + speciesTag + outSuffix + ".root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(geoManFile.Data())) {
      std::cout << "\033[1;31mERROR: geometry not found: " << geoManFile << "\033[0m\n";
      return;
   }

   std::cout << "\033[1;33m=== fitGenfit_C15d (15C + d, D2 300 torr, ejectile '" << speciesTag << "' pdg=" << pdg
             << " m=" << massAmu << " amu Z=" << Z << ") ===\033[0m\n"
             << "  in  : " << inputFile << "\n  out : " << outputFile << "\n  geo : " << geoName
             << "\n  par : " << parName << "\n  B=" << bField << " T   iter " << minIter << "-" << maxIter
             << "   theta[" << thetaMinDeg << "," << thetaMaxDeg << "]\n";

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

   // Empty eloss file -> genfit's internal model; CATIMA replaces it per step below.
   auto fitter =
      std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, "", !matEffects, minIter, maxIter);

   if (catimaMSC || catimaStraggling) {
      fitter->SetCatimaMaterial(catimaMSC, catimaStraggling);
      std::cout << "  \033[1;35mCATIMA material model: MSC " << (catimaMSC ? "ON" : "off") << ", straggling "
                << (catimaStraggling ? "ON" : "off") << "\033[0m\n";
      if (!matEffects)
         std::cout << "\033[1;31m  WARNING: matEffects is OFF -- the CATIMA material flags are INERT.\033[0m\n";
      if (catimaMSC != catimaStraggling)
         std::cout << "\033[1;31m  WARNING: MSC and straggling disagree; genfit's own model still handles "
                      "the other term.\033[0m\n";
   }

   fitter->SetCatimaELoss(catimaELoss, catimaELossFull);
   if (catimaELoss) {
      if (!matEffects)
         std::cout << "\033[1;31m  WARNING: catimaELoss set but matEffects is OFF -- it is inert.\033[0m\n";
      else
         std::cout << "  \033[1;32mdE/dx from CATIMA"
                   << (catimaELossFull ? " over the FULL range (Bethe-Bloch replaced too)"
                                       : " below beta*gamma=0.05 (Bethe-Bloch kept above)")
                   << ", per-step material\033[0m\n";
   }

   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetSeedFromSpyral(seedFromSpyral);
   fitter->SetMatEffectsFallback(matFallback);
   if (matEffects && !matFallback)
      std::cout << "  \033[1;33mmatFX fallback DISABLED: failed material-effects fits are dropped, not "
                   "retried without material effects\033[0m\n";
   if (rangeConstraint) {
      fitter->SetRangeConstraint(kTRUE, gasDensity, matA);
      std::cout << "  \033[1;32mRANGE CONSTRAINT ON: stopping tracks get a FullMeasurement on |p| from "
                   "their path length\033[0m\n";
   }
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);
   fitter->SetBackwardSeedFix(backwardSeedFix);
   fitter->SetBackExtrapToAxis(backExtrap);
   if (backwardSeedFix)
      std::cout << "  \033[1;36mBackward seed-fix: ON\033[0m\n";

   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "\n";
      }
   } else {
      std::cout << "  PID gate: NONE (fit all tracks; build the gate from the PID plane afterwards)\n";
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

   // Opt-in per-run gain matching of the PID dE/dx, the counterpart of Spyral's GainMatchPhase.
   // Must be added AFTER AtPIDTask -- it rescales that task's output in place.
   AtGainMatchTask *gainTask = nullptr;
   if (gainMatch) {
      TString tbl = gainTable.Length() ? gainTable : (dir + "/macro/Unpack_HDF5/C15d/gainmatch_C15d.csv");
      const Int_t runNo = AtGainMatchTask::RunNumberFromName(fileName);
      if (runNo < 0) {
         std::cout << "\033[1;31mERROR: cannot parse a run number from '" << fileName
                   << "' -- gain matching would use the wrong factor. Aborting.\033[0m\n";
         return;
      }
      gainTask = new AtGainMatchTask(tbl.Data(), runNo);
      std::cout << "  \033[1;32mGain match: ON  (run " << runNo << ", table " << tbl << ")\033[0m\n";
      // AtGenfitter evaluates SetPIDGate with its OWN internal AtSpyralPID; it does not read the
      // AtPIDEvent branch. So an in-fit gate sees UNMATCHED dE/dx while the persisted plane is
      // matched -- two different scales, and the gate silently selects the wrong band.
      if (pidGate.Length())
         std::cout << "\033[1;31m  WARNING: gainMatch AND an in-fit PID gate are both set. The fitter's gate "
                      "is evaluated on UNMATCHED dE/dx (its own AtSpyralPID), while the persisted AtPIDEvent "
                      "IS matched. Either draw that gate on unmatched values, or drop the in-fit gate and "
                      "gate downstream on the matched plane.\033[0m\n";
   }

   run->AddTask(pidTask);
   if (gainTask)
      run->AddTask(gainTask);
   run->AddTask(fitterTask);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   t.Stop();
   std::cout << "\n\033[1;32mDone.\033[0m " << outputFile << "  (Real " << t.RealTime() << " s)\n";
}
