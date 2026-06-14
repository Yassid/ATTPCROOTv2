/// @file fitGenfitter_a1975.C
/// @brief Fit a pre-reconstructed a1975 file with the CLEAN EventFit::AtGenfitter
/// (modern AtFitter base + AtFitterTask, geometric closest-to-axis seed, no merging).
/// Reads <run>_reco.root (AtPatternEvent), writes <run>_genfitter<suffix>.root with
/// AtTrackingEvent + AtPIDEvent. Sibling of fitGenfit_a1975.C (the legacy fallback).
///
///   root -b -q 'pipeline/fitGenfitter_a1975.C("run_0106", 500, "/mnt/f/a1975/reco/", "_test", "/tmp/")'

void fitGenfitter_a1975(TString fileName = "run_0106", Long64_t nEvents = -1, TString ioDir = "/mnt/f/a1975/reco/",
                        TString outSuffix = "", TString outDir = "", Double_t bField = -2.85, Int_t minIter = 2,
                        Int_t maxIter = 5, TString elossName = "deuteron_H2_catima.txt",
                        TString pidGate = "pid/deuteron_band.json", Double_t measSigma = 4.0,
                        Bool_t matEffects = kFALSE, Int_t pdg = 1000010020, Double_t massAmu = 2.0135532,
                        Int_t Z = 1, Bool_t mergeContinuity = kFALSE)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outBase = (outDir.Length() ? outDir : ioDir);
   TString outputFile = outBase + fileName + "_genfitter" + outSuffix + ".root";
   // empty elossName -> empty eloss file: rely on genfit's internal Bethe-Bloch (material effects)
   std::string elossFile = elossName.Length() ? ((std::string)dir.Data() + "/resources/energy_loss/" + elossName.Data()) : std::string("");

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   std::cout << "\033[1;33m=== fitGenfitter_a1975 (clean) ===\033[0m\n  in  : " << inputFile << "\n  out : "
             << outputFile << "\n  B=" << bField << " iter " << minIter << "-" << maxIter << "\n";

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

   // particle hypothesis from args (6th ctor arg is noMatEffects -> invert matEffects)
   std::cout << "  particle: pdg=" << pdg << " massAmu=" << massAmu << " Z=" << Z << "\n";
   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, elossFile, !matEffects, minIter,
                                                         maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   if (mergeContinuity) {
      fitter->SetMergeContinuity(kTRUE); // merge PRA-split fragments of one track before fitting
      std::cout << "  continuity merging: ON\n";
   }
   if (pidGate.Length()) {
      if (gSystem->AccessPathName(pidGate.Data()))
         std::cout << "\033[1;31mWARNING: PID gate " << pidGate << " not found; gating will be disabled.\033[0m\n";
      else {
         fitter->SetPIDGate(pidGate.Data());
         std::cout << "  PID gate: " << pidGate << "\n";
      }
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
}
