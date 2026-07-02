// Fit one branch of the A/B reco (proton hypothesis), parametrized input branch.
//   root -b -q 'fitGenfitter_ab.C("run_0016","AtPatternEvent","A","/tmp/abfull/")'
//   root -b -q 'fitGenfitter_ab.C("run_0016","AtPatternEventClean","B","/tmp/abfull/")'
void fitGenfitter_ab(TString fileName, TString inputBranch, TString outTag, TString ioDir = "/tmp/abfull/",
                     Long64_t nEvents = -1, Double_t bField = -2.85, Int_t minIter = 2, Int_t maxIter = 5,
                     Double_t measSigma = 4.0, Double_t thetaMinDeg = 10.0, Double_t thetaMaxDeg = 170.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   const Int_t pdg = 2212; const Double_t massAmu = 1.00782503207; const Int_t Z = 1;

   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";
   TString inputFile = ioDir + fileName + "_ab_reco.root";
   TString outputFile = ioDir + fileName + "_genfit_" + outTag + ".root";
   if (gSystem->AccessPathName(inputFile.Data())) { std::cout << "ERROR: no " << inputFile << "\n"; return; }
   std::cout << "FIT branch " << inputBranch << " -> " << outputFile << "\n";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   if (gROOT->FindObject("FAIRGeom") == nullptr) { TFile *gf = TFile::Open(geoManFile); gf->Get("FAIRGeom"); }

   auto fitter = std::make_unique<EventFit::AtGenfitter>(bField, pdg, massAmu, Z, "", true, minIter, maxIter);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetThetaWindow(thetaMinDeg, thetaMaxDeg);
   fitter->SetBackwardSeedFix(kTRUE);

   AtFitterTask *fitterTask = new AtFitterTask(std::move(fitter));
   fitterTask->SetInputBranch(inputBranch);
   fitterTask->SetOutputBranch("AtTrackingEvent");
   fitterTask->SetFitMetadataBranch("AtFitMetadata");
   fitterTask->SetPersistence(kTRUE);

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch(inputBranch);
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);

   run->AddTask(pidTask);
   run->AddTask(fitterTask);

   TStopwatch t; t.Start();
   run->Init();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   t.Stop();
   std::cout << "\033[1;32mFIT " << outTag << " done\033[0m " << outputFile << "  (" << t.RealTime() << " s)\n";
}
