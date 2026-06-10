/// @file runPID_a1975.C
/// @brief Run AtPIDTask on a pre-reconstructed a1975 file, writing both PID methods.
///
/// Reads <run>_reco.root (AtPatternEvent) and runs AtPIDTask, which computes PID
/// observables for every track with BOTH the legacy AtPIDEstimator (charge/arclen)
/// and the Spyral-style AtSpyralPID (first-arc + spline). Output <run>_pid.root has
/// an AtPIDEvent branch (parallel per-track vectors GetClassic()/GetSpyral()) so the
/// two methods can be compared at full statistics.
///
/// Note: brho/dEdx are magnitude-based (sign-agnostic), so B is set positive here —
/// no handedness flip is needed for PID (unlike the UKF fit).
///
/// Run: root -b -q 'runPID_a1975.C("run_0116")'

void runPID_a1975(TString fileName = "run_0116", Long64_t nEvents = -1, Double_t bField = 2.85, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = outDir + fileName + "_reco.root";
   TString outputFile = outDir + fileName + "_pid.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m" << std::endl;
      return;
   }

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetBField(bField);
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetPersistence(kTRUE);
   fRun->AddTask(pidTask);

   fRun->Init();
   TStopwatch timer;
   timer.Start();
   fRun->Run(0, nEvents < 0 ? 0 : nEvents);
   timer.Stop();

   std::cout << "\n\033[1;32mDone.\033[0m  Output: " << outputFile << "  (Real " << timer.RealTime() << " s)"
             << std::endl;
}
