/// @file pidPass_a1975.C
/// @brief Fast PID-only retrofit: read <run>_reco.root (AtPatternEvent) and run
/// AtPIDTask to compute + PERSIST the Spyral PID (brho, sqrt(dEdx) per track) in
/// an AtPIDEvent branch -> <run>_pid.root. Lets analysis read the gate observables
/// back instead of re-running AtSpyralPID::Estimate() every pass.
///
///   root -b -q 'pipeline/pidPass_a1975.C("run_0106","/mnt/f/a1975/reco/","/mnt/f/a1975/reco_gf/")'

void pidPass_a1975(TString fileName = "run_0106", TString inDir = "/mnt/f/a1975/reco/",
                   TString outDir = "/mnt/f/a1975/reco_gf/", Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = inDir + fileName + "_reco.root";
   TString outputFile = outDir + fileName + "_pid.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);

   AtPIDTask *pid = new AtPIDTask();
   pid->SetBField(bField);
   pid->SetInputBranch("AtPatternEvent");
   pid->SetOutputBranch("AtPIDEvent");
   pid->SetPersistence(kTRUE);
   run->AddTask(pid);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, 0);
   t.Stop();
   std::cout << "\033[1;32mDone.\033[0m -> " << outputFile << "  (Real " << t.RealTime() << " s)\n";
}
