/// @file run_pr_only.C
/// @brief Run pattern recognition only on already-digitized+PSA data.
///
/// Reads output with AtEventH (from run_digi_noPR.C) and runs AtPRAtask.
/// Allows quick scanning of PR parameters without re-digitizing.
///
/// Run: root -b -q 'run_pr_only.C(1000, 6.0)' // 1000 events, t=6.0

void run_pr_only(int nEvents = 10000, float tCluster = 4.0)
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_psa.root";
   TString outputFile = inOutDir + "output_digi.root"; // overwrite with new PR
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/Lookup20150611.xml";

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   AtPRAtask *praTask = new AtPRAtask();
   praTask->SetTcluster(tCluster);
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(praTask);
   fRun->Init();

   timer.Start();
   fRun->Run(0, nEvents);
   timer.Stop();

   std::cout << "PR done: t=" << tCluster << ", " << timer.RealTime() << " s" << std::endl;
}
