/// @file rerun_pra_dp.C
/// @brief Re-runs ONLY cleaning + pattern recognition from the AtEventH stored in a reco file (2026-11-02), so that a finder option can be tested
///        without unpacking or PSA. Settings are those of the production of that file: data dv_dt.sh mcs 0, ms 3; sim run_reco_C16dt.C defaults
///        mcs 20, ms 8. stray = 1 turns on AtPRA::SetRemoveStrayHits. join: HDBSCAN join method ("mover" in production, "none" to test whether
///        the joining attaches the stray hits).
///   root -l -b -q 'rerun_pra_dp.C("<reco>","out.root",0,3,"mover",1,8000)'
void rerun_pra_dp(TString in, TString out, int mcs, int ms, TString join = "mover", int stray = 0, Long64_t nEvents = 0,
                  TString par = "ATTPC.a1975_deuterium_dv1104.par")
{
   gSystem->Load("libAtReconstruction.so");
   TString dir = getenv("VMCWORKDIR");
   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(in));
   run->SetSink(new FairRootFileSink(out));
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/" + par).Data(), "in");
   rtdb->setFirstInput(parIo);
   auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>();
   auto *cleanTask = new AtDataCleaningTask(std::move(cleaner));
   cleanTask->SetInputBranch("AtEventH");
   cleanTask->SetOutputBranch("AtEventClean");
   cleanTask->SetPersistence(kFALSE);
   auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
   p->SetMinClusterSize(mcs);
   p->SetMinSamples(ms);
   p->SetClusterSelectionEpsilon(10.0);
   p->SetJoinMethod(join.Data());
   p->SetMinClusterSizeJoin(15);
   p->SetCircleOverlapRatio(0.25);
   p->SetMotionGapTol(40);
   p->SetMotionAngleTol(35);
   p->SetRemoveStrayHits(stray != 0);
   AtPRAtask *praTask = new AtPRAtask(std::move(p));
   praTask->SetInputBranch("AtEventClean");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);
   run->AddTask(cleanTask);
   run->AddTask(praTask);
   run->Init();
   run->Run(0, nEvents);
   std::cout << "Done " << out << std::endl;
}
