// PRA-only re-reconstruction: read AtEventCorrected from an existing reco file and
// re-run ONLY AtPRAtask (AtTrackFinderTC) with chosen triplet params. Skips the
// expensive unpack/PSA/SC. Output is a drop-in <fileName><tag>_reco.root that the
// genfit fitter can consume (AtEventCorrected passthrough + fresh AtPatternEvent).
//
//   root -l -b -q 'reReco_PRA.C("run_0016","loose",0.20,12.0,-1)'
//   root -l -b -q 'reReco_PRA.C("run_0016","def",0.03,4.0,-1)'
void reReco_PRA(TString fileName = "run_0016", TString tag = "loose",
                float aTrip = 0.20, float tClus = 12.0, Long64_t nEvents = -1,
                TString ioDir = "/mnt/f/a1975/reco_d2/") {
   TString inputFile = ioDir + fileName + "_reco.root";
   TString outputFile = ioDir + fileName + tag + "_reco.root";
   if (gSystem->AccessPathName(inputFile.Data())) {
      printf("\033[1;31mERROR: %s not found\033[0m\n", inputFile.Data()); return;
   }
   printf("re-PRA  in=%s  out=%s  a=%.2f t=%.1f\n", inputFile.Data(), outputFile.Data(), aTrip, tClus);

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);

   // PRA only — same cluster radius/distance as production, tunable triplet a/t
   auto praAlgo = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   praAlgo->SetClusterRadius(15.0);
   praAlgo->SetClusterDistance(7.5);
   praAlgo->SetAtriplet(aTrip);
   praAlgo->SetTcluster(tClus);
   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch("AtEventCorrected");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(true);
   run->AddTask(praTask);

   run->Init();
   TStopwatch t; t.Start();
   run->Run(0, nEvents < 0 ? 0 : nEvents);
   printf("\033[1;32mDone\033[0m %s  (Real %.1fs)\n", outputFile.Data(), t.RealTime());
}
