/// @brief Quick regression test: run PRA on existing PSA data and report
/// track yields. Used to verify the AtPRA::SelectAndMergeTracks vertex-by-xy
/// fix doesn't break the forward-stopping-proton pipeline.
///
/// Run: root -b -q 'run_pr_regression.C(200)'

void run_pr_regression(int nEvents = 200, float tCluster = 4.0)
{
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");

   TString inOutDir = "./data/";
   TString inputFile = inOutDir + "output_psa.root";
   TString outputFile = inOutDir + "output_pr_regression.root";
   TString paramFile = "ATTPC.e20009_sim.par";

   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;

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

   std::cout << "PR regression done: " << timer.RealTime() << " s\n";

   // Quick analysis of yields
   TFile fOut(outputFile);
   auto *tOut = (TTree *)fOut.Get("cbmsim");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tOut->SetBranchAddress("AtPatternEvent", &patArr);
   int nWith = 0, nTotalTracks = 0;
   double sumR = 0, sumTh = 0;
   int nValidGeo = 0;
   Long64_t n = tOut->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      tOut->GetEntry(i);
      if (patArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patArr->At(0);
      auto &cands = pat->GetTrackCand();
      if (!cands.empty()) {
         ++nWith;
         nTotalTracks += cands.size();
         for (auto &t : cands) {
            double R = t.GetGeoRadius();
            double th = t.GetGeoTheta() * 180. / M_PI;
            if (std::isfinite(R) && R > 0 && std::isfinite(th)) {
               sumR += R;
               sumTh += th;
               ++nValidGeo;
            }
         }
      }
   }
   std::cout << "\n=== 16C+pp PRA regression ===\n";
   std::cout << "Events:               " << n << "\n";
   std::cout << "Events with tracks:   " << nWith << " (" << 100. * nWith / n << "%)\n";
   std::cout << "Total tracks:         " << nTotalTracks << "\n";
   std::cout << "Tracks/evt (avg):     " << (double)nTotalTracks / std::max(1, nWith) << "\n";
   std::cout << "Valid Geo* tracks:    " << nValidGeo << "\n";
   if (nValidGeo > 0) {
      std::cout << "  <GeoRadius> [mm]:  " << sumR / nValidGeo << "\n";
      std::cout << "  <GeoTheta> [deg]:  " << sumTh / nValidGeo << "\n";
   }
}
