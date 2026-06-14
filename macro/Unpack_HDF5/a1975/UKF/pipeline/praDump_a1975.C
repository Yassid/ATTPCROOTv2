/// @file praDump_a1975.C
/// @brief Emit a deterministic per-event signature of the LEGACY AtPRAtask PRA output
/// (track count, total cluster count, and a position checksum over all clusters), so the
/// same dump can be produced on the pre-refactor and post-refactor builds and diffed to
/// prove the AtPRAtask injection refactor changed nothing on the legacy path.
/// Uses only the legacy API (new AtPRAtask()), so it compiles on BOTH builds.
///   root -l -b -q 'praDump_a1975.C("run_0106", 300, "/mnt/f/a1975/reco/", "/tmp/pra_new.txt")'

void praDump_a1975(TString runName = "run_0106", Long64_t nEv = 300, TString recoDir = "/mnt/f/a1975/reco/",
                   TString outTxt = "/tmp/pra_dump.txt")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = getenv("VMCWORKDIR");
   TString inputFile = recoDir + runName + "_reco.root";
   TString outputFile = "/tmp/" + runName + "_pradump.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);

   auto *pra = new AtPRAtask();
   pra->SetInputBranch("AtEventCorrected");
   pra->SetOutputBranch("AtPatternEvent");
   pra->SetPersistence(true);
   pra->SetClusterRadius(15.0);
   pra->SetClusterDistance(7.5);
   run->AddTask(pra);
   run->Init();
   run->Run(0, nEv);

   // read back and write a deterministic signature per event
   TFile *f = TFile::Open(outputFile);
   auto *t = (TTree *)f->Get("cbmsim");
   auto *arr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &arr);
   FILE *out = fopen(outTxt.Data(), "w");
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)arr->At(0);
      int nTr = pe ? pe->GetTrackCand().size() : 0;
      int nClus = 0;
      double cks = 0.0; // position checksum over all clusters (deterministic: TriplClust + ClusterizeSmooth3D)
      double rsum = 0.0; // sum of PRA circle radii (may jitter if the radius fit uses RANSAC)
      if (pe)
         for (auto &tr : pe->GetTrackCand()) {
            auto *cl = tr.GetHitClusterArray();
            nClus += cl->size();
            for (auto &c : *cl) {
               auto p = c.GetPosition();
               cks += p.X() + p.Y() + p.Z();
            }
            rsum += tr.GetGeoRadius();
         }
      fprintf(out, "ev %lld: tracks=%d clusters=%d poschk=%.4f radsum=%.4f\n", i, nTr, nClus, cks, rsum);
   }
   fclose(out);
   std::cout << "wrote signature -> " << outTxt << "\n";
}
