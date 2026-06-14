/// @file testPRAinject_a1975.C
/// @brief Validate the new AtPRAtask dependency-injection constructor is behaviour-
/// preserving: run PRA twice on the same AtEventCorrected hits — once via the LEGACY
/// path (new AtPRAtask() + SetPRAlgorithm default TC) and once via INJECTION
/// (new AtPRAtask(make_unique<AtTrackFinderTC>())) configured identically — and compare
/// the per-event track counts. They must match exactly.
///   root -l -b -q 'testPRAinject_a1975.C("run_0106", 300, "/mnt/f/a1975/reco/")'

void testPRAinject_a1975(TString runName = "run_0106", Long64_t nEv = 300, TString recoDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = getenv("VMCWORKDIR");
   TString inputFile = recoDir + runName + "_reco.root"; // has AtEventCorrected
   TString outputFile = "/tmp/" + runName + "_prainject.root";
   TString digiParFile = dir + "/parameters/ATTPC.a1954.par";

   const double clusterRadius = 15.0, clusterDistance = 7.5;

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);

   // LEGACY path
   auto *praLegacy = new AtPRAtask();
   praLegacy->SetInputBranch("AtEventCorrected");
   praLegacy->SetOutputBranch("AtPatternEventLegacy");
   praLegacy->SetPersistence(true);
   praLegacy->SetClusterRadius(clusterRadius);
   praLegacy->SetClusterDistance(clusterDistance);

   // INJECTION path — AtTrackFinderTC already carries the standard TC HC defaults
   // (s=0.3,k=19,n=2,m=15,r=2,a=0.03,t=4.0, selectAndMerge on) that the legacy task
   // pushed, so only the cluster radius/distance need setting. This is the CLEAN
   // migration form used by the production macros.
   auto tc = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   tc->SetClusterRadius(clusterRadius);
   tc->SetClusterDistance(clusterDistance);
   auto *praInject = new AtPRAtask(std::move(tc));
   praInject->SetInputBranch("AtEventCorrected");
   praInject->SetOutputBranch("AtPatternEventInject");
   praInject->SetPersistence(true);

   run->AddTask(praLegacy);
   run->AddTask(praInject);
   run->Init();
   run->Run(0, nEv);

   // --- compare the two output branches ---
   TFile *f = TFile::Open(outputFile);
   auto *t = (TTree *)f->Get("cbmsim");
   auto *legArr = new TClonesArray("AtPatternEvent");
   auto *injArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEventLegacy", &legArr);
   t->SetBranchAddress("AtPatternEventInject", &injArr);
   Long64_t n = t->GetEntries();
   int nev = 0, mismatch = 0, legTot = 0, injTot = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *le = (AtPatternEvent *)legArr->At(0);
      auto *in = (AtPatternEvent *)injArr->At(0);
      int lc = le ? le->GetTrackCand().size() : 0;
      int ic = in ? in->GetTrackCand().size() : 0;
      legTot += lc; injTot += ic; nev++;
      if (lc != ic) { mismatch++; if (mismatch <= 5) std::cout << "  MISMATCH ev " << i << ": legacy " << lc << " vs inject " << ic << "\n"; }
   }
   std::cout << "\n=== AtPRAtask injection vs legacy, " << nev << " events ===\n";
   std::cout << "legacy total tracks = " << legTot << ", inject total tracks = " << injTot << "\n";
   std::cout << (mismatch == 0 ? "\033[1;32mIDENTICAL — injection is behaviour-preserving.\033[0m\n"
                               : ("\033[1;31m" + std::to_string(mismatch) + " events differ!\033[0m\n"));
}
