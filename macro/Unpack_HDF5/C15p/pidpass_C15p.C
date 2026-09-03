/// @file pidpass_C15p.C
/// @brief 15C + p  --  add an AtPIDEvent branch to the EXISTING a2091 reconstruction.
///
///   root -b -q 'pidpass_C15p.C("run_0138")'
///
/// ★ WHY THIS EXISTS AND WHY THERE IS NO RECO STAGE HERE.
/// The a2091 proton runs were already reconstructed by the older ATTPCROOTv2 tree -- 41 runs,
/// 86 GB, in /home/yassid/a2091_C15_reco. Those files were checked against THIS build before
/// anything was ported: run_0138 gives 333 tracks and 95381 hits over 300 events, so there is no
/// cross-branch ClassDef mismatch and they can be reused as they stand. That saves the entire
/// reconstruction stage, which is the expensive one.
///
/// What they do NOT carry is AtPIDEvent -- the old pipeline computed PID in its own analysis
/// macros instead of persisting it. Everything downstream in this workspace (pidntuple, the PID
/// plane, the gain match, the gates) reads AtPIDEvent, exactly as the C15d workspace does. So this
/// macro runs ONLY AtPIDTask over the existing AtPatternEvent and writes a small companion file.
///
/// ⚠ The output is a NEW file, <run>_pid_reco.root, NOT an edit of the 86 GB original. The
/// originals belong to the older analysis and are left untouched.
///
/// ⚠ gainMatch is deliberately absent: as in C15d the reco persists RAW dE/dx and the per-run
/// gain is measured from this workspace's own plane, applied at read time. Baking it in here would
/// freeze a number that has not been determined yet.

#define cRED "\033[1;31m"
#define cGREEN "\033[1;32m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"

void pidpass_C15p(TString fileName = "run_0138", Long64_t nEvents = -1,
                  TString inDir = "/home/yassid/a2091_C15_reco/",
                  TString outDir = "/home/yassid/C15p_reco/",
                  TString parName = "ATTPC.a2091_C15.par",
                  TString geoName = "ATTPC_H300torr_RT_geomanager.root",
                  // 30 = the same min-points as C15d, so the two planes are cut the same way and
                  // a gate drawn on one is at least comparable to the other. WHATEVER VALUE A
                  // PRODUCTION PERSISTS IS THE PLANE ITS GATES MUST BE DRAWN ON.
                  Int_t pidMinPoints = 30, Double_t pidZTieTol = 0.0, Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   if (dir.Length() == 0) {
      std::cout << cRED << "ERROR: VMCWORKDIR unset -- source build/config.sh first." << cNORMAL << std::endl;
      return;
   }
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());

   TString inputFile = inDir + fileName + "_reco.root";
   TString outputFile = outDir + fileName + "_reco.root";
   TString digiParFile = dir + "/parameters/" + parName;
   TString geoManFile = dir + "/geometry/" + geoName;

   std::cout << cYELLOW << "=== pidpass_C15p (15C + p, H2 300 torr) ===" << cNORMAL << "\n"
             << "  in  : " << inputFile << "\n"
             << "  par : " << digiParFile << "\n"
             << "  geo : " << geoManFile << "\n"
             << "  out : " << outputFile << "\n";

   for (const TString &f : {inputFile, digiParFile, geoManFile})
      if (gSystem->AccessPathName(f.Data())) {
         std::cout << cRED << "ERROR: not found: " << f << cNORMAL << std::endl;
         return;
      }
   gSystem->mkdir(outDir, kTRUE);

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   // ★ PAR SANITY CHECK, same as C15d. A standalone comment line inside [AtDigiPar] makes
   // FairParAsciiFileIo stop parsing there and everything after reads as 0 -- which shows up as
   // theta exactly 90 deg and absurd energies rather than as an error. Refuse to run instead.
   {
      auto *dp = dynamic_cast<AtDigiPar *>(rtdb->getContainer("AtDigiPar"));
      rtdb->initContainers(0);
      if (dp == nullptr || dp->GetTBTime() <= 0 || dp->GetDriftVelocity() <= 0 || dp->GetTBEntrance() <= 0 ||
          dp->GetZPadPlane() <= 0) {
         std::cout << cRED << "ERROR: " << digiParFile << " did not parse.\n"
                   << "  TBTime=" << (dp ? dp->GetTBTime() : -1) << " dv=" << (dp ? dp->GetDriftVelocity() : -1)
                   << " TBEntrance=" << (dp ? dp->GetTBEntrance() : -1)
                   << " ZPadPlane=" << (dp ? dp->GetZPadPlane() : -1) << "\n"
                   << "  Check for standalone comment lines INSIDE the [AtDigiPar] block." << cNORMAL << std::endl;
         return;
      }
      std::cout << cGREEN << "  par OK: " << dp->GetTBTime() << " ns, dv " << dp->GetDriftVelocity()
                << " cm/us, entrance TB " << dp->GetTBEntrance() << cNORMAL << "\n";
   }

   AtPIDTask *pidTask = new AtPIDTask();
   pidTask->SetInputBranch("AtPatternEvent");
   pidTask->SetOutputBranch("AtPIDEvent");
   pidTask->SetBField(bField);
   pidTask->SetMinPoints(pidMinPoints);
   pidTask->SetZTieTolerance(pidZTieTol);
   pidTask->SetPersistence(kTRUE);
   run->AddTask(pidTask);
   std::cout << "  PID  : minPoints=" << pidMinPoints << ", zTieTol=" << pidZTieTol << ", B=" << bField << "\n";

   run->Init();
   run->Run(0, nEvents > 0 ? nEvents : 0);

   timer.Stop();
   std::cout << "\n" << cGREEN << "Done." << cNORMAL << " " << outputFile << "\n"
             << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}
