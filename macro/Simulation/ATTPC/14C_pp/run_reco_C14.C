/// @file run_reco_C14.C
/// @brief Digitize + reconstruct the 14C(p,p') simulation with EXACTLY the chain used on the
/// a1954 data, so the output can be fed to the same fitter macro.
///
/// Mirrors pipeline/unpackReco_C14.C as invoked by reco_hdb_C14.sh:
///     AtPSAMultiFit(thr = 20, peakingTime 0.720, maxPeaks 4, minSeparation 4)
///     AtDirDeDxCleaner
///     AtTrackFinderHDBSCAN(mcs = 20, ms = 8, eps 10, mover, join 15, overlap 0.25,
///                          gapTol 40, angleTol 35)
/// and uses parameters/ATTPC.a1954_C14.par -- the SAME par as the data reconstruction, so
/// the drift velocity (1.30 cm/us), sampling rate and TB entrance all match.
///
/// CLOSURE TEST, not a drift-velocity measurement. Digitizing and reconstructing with the
/// same dv is an identity: z -> TB -> z cancels exactly, so this run cannot say whether
/// 1.30 is the right dv. What it CAN do is show whether the chain returns the ground state
/// at 0 and flat in angle when every input is correct by construction. A residual drift
/// here is a chain artefact, not a dv error -- which is exactly what needs ruling out
/// before blaming dv for the residual drift seen in the data.
///
/// Output: ./data/sim_reco.root (AtPatternEvent), fit it with the data's own fitter:
///   root -l 'fitUKF_C14.C("sim",-1,"proton",+1,2.85,3.553e-5,"","./data/",0.5,0.1,1,10,"./data/")'
///   NOTE bFieldSign = +1 for simulation. The data uses -1 because the sim reverses the
///   drift-z handedness in digitization while the experiment does not.
///
///   root -l 'run_reco_C14.C()'

#include "../AtBeamHole.h"

void run_reco_C14(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                  TString paramFile = "ATTPC.a1954_C14.par", Double_t thr = 20, Int_t hdMcs = 20, Int_t hdMs = 8,
                  Int_t nEvents = 0, Double_t holeR = 30.0, TString joinMethod = "mover")
{
   TString scriptfile = "Lookup20150611.xml";
   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();

   // ---- BEAM HOLE ---------------------------------------------------------
   // a1954 has a 3 cm beam hole. Without it the sim digitizes the beam over a region the real
   // detector does not read out, and HDBSCAN merges the beam into the recoil-proton cluster
   // (67 % of protons -> 0.2 % once inhibited; see cluster_eval_C14.C). Shared helper, because
   // every AT-TPC simulation needs this -- 15C and 16C included. holeR = 0 disables it.
   AtSim::InhibitBeamHole(mapping, holeR);

   // ---- digitization ------------------------------------------------------
   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);
   // Propagate MC truth (A, Z, trackID) onto every hit. Needed by gate_truth_C14.C to pick
   // out the recoil proton: without it the fitter treats the beam and the scattered 14C as
   // proton candidates too, which is what swamped the first closure attempt (6605 of 6771
   // fits were KE<5 MeV junk at theta~13 deg, only 139 real protons survived).
   pulse->SetSaveMCInfo();

   // ---- PSA : same as the data ------------------------------------------
   auto p = new AtPSAMultiFit();
   p->SetThreshold(thr);
   p->SetPeakingTime(0.720);
   p->SetMaxPeaks(4);
   p->SetMinSeparation(4);
   p->SetProminenceOnPrimary(false);
   AtPSAtask *psaTask = new AtPSAtask(p);
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   // ---- cleaning : same as the data --------------------------------------
   auto cleaner = std::make_unique<AtTools::DataCleaning::AtDirDeDxCleaner>();
   auto *cleanTask = new AtDataCleaningTask(std::move(cleaner));
   cleanTask->SetInputBranch("AtEventH");
   cleanTask->SetOutputBranch("AtEventClean");
   cleanTask->SetPersistence(kFALSE);

   // ---- pattern recognition : same as the data ---------------------------
   auto hdb = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
   hdb->SetMinClusterSize(hdMcs);
   hdb->SetMinSamples(hdMs);
   hdb->SetClusterSelectionEpsilon(10.0);
   hdb->SetJoinMethod(joinMethod.Data()); // "mover" (default) | "overlap" | "motion" | "both" | "none"
   hdb->SetMinClusterSizeJoin(15);
   hdb->SetCircleOverlapRatio(0.25);
   hdb->SetMotionGapTol(40);
   hdb->SetMotionAngleTol(35);

   AtPRAtask *praTask = new AtPRAtask(std::move(hdb));
   praTask->SetInputBranch("AtEventClean");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(cleanTask);
   fRun->AddTask(praTask);

   fRun->Init();
   timer.Start();
   if (nEvents > 0)
      fRun->Run(0, nEvents);
   else
      fRun->Run();
   timer.Stop();

   std::cout << "\nsim reco done -> " << outputFile << "   par = " << paramFile << std::endl;
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n" << std::endl;
}
