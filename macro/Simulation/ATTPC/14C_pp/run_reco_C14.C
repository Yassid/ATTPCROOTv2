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

/// PAD PLANE IS SWITCHABLE (added 2026-08-28 for the field x pitch study). padSize_mm <= 0, the
/// default, keeps the real AT-TPC pad plane from Lookup20150611.xml and every existing caller
/// therefore behaves exactly as before. A positive value builds a uniform square pad plane of
/// that pitch with AtTpcSquareMap, sized to cover activeExtent_mm -- 2 mm over 500 mm gives
/// 250 x 250 = 62500 pads against the AT-TPC's ~10240. Everything downstream goes through the
/// AtMap interface, the beam-hole inhibition included, so nothing else changes.
///
/// !! THE HDBSCAN PARAMETERS ARE NOT RE-TUNED FOR A DIFFERENT PITCH !! epsilon = 10 mm and
/// minClusterSize = 20 were set for 8 x 12 mm pads. At 2 mm pitch a track deposits several times
/// as many hits and 10 mm spans five pads instead of one, so the finder sees a much denser cloud
/// under an unchanged neighbourhood definition. Compare two pitches only after checking that
/// pattern recognition is not the thing that changed (cluster_eval_C14.C).
///
/// @param persistRaw  write the AtRawEvent pad traces. kTRUE reproduces the historical behaviour
///                    and costs ~4.7 GB per 8000 events; nothing downstream of PSA reads them, so
///                    a production run that only needs AtPatternEvent should pass kFALSE.
void run_reco_C14(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                  TString paramFile = "ATTPC.a1954_C14.par", Double_t thr = 20, Int_t hdMcs = 20, Int_t hdMs = 8,
                  Int_t nEvents = 0, Double_t holeR = 30.0, TString joinMethod = "mover",
                  Double_t padSize_mm = -1.0, Double_t activeExtent_mm = 500.0, Bool_t persistRaw = kTRUE,
                  Bool_t useArcWalk = kFALSE, Int_t arcWalkKNN = 10)
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

   std::shared_ptr<AtMap> mapping;
   if (padSize_mm > 0) {
      int nPad = (int)std::lround(activeExtent_mm / padSize_mm);
      auto sq = std::make_shared<AtTpcSquareMap>(padSize_mm, nPad, nPad);
      sq->GeneratePadPlane();
      mapping = sq;
      std::cout << "PADS  : square " << padSize_mm << " mm, " << nPad << " x " << nPad << " = " << nPad * nPad
                << " pads over " << activeExtent_mm << " mm" << std::endl;
   } else {
      auto at = std::make_shared<AtTpcMap>();
      at->ParseXMLMap(mapParFile.Data());
      at->GeneratePadPlane();
      mapping = at;
      std::cout << "PADS  : AT-TPC pad plane from " << scriptfile << std::endl;
   }

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
   pulse->SetPersistence(persistRaw);
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
   // arc-walk orders clusters by a kNN walk over the hits instead of by the drift coordinate.
   // Meant for tracks nearly perpendicular to the beam, where the z advance per cluster drops
   // below the z noise. NOTE: AtGenfitter must ALSO be told to keep that order
   // (SetUseClusterOrder), otherwise it re-sorts the measurements by z and discards it.
   if (useArcWalk) {
      praTask->SetUseArcWalk(kTRUE);
      praTask->SetArcWalkKNN(arcWalkKNN);
   }

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
