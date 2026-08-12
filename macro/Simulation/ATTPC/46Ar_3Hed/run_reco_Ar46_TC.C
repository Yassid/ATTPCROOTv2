/// @file run_reco_Ar46_TC.C
/// @brief Digitise + reconstruct the 46Ar(3He,d)47K simulation. Stops at AtPatternEvent.
///
/// Mirrors run_reco_C16_TC.C (a1975), which was itself matched stage by stage to the a1975 data
/// chain. Nothing about the chain is re-tuned here -- only the gas parameters change:
///     AtPSAMax(threshold 20)
///     no cleaning
///     AtTrackFinderTC(clusterRadius 15.0, clusterDistance 7.5)
///
/// WHY AtTrackFinderTC AND NOT HDBSCAN. On a1975, reconstructing simulation with HDBSCAN at
/// minClusterSize = 20 discarded short tracks by construction and collapsed the acceptance
/// exactly where the data was richest. That argument applies here with more force, not less: at
/// theta_cm = 15 deg the deuteron carries 5.1 MeV and its helix is only 24 cm across, so it is a
/// short, tightly curled track -- the kind a minimum-cluster-size cut removes first.
///
/// NO FIT, BY DESIGN. The chain ends at AtPatternEvent and the PID observables computed from it
/// (pidPass_Ar46.C). The proton/deuteron separation gets drawn by hand on that plane before any
/// fitting is done, so nothing here selects a species.
///
/// GEOMETRY OF THE EJECTILE, worth knowing before reading the reconstruction. The deuteron helix
/// diameter runs 24 cm at theta_cm = 15 deg to 93 cm at 80 deg, against a 50 cm chamber. Forward
/// CM curls up entirely inside the volume; backward CM leaves radially after part of a turn.
/// Both are real, and the turn-over between them is the acceptance this simulation exists to
/// measure -- so no angular preselection is applied anywhere in this chain.
///
///   root -l 'run_reco_Ar46_TC.C("./data/test_gs.root","./data/test_gs_reco.root")'

#include "../AtBeamHole.h"

/// @param praType  "tc" (default, matches the a1975 data chain) or "hdbscan" for comparison only.
/// @param hdMcs    HDBSCAN min_cluster_size, ignored for "tc".
void run_reco_Ar46_TC(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                      TString paramFile = "ATTPC.46Ar_3Hed_sim.par", Double_t thr = 20, Int_t nEvents = 0,
                      Double_t holeR = 20.0, Double_t clusterRadius = 15.0, Double_t clusterDistance = 7.5,
                      TString praType = "tc", Int_t hdMcs = 20)
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
   // Without it the sim digitises the beam over a region the real detector does not read out, and
   // the clustering merges the beam into the ejectile track. For a1954 that cost 67 percent of
   // the recoil protons until it was inhibited.
   //
   // 20 mm is the a1975 measured radius, CARRIED OVER rather than measured for this setup -- there
   // is no 46Ar data to measure it in. It matters more here than it did there: 46Ar is Z = 18, so
   // the beam deposits (18/6)^2 = 9 times the charge per unit length that the a1975 16C beam did,
   // and whatever leaks out of the inhibited region is that much brighter next to a 5 MeV
   // deuteron. If the deuteron track is being eaten near the axis, this radius is the first thing
   // to vary.
   AtSim::InhibitBeamHole(mapping, holeR);

   // ---- digitization ------------------------------------------------------
   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   // The pad traces are not written: nothing downstream of PSA reads them, and they were the bulk
   // of a multi-GB output file per seed. Persistence kFALSE only drops the branch; the task still
   // hands AtRawEvent to PSA in memory.
   pulse->SetPersistence(kFALSE);
   // Propagate MC truth (A, Z, trackID) onto every hit. This is what lets the deuteron be
   // identified from truth INDEPENDENTLY of the hand-drawn PID gate, which is the only way to tell
   // whether the gate is keeping what it should.
   pulse->SetSaveMCInfo();

   // ---- PSA : AtPSAMax, as the data uses ----------------------------------
   auto psa = new AtPSAMax();
   psa->SetThreshold(thr);
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   // ---- NO cleaning -------------------------------------------------------
   // The a1975 data chain has no cleaning stage; AtDirDeDxCleaner was removing hits the data
   // keeps. PRA runs straight off AtEventH.

   // ---- pattern recognition -----------------------------------------------
   std::unique_ptr<AtPATTERN::AtPRA> praAlgo;
   if (praType == "hdbscan") {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
      p->SetMinClusterSize(hdMcs);
      p->SetMinSamples(3);
      p->SetClusterSelectionEpsilon(10.0);
      p->SetJoinMethod("mover");
      p->SetMinClusterSizeJoin(15);
      p->SetCircleOverlapRatio(0.25);
      p->SetMotionGapTol(40);
      p->SetMotionAngleTol(35);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderHDBSCAN mover (mcs " << hdMcs << ")  -- COMPARISON ONLY" << std::endl;
   } else {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderTC>();
      p->SetClusterRadius(clusterRadius);
      p->SetClusterDistance(clusterDistance);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderTC (as the a1975 data chain)" << std::endl;
   }

   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch("AtEventH");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
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
