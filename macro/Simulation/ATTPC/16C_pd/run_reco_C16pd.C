/// @file run_reco_C16pd.C
/// @brief Digitise + reconstruct the a1975 16C(p,d)15C simulation, so the output can be fed to
/// the same fitter macro the data uses.
///
/// Mirrors pipeline/unpackReco_C14.C as invoked by reco_hdb_C14.sh:
///     AtPSAMultiFit(thr = 20, peakingTime 0.720, maxPeaks 4, minSeparation 4)
///     AtDirDeDxCleaner
///     AtTrackFinderHDBSCAN(mcs = 20, ms = 8, eps 10, mover, join 15, overlap 0.25,
///                          gapTol 40, angleTol 35)
/// It uses parameters/ATTPC.a1975_C16_sim.par. NOTE that this is NOT the par the a1975 data
/// reconstruction loads: that one is ATTPC.a1954.par, which declares 600 torr and a density of
/// 0.197 mg/cm3, neither of which is H2 at 1 bar. The two agree on what matters for the geometry
/// of the reconstruction -- drift velocity 1.30 cm/us, sampling rate, B field, drift length --
/// and differ in the gas description and in the digitisation gain and diffusion, which the data
/// par does not need but the simulation does. That inconsistency in the data par is unresolved
/// and is worth settling before any absolute use of these acceptances.
///
/// CLOSURE TEST, not a drift-velocity measurement. Digitising and reconstructing with the same
/// drift velocity is an identity: z -> TB -> z cancels exactly, so this cannot say whether
/// 1.30 cm/us is right. What it can show is whether the chain returns the generated state at the
/// energy it was generated at, once the beam energy loss along the vertex is accounted for -- and
/// that loss is NOT small here: the truth alone puts the apparent Ex 1.3 MeV higher at the far
/// end of the chamber than at the entrance, which is larger than the 0.74 MeV spacing between the
/// first two states.
///
/// Output: ./data/sim_reco.root (AtPatternEvent), fit it with the data's own fitter:
///   root -l 'fitUKF_C14.C("sim",-1,"proton",+1,2.85,3.553e-5,"","./data/",0.5,0.1,1,10,"./data/")'
///   NOTE bFieldSign = +1 for simulation. The data uses -1 because the sim reverses the
///   drift-z handedness in digitization while the experiment does not.
///
///   root -l 'run_reco_C16pd.C()'

#include "../AtBeamHole.h"

void run_reco_C16pd(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                  TString paramFile = "ATTPC.a1975_C16_sim.par", Double_t thr = 20, Int_t hdMcs = 20, Int_t hdMs = 8,
                  Int_t nEvents = 0, Double_t holeR = 20.0, TString joinMethod = "mover")
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
   // THE a1975 RADIUS IS 20 mm, NOT THE 30 mm USED FOR a1954, and it is measured rather than
   // assumed: in run_0106 the hit density (counts per bin divided by r, since a uniform density
   // gives counts proportional to r) runs 20, 20, 72, 46 over the first 20 mm against 180, 175,
   // 111, 120 beyond it. The suppression sets in below about 20 mm and is near-total below 10.
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
