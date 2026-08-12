/// @file run_reco_C16_TC.C
/// @brief Digitise + reconstruct the a1975 simulation with the SAME track finder as the data.
///
/// !! THE TRACK FINDER MUST MATCH THE DATA, AND HDBSCAN DOES NOT !!
/// The a1975 data reconstruction uses AtTrackFinderTC (triplclust) -- all three pipeline macros,
/// unpackReco_a1975_UKF.C, unpackNFit_a1975_UKF.C and testPRAinject_a1975.C, instantiate it and
/// none uses HDBSCAN. Reconstructing the simulation with HDBSCAN at minClusterSize = 20 discards
/// short tracks by construction, and the recoil proton at forward theta_cm carries 0.5-2 MeV and
/// travels a couple of centimetres. The result was a simulated acceptance that collapsed to zero
/// exactly where the data is richest: only 7 % of simulated tracks fell below 2 MeV against 58 %
/// in the data, and a gain scan over a factor 33 moved that by 2 percentage points, ruling out
/// the digitisation as the cause.
///
/// So this mirrors the data chain instead:
///     AtPSAMax(threshold 20)                                   <-- as in the data
///     no cleaning                                              <-- as in the data
///     AtTrackFinderTC(clusterRadius 15.0, clusterDistance 7.5)  <-- as in the data
///
/// The data also applies a space-charge correction between PSA and PRA (AtEDistortionModel with
/// z/radial/transverse LUTs). That is deliberately NOT reproduced here.
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
///   root -l 'run_reco_C16_TC.C()'

#include "../AtBeamHole.h"

/// @param praType  "tc" (default, and the only one that matches the data) or "hdbscan". The
///                 alternative exists so the two finders can be compared on the SAME simulation
///                 and the same digitisation; it is not an endorsement of running HDBSCAN for
///                 physics here -- see the warning at the top of this file. Leaving praType at its
///                 default reproduces the validated chain exactly.
/// @param hdMcs    HDBSCAN min_cluster_size, ignored for "tc".
void run_reco_C16_TC(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                  TString paramFile = "ATTPC.a1975_C16_sim.par", Double_t thr = 20,
                  Int_t nEvents = 0, Double_t holeR = 20.0, Double_t clusterRadius = 15.0,
                    Double_t clusterDistance = 7.5, TString praType = "tc", Int_t hdMcs = 20)
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
   // The pad traces are not written. Nothing downstream of PSA reads them -- the fit stage takes
   // AtPatternEvent and AtSpyralPID works off the hits' charge -- and they were the bulk of a
   // ~4 GB output file per seed, all of it crossing a 9p mount. Persistence kFALSE only drops the
   // branch from the tree; the task still hands AtRawEvent to PSA in memory, and SetSaveMCInfo
   // below still propagates the MC truth onto the hits, which is what the truth matching needs.
   pulse->SetPersistence(kFALSE);
   // Propagate MC truth (A, Z, trackID) onto every hit. Needed by gate_truth_C14.C to pick
   // out the recoil proton: without it the fitter treats the beam and the scattered 14C as
   // proton candidates too, which is what swamped the first closure attempt (6605 of 6771
   // fits were KE<5 MeV junk at theta~13 deg, only 139 real protons survived).
   pulse->SetSaveMCInfo();

   // ---- PSA : AtPSAMax, which is what the DATA uses -----------------------
   // The a1975 data reconstruction runs AtPSAMax with threshold 20 and nothing else set
   // (unpackReco_a1975_UKF.C:69-70). It was AtPSAMultiFit here, inherited from the a1954/(p,d)
   // macros, and that is not a cosmetic difference: MultiFit fits a sum of GET response pulses
   // and deconvolves overlapping peaks, so it produces a different hit multiplicity and a
   // different z per hit -- it reports the model maximum rather than t0, which carries an ~11 mm
   // z shift. Reconstructing the simulation with a different PSA from the data means the
   // acceptance describes a detector the data was never analysed with.
   auto psa = new AtPSAMax();
   psa->SetThreshold(thr);
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   // ---- NO cleaning ------------------------------------------------------
   // The data chain has no cleaning stage at all; AtDirDeDxCleaner was inherited from the a1954
   // macros and was removing hits the data keeps. PRA therefore runs straight off AtEventH.
   //
   // The data DOES apply a space-charge correction (AtEDistortionModel, AtEventH ->
   // AtEventCorrected) which is deliberately not reproduced here.

   // ---- pattern recognition : same as the data ---------------------------
   // Exactly the data's configuration: AtTrackFinderTC carries the standard TC defaults and only
   // the cluster radius and distance are set (unpackReco_a1975_UKF.C:83-92).
   std::unique_ptr<AtPATTERN::AtPRA> praAlgo;
   if (praType == "hdbscan") {
      // The a1954/a2091 "mover" configuration, the one actually in service elsewhere in this repo
      // (unpackReco_C14.C:109-120). Copied rather than re-tuned so the comparison is against a
      // known-good setup and not against a configuration invented for this test.
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
      std::cout << "PRA   : AtTrackFinderTC (as the data)" << std::endl;
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
