/// @file run_reco_C16dt.C
/// @brief Digitise + reconstruct the a1975 16C(d,t)15C simulation with the SAME chain the data
/// uses, so the acceptance it yields describes the data production rather than a different one.
///
/// WHY THIS EXISTS WHEN A (d,t) ACCEPTANCE ALREADY DOES. The acceptance in the analysis repo was
/// measured with attpc_engine + Spyral. The DATA is reconstructed with ATTPCROOT
/// (AtPSAMultiFit + HDBSCAN) and fitted with genfit. An acceptance measured on a different
/// reconstruction chain corrects for that chain's losses, not for the one being corrected -- the
/// same mismatch already flagged for (p,d), where the simulation ran AtPSAMultiFit + HDBSCAN
/// while the data ran AtPSAMax + TC. This macro closes that gap for (d,t).
///
/// THE PAR IS THE DATA'S OWN, ATTPC.a1975_deuterium_dv1104.par, and that is deliberate: drift
/// velocity 1.10424 cm/us, TBEntrance 560 and ZPadPlane 971.7312 mm must match the production
/// exactly or the z scale of the simulation is not the z scale of the data. Note this is the
/// opposite choice from run_reco_C16pd.C, which uses a sim-only par because the (p,d) data par
/// (ATTPC.a1954.par) describes neither the right gas nor the right pressure.
///
/// *** GAIN IS THE OPEN QUESTION, AND IT IS NOT COSMETIC ***
///     ATTPC.a1975_deuterium_dv1104.par declares Gain = 10000; the (p,d) simulation par declares
///     150000, and that one was clearly tuned (g400000/g1000000/g2500000/g5000000 variants of it
///     exist). Gain sets how many electrons land on a pad, hence whether the pad crosses the PSA
///     threshold at all, hence how many hits a track keeps -- which is precisely what an
///     ACCEPTANCE measures. Too low a gain silently deletes short and low-dE/dx tracks and the
///     acceptance comes out too small, with nothing in the output looking wrong.
///     DO NOT ADOPT AN ACCEPTANCE FROM THIS UNTIL THE GAIN IS VALIDATED against the data: compare
///     hits-per-track and pads-per-track between this output and reco_d2_dv1104. The (p,d) work
///     already found sim/data divergence of that kind (clusters per track 54.4 sim vs 34.5 data),
///     so it is a known failure mode here, not a hypothetical.
///     Gain is scanned by PAR FILE, not by an argument: AtDigiPar exposes GetGain() but no
///     setter, which is also how the (p,d) simulation did it (ATTPC.a1975_C16_sim_g*.par).
///     make_gain_pars.sh writes ATTPC.a1975_deuterium_dv1104_g<N>.par variants that differ from
///     the production par in the Gain line and nothing else.
///
/// CLOSURE, not calibration: digitising and reconstructing at one drift velocity is an identity,
/// z -> TB -> z cancels, so this cannot test whether 1.10424 is right. What it tests is whether
/// the chain returns the state at the energy it was generated at.
///
/// FIELD SIGN: a sim generated at -28.5 kG fits at -2.85, the same as the data, and the fitted
/// polar angle then compares directly with truth. Measured with sim_fit_sign_check.C.
///
///   root -l 'run_reco_C16dt.C("/mnt/f/a1975_C16_dt_sim/gs_s3001_sim.root","/mnt/f/a1975_C16_dt_sim/gs_s3001_reco.root")'

#include "../AtBeamHole.h"

void run_reco_C16dt(TString mcFile = "./data/dt_test.root", TString outputFile = "./data/dt_sim_reco.root",
                    TString paramFile = "ATTPC.a1975_deuterium_dv1104.par", Double_t thr = 20, Int_t hdMcs = 20,
                    Int_t hdMs = 8, Int_t nEvents = 0, Double_t holeR = 20.0, TString joinMethod = "mover",
                    Double_t adcCeiling = 0)
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
   // The real detector does not read out the beam region; without inhibiting it the simulation
   // digitises a beam the data never sees and the clustering merges it into the ejectile track.
   // The a1975 radius is 20 mm, measured from run_0106's hit density rather than assumed.
   AtSim::InhibitBeamHole(mapping, holeR);

   // ---- digitization ------------------------------------------------------
   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   // ADC CEILING (2026-11-02): the GET electronics saturate and the data lose the charge above it; without the clip
   // the simulated high-dE/dx tracks keep charge the data never record. 0 = off (every existing production).
   auto pulseModel = std::make_shared<AtPulse>(mapping);
   if (adcCeiling > 0) {
      pulseModel->SetADCCeiling(adcCeiling);
      std::cout << "ADC ceiling: traces clipped at " << adcCeiling << " ADC" << std::endl;
   }
   AtPulseTask *pulse = new AtPulseTask(pulseModel);
   pulse->SetPersistence(kTRUE);
   // MC truth (A, Z, trackID) on every hit: without it the acceptance cannot tell the triton from
   // the beam or the recoil, and every downstream truth match becomes a guess.
   pulse->SetSaveMCInfo();

   // ---- PSA : same as the data -------------------------------------------
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
   hdb->SetJoinMethod(joinMethod.Data());
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

   std::cout << "\n(d,t) sim reco done -> " << outputFile << "   par = " << paramFile << std::endl;
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n" << std::endl;
}
