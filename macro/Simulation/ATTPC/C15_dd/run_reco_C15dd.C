/// @file run_reco_C15dd.C
/// @brief Digitise + reconstruct the a2091 15C(d,d') simulation with the SAME chain the C15d DATA
///        uses, so the acceptance it yields corrects the production actually being corrected.
///        Adapted from 16C_dt/run_reco_C16dt.C (a1975).
///
/// THE CHAIN IS COPIED FROM THE DATA, not from the a1975 simulation. `rereco_mfi_C15d.sh` runs
///   unpackReco_C15d.C(..., "ATTPC.C15d_a2091_D2.par", "ATTPC_D300torr_v2_geomanager.root",
///                     psaThreshold 20.0, clusterRadius 15.0, clusterDistance 7.5,
///                     doPID true, gainMatch false, "", pidMinPoints 30, pidZTieTol 0.0,
///                     bField 2.85, 0, "multifit", "tc", doSC false, floatTau false, pidInt true)
/// so this macro uses AtPSAMultiFit + AtTrackFinderTC, NOT HDBSCAN, and runs NO cleaning task --
/// the data's PRA reads AtEventH directly (`praTask->SetInputBranch(doSC ? ... : "AtEventH")`).
/// Using the a1975 chain here would measure a different reconstruction's losses.
///
/// *** GAIN IS THE OPEN QUESTION, AND IT IS NOT COSMETIC ***  (inherited verbatim from a1975)
///     ATTPC.C15d_a2091_D2.par declares Gain = 10000. Gain sets how many electrons land on a pad,
///     hence whether the pad crosses the PSA threshold, hence how many hits a track keeps -- which
///     is precisely what an ACCEPTANCE measures. Too low a gain silently deletes short and low
///     dE/dx tracks and the acceptance comes out too small, with nothing looking wrong.
///     DO NOT ADOPT AN ACCEPTANCE FROM THIS UNTIL THE GAIN IS VALIDATED against the data:
///     compare hits-per-track and pads-per-track against /home/yassid/Data/C15d_reco_mfi.
///     This matters more here than for a1975: the whole point is the LOW-ENERGY deuterons
///     (theta_CM bin 0, KE ~ 1.4 MeV), which are exactly the tracks a low gain deletes.
///
/// ⚠ FIELD SIGN. The simulation generates at -28.5 kG, the a1975 data convention; the C15d fit is
///     given bField = +2.85 T. If the sense of rotation disagrees, anything inferring direction
///     from curvature rejects most simulated tracks and the acceptance measures the sign error
///     instead of the detector. VERIFY before trusting any number from this.
///
/// ⚠⚠ DO NOT COMBINE doPID = kTRUE WITH nEvents > 0.  `fRun->Run(0, nEvents)` together with
///     AtPIDTask writes an EMPTY AtPatternEvent -- zero tracks, silently, with the job exiting 0
///     and a well-formed file.  Measured on the same 200 events of ex00.00:
///         nEvents = 200, doPID = kFALSE  ->  96 tracks
///         nEvents = 200, doPID = kTRUE   ->   0 tracks
///         nEvents =   0, doPID = kTRUE   ->  96 tracks   (the production; 97 % valid PID)
///     So the full production is correct and only TRUNCATED test runs are affected. For a quick
///     test either pass doPID = kFALSE or run the whole file.
///
///   root -l -b -q 'run_reco_C15dd.C("/home/yassid/Data/C15dd_sim/ex00.00_sim.root","/home/yassid/Data/C15dd_sim/ex00.00_reco.root")'

#include "../AtBeamHole.h"

void run_reco_C15dd(TString mcFile = "./data/C15dd_sim.root", TString outputFile = "./data/C15dd_reco.root",
                    TString paramFile = "ATTPC.C15d_sim_g150000.par", Double_t thr = 20,
                    Double_t clusterRadius = 15.0, Double_t clusterDistance = 7.5,
                    Int_t nEvents = 0, Double_t holeR = 20.0, Double_t adcCeiling = 0,
                    Bool_t doPID = kTRUE, Double_t bField = 2.85, Int_t pidMinPoints = 30,
                    Double_t pidZTieTol = 0.0, Bool_t pidUseIntegral = kTRUE)
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
   // ⚠ BOTH this and psaTask MUST stay kTRUE.  The reco files are ~7 GB per 20k-event point
   // because AtRawEvent and AtEventH are persisted, and it is tempting to switch them off since
   // nothing downstream reads them back -- but the chain passes data through the persisted
   // branches, and with either one kFALSE the PRA produces ZERO tracks.  Silently: the job exits
   // 0 and writes a well-formed file with an empty AtPatternEvent.  Tested both ways.
   pulse->SetPersistence(kTRUE);
   // MC truth (A, Z, trackID) on every hit: without it the acceptance cannot tell the triton from
   // the beam or the recoil, and every downstream truth match becomes a guess.
   pulse->SetSaveMCInfo();

   // ---- PSA : same as the data -------------------------------------------
   auto p = new AtPSAMultiFit();
   // EXACTLY what the data sets, and nothing else.  unpackReco_C15d.C touches only the threshold
   // and the floated peaking time; everything else is the class default.  The a1975 macro this was
   // copied from also set SetPeakingTime(0.720), SetMaxPeaks(4), SetMinSeparation(4) and
   // SetProminenceOnPrimary(false) -- keeping those here would reconstruct the simulation with a
   // DIFFERENT pulse analysis from the data and the acceptance would correct the wrong chain.
   p->SetThreshold(thr);
   // pidUseIntegral = true in the data production, which implies a floated peaking time: with tau
   // fixed the trace integral is exactly 2.9618x the amplitude and carries no extra information
   // (see the AtPSAMultiFit trace-integral fix, commit 9eef5aad).
   p->SetFloatPeakingTime(kTRUE);
   AtPSAtask *psaTask = new AtPSAtask(p);
   // ⚠ MUST stay kTRUE.  With it kFALSE the PRA reads an empty AtEventH and produces ZERO
   // tracks -- silently, with the job exiting 0 and a well-formed output file.  Tested.
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   // ---- pattern recognition : AtTrackFinderTC, as the C15d data ------------
   // NO cleaning task: the data's PRA reads AtEventH directly.
   auto tc = std::make_unique<AtPATTERN::AtTrackFinderTC>();
   tc->SetClusterRadius(clusterRadius);
   tc->SetClusterDistance(clusterDistance);

   AtPRAtask *praTask = new AtPRAtask(std::move(tc));
   praTask->SetInputBranch("AtEventH");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);

   // ---- PID : same as the data -------------------------------------------
   // AtGenfitter gates on AtSpyralResult::valid, so without AtPIDTask the fit stage has nothing
   // to select on and the acceptance would stop at the pattern-recognition level.  The data's
   // kinematics tree (dd/kin_dd_C15d.C) takes ke and theta from the FITTED sample, so the
   // response function needs this branch to exist.
   AtPIDTask *pidTask = nullptr;
   if (doPID) {
      pidTask = new AtPIDTask();
      pidTask->SetInputBranch("AtPatternEvent");
      pidTask->SetOutputBranch("AtPIDEvent");
      pidTask->SetBField(bField);
      pidTask->SetMinPoints(pidMinPoints);
      pidTask->SetZTieTolerance(pidZTieTol);
      if (pidUseIntegral)
         pidTask->SetUseTraceIntegral(kTRUE);
      pidTask->SetPersistence(kTRUE);
   }

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);
   if (pidTask)
      fRun->AddTask(pidTask);

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
