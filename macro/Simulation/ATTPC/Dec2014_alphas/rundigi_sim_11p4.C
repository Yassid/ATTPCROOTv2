// Digitisation of the simulated Dec 2014 alpha+alpha elastic events.
//
// Chain: AtClusterizeTask -> AtPulseTask -> AtPSAtask, i.e. the same PSA the experimental
// data goes through, so simulated and measured events are reconstructed identically.
//
// The point of running this is not only acceptance. AtClusterizeTask now applies the
// Langevin E x B lateral drift in the forward direction, and AtPSA removes it in the
// reverse direction. Since the MC truth position of every ionisation point is known, the
// residual between truth and the reconstructed AtHit is a direct measurement of how well
// the Lorentz/Langevin correction works -- something no amount of experimental data can
// give, because there the truth is never known.
//
// Parameters are the calibrated ones from ATTPC.alpha_150torr_sim.par (v_D = 2.251 cm/us,
// tilt 6.47 deg, ThetaRot -161.9, ThetaPad 110.9, E = 12 kV/m, B = 0.5691 T).
//
//   root -l 'rundigi_sim.C("./data11p4/attpcsim_in.root")'
void rundigi_sim_11p4(TString mcFile = "./data11p4/attpcsim_in.root",
                 TString outFile = "./data11p4/output_digi.root",
                 Bool_t keepElectrons = kTRUE, Bool_t keepTraces = kTRUE,
                 Int_t psaThreshold = 3)
{
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString scriptdir = dir + "/scripts/Lookup20141208.xml";   // the Dec 2014 pad map
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());

   FairRunAna *fRun = new FairRunAna();
   // FairRoot 18.6 dropped FairRunAna::SetInputFile -- input comes from a FairSource now.
   fRun->SetSource(new FairFileSource(mcFile));
   fRun->SetGeomFile(geomDir + "ATTPC_HeCO2_150torr_geomanager.root");
   fRun->SetOutputFile(outFile);

   TString digiParFile = dir + "/parameters/ATTPC.Dec2014.P300_B0_sim.par";
   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   // ---- ionisation, diffusion and the Langevin drift -----------------------
   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   // The drifted-electron collection is ~90% of the output volume and is only needed for
   // the truth-residual study. Bulk production for charge-spectrum work turns it off.
   clusterizer->SetPersistence(keepElectrons);

   // ---- pad response and electronics --------------------------------------
   AtPulseTask *pulse = new AtPulseTask();
   // The pad TRACES are the other bulk term.  AtPSAtask consumes them in memory inside the
   // same run, so they only need to reach disk if you intend to look at waveforms later.
   // With keepElectrons and keepTraces both off, a 6000-event digi.root goes from 50 GB to
   // a few tens of MB -- which is what makes high-statistics efficiency runs possible.
   pulse->SetPersistence(keepTraces);
   pulse->SetSaveMCInfo();               // keeps the MC point each pad signal came from

   // ---- exactly the PSA used on the experimental data ---------------------
   AtPSASimple2 *psa = new AtPSASimple2();
   psa->SetThreshold(psaThreshold);   // NOT the 20 used on real data -- the simulated ADC scale is
                           // different, and the working sims on develop use 3 (or 0).
                           // With 20 every simulated event comes out with zero hits.
   psa->SetMaxFinder();
   psa->SetBaseCorrection(kTRUE);
   psa->SetTimeCorrection(kFALSE);

   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(kTRUE);

   // NB: no AtPRAtask here. Pattern recognition is one of the PCL-only classes and is not
   // built on this branch; clustering is done offline instead (calib/cluster_tracks.py).

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);

   fRun->Init();
   fRun->Run(0, 0);

   timer.Stop();
   std::cout << "Digitisation done. Real time " << timer.RealTime()
             << " s, CPU " << timer.CpuTime() << " s" << std::endl;
}

/* STATUS -- the chain runs end to end; the REACTION RATE is the open problem.
 *
 *   Geant4 -> AtClusterizeTask -> AtPulseTask -> AtPSAtask -> AtEvent -> tracker
 *
 * all works, and produces the same AtEvent class the real data does, so the tracker runs
 * on both identically.  Settings that had to be right before anything appeared at all:
 *
 *   PSA threshold 3, NOT the 20 used on real data.  The simulated ADC scale is different;
 *   the working sims on develop use 3 or 0.  At 20 every event comes out with zero hits
 *   and nothing reports an error.
 *
 *   Gain 20000 and CoefL/CoefT 0.00002/0.0002, from ATTPC.e20020_sim.par.  Taking these
 *   from the EXPERIMENTAL par file instead gives Gain 0.28 -- a different quantity, the
 *   measured pulses being already amplified -- and pulses peak at 0.03 ADC.
 *
 * REMAINING PROBLEM: the two-body reaction fires in only 32% of attempts, and when it
 * does, only at shallow vertices.  Measured over 200 events:
 *
 *     reaction attempts                  200
 *     vertex energy <= 0, no reaction    137  (68%)
 *     vertex energy when it fires        7.0 - 10.6 MeV, median 8.2
 *
 * So the simulated sample covers E_cm 3.5-5.3 MeV only -- the low-energy half, where the
 * 8Be 2+ sits, is never populated.  AtTPC2Body.cxx:225 requires a non-zero vertex energy
 * from AtVertexPropagator and its comment says "not punch through".  This beam HAS a
 * 1071 mm range in a 1000 mm chamber, so it punches through with ~2.1 MeV in hand -- which
 * is physically correct and is the same fact that sets the E_cm >= 1.07 MeV floor in the
 * data.  The propagator appears to treat a beam that does not stop inside as having no
 * valid vertex.  That is what to fix before any efficiency is taken from this.
 */
