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
//   root -l 'rundigi_sim.C("./data/attpcsim_in.root")'
void rundigi_sim(TString mcFile = "./data/attpcsim_in.root",
                 TString outFile = "./data/output_digi.root",
                 Bool_t keepElectrons = kTRUE)
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

   TString digiParFile = dir + "/parameters/ATTPC.alpha_150torr_sim.par";
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
   pulse->SetPersistence(kTRUE);
   pulse->SetSaveMCInfo();               // keeps the MC point each pad signal came from

   // ---- exactly the PSA used on the experimental data ---------------------
   AtPSASimple2 *psa = new AtPSASimple2();
   psa->SetThreshold(20);                // same as run_unpack_Dec2014_alphas.C
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
