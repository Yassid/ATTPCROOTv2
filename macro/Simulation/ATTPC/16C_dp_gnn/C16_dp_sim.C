// 16C(d,p)17C simulation for GNN training data (labeled point clouds).
// Adapted from 15C_d/C15_dp_sim.C: beam 16C, products 17C + p, D2 gas, B = 2.85 T
// (matches a1975 deuterium runs). Digitize with SaveMCInfo to keep per-pad MC track IDs.
// FLAG: thetaCMSmin/max select the CM angle range of the (d,p) reaction. Current focus
// = BACKWARD-LAB protons (lab 90-180). In this inverse-kinematics reaction FORWARD CM
// maps to BACKWARD lab (empirical: CM[2,60]->lab 76-168, 61% backward), so the default
// band is CM[2,55]. The definitive lab>=90 cut is applied in package_training.py.
// To revert to ALL angles later: call C16_dp_sim(N,"TGeant4", 2, 178).
void C16_dp_sim(Int_t nEvents = 5000, TString mcEngine = "TGeant4",
                Double_t thetaCMSmin = 2.0, Double_t thetaCMSmax = 55.0, UInt_t seedArg = 0)
{
   // seedArg>0 => deterministic distinct seed (for parallel chunks); 0 => time-based
   UInt_t seed = seedArg;
   if (seed == 0) { srand((unsigned)time(NULL)); seed = (float)rand() / RAND_MAX * 100000; }
   gRandom->SetSeed(seed);

   TString dir = getenv("VMCWORKDIR");

   // Output file name
   TString outFile = "./data/attpcsim.root";

   // Parameter file name
   TString parFile = "./data/attpcpar.root";

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   // ------------------------------------------------------------------------

   // gSystem->Load("libAtGen.so");

   // -----   Create simulation run   ----------------------------------------
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      // Transport engine
   run->SetOutputFile(outFile); // Output file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   // ------------------------------------------------------------------------

   // -----   Create media   -------------------------------------------------
   run->SetMaterials("media.geo"); // Materials
   // ------------------------------------------------------------------------

   // -----   Create geometry   ----------------------------------------------

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   // FairModule* magnet = new AtMagnet("Magnet");
   // run->AddModule(magnet);

   /*FairModule* pipe = new AtPipe("Pipe");
   run->AddModule(pipe);*/

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("ATTPC_D300torr_v2.root");
   // ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(ATTPC);

   // ------------------------------------------------------------------------

   // -----   Magnetic field   -------------------------------------------
   // Constant Field
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 28.5);                     // values are in kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // values are in cm
                                                          //  (xmin,xmax,ymin,ymax,zmin,zmax)
   run->SetField(fMagField);
   // --------------------------------------------------------------------

   // -----   Create PrimaryGenerator   --------------------------------------
   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // Beam Information  (16C @ ~12 MeV/u to match a1975; total p ~2.43 GeV/c)
   Int_t z = 6;  // Atomic number
   Int_t a = 16; // Mass number
   Int_t q = 0;  // Charge State
   Int_t m = 1;  // Multiplicity  NOTE: Due the limitation of the TGenPhaseSpace accepting only pointers/arrays the
                 // maximum multiplicity has been set to 10 particles.
   Double_t px = 0.000 / a; // X-Momentum / per nucleon!!!!!!
   Double_t py = 0.000 / a; // Y-Momentum / per nucleon!!!!!!
   Double_t pz = 2.432 / a; // Z-Momentum / per nucleon!!!!!!
   Double_t BExcEner = 0.0;
   Double_t Bmass = 16.0147;
   Double_t NomEnergy = 11.77;

   AtTPCIonGenerator *ionGen = new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy);
   ionGen->SetSpotRadius(0, -100, 0);
   // add the ion generator

   primGen->AddGenerator(ionGen);

   // primGen->SetBeam(1,1,0,0); //These parameters change the position of the vertex of every track added to the
   // Primary Generator
   //  primGen->SetTarget(30,0);

   // Variables for 2-Body kinematics reaction
   std::vector<Int_t> Zp;      // Zp
   std::vector<Int_t> Ap;      // Ap
   std::vector<Int_t> Qp;      // Electric charge
   Int_t mult;                 // Number of particles
   std::vector<Double_t> Pxp;  // Px momentum X
   std::vector<Double_t> Pyp;  // Py momentum Y
   std::vector<Double_t> Pzp;  // Pz momentum Z
   std::vector<Double_t> Mass; // Masses
   std::vector<Double_t> ExE;  // Excitation energy
   Double_t ResEner;           // Energy of the beam (Useless for the moment)

   // Note: Momentum will be calculated from the phase Space according to the residual energy of the beam

   mult = 4; // Number of Nuclei involved in the reaction (Should be always 4) THIS DEFINITION IS MANDATORY (and the
             // number of particles must be the same)
   ResEner = 11.77; // MeV/u (vertex)

   // ---- Beam ----  16C  TRACKID=0
   Zp.push_back(z);
   Ap.push_back(a); //
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(16.0147); // uma
   ExE.push_back(BExcEner);

   // ---- Target ----
   Zp.push_back(1); // d
   Ap.push_back(2); //
   Qp.push_back(0); //
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.01410177785); // uma
   ExE.push_back(0.0);            // In MeV

   //--- Scattered (heavy recoil 17C) -----  TRACKID=2
   Zp.push_back(6);  //
   Ap.push_back(17); //
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(17.0226); // uma
   ExE.push_back(0.0);      //

   // ---- Recoil (light ejectile proton) -----  TRACKID=3
   Zp.push_back(1); // p
   Ap.push_back(1); //
   Qp.push_back(0); //
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.00727646); // uma
   ExE.push_back(0.0);         // In MeV

   // CM range (FLAG: set via macro args; default = backward-lab focused band)
   Double_t ThetaMinCMS = thetaCMSmin;
   Double_t ThetaMaxCMS = thetaCMSmax;

   AtTPC2Body *TwoBody =
      new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE, ResEner, ThetaMinCMS, ThetaMaxCMS);
   primGen->AddGenerator(TwoBody);

   run->SetGenerator(primGen);

   // ------------------------------------------------------------------------

   //---Store the visualiztion info of the tracks, this make the output file very large!!
   //--- Use it only to display but not for production!
   run->SetStoreTraj(kTRUE);

   // -----   Initialize simulation run   ------------------------------------
   run->Init();
   // ------------------------------------------------------------------------

   // -----   Runtime database   ---------------------------------------------

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();
   // ------------------------------------------------------------------------

   // -----   Start run   ----------------------------------------------------
   run->Run(nEvents);

   // You can export your ROOT geometry ot a separate file
   run->CreateGeometryFile("./data/geofile_full.root");
   // ------------------------------------------------------------------------

   // -----   Finish   -------------------------------------------------------
   timer.Stop();
   Double_t rtime = timer.RealTime();
   Double_t ctime = timer.CpuTime();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << rtime << " s, CPU time " << ctime << "s" << endl << endl;
   // ------------------------------------------------------------------------
}
