// seed = 0 keeps the original time-seeded behaviour. Pass a NON-ZERO seed when running
// several streams in parallel: the time-based seed has one-second granularity, so two
// streams launched in the same second draw byte-identical events and the extra CPU buys
// no extra statistics at all. An explicit seed also makes a run reproducible.
void He4He4_sim_el(Int_t nEvents = 10000, TString mcEngine = "TGeant4", UInt_t seed = 0)
{

   if (seed == 0) {
      srand((unsigned)time(NULL));
      seed = (float)rand() / RAND_MAX * 100000;
   }
   gRandom->SetSeed(seed);
   std::cout << " ==== Generator seed : " << seed << std::endl;


   TString dir = getenv("VMCWORKDIR");

   // Output file name
   TString outFile = "./data/attpcsim_in.root";

   // Parameter file name
   TString parFile = "./data/attpcpar_in.root";

   // -----   Timer   --------------------------------------------------------
   TStopwatch timer;
   timer.Start();
   // ------------------------------------------------------------------------

   // gSystem->Load("libAtGen.so");

   AtVertexPropagator *vertex_prop = new AtVertexPropagator();

   // -----   Create simulation run   ----------------------------------------
   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);      // Transport engine
   run->SetOutputFile(outFile); // Output file
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   // ------------------------------------------------------------------------

   // -----   Create media   -------------------------------------------------
   run->SetMaterials("media.geo"); // Materials (heco2_150 = He:CO2 90/10 at 150 torr)
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
   ATTPC->SetGeometryFileName("ATTPC_HeCO2_150torr.root");
   // ATTPC->SetModifyGeometry(kTRUE);
   run->AddModule(ATTPC);

   // ------------------------------------------------------------------------

   // -----   Magnetic field   -------------------------------------------
   // Constant Field
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., 5.691);                    // 0.5691 T, Dec 2014 (kG)
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // values are in cm
                                                          //  (xmin,xmax,ymin,ymax,zmin,zmax)
   run->SetField(fMagField);
   // --------------------------------------------------------------------

   // -----   Create PrimaryGenerator   --------------------------------------
   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // Beam Information
   Int_t z = 2;  // Atomic number
   Int_t a = 4;  // Mass number
   Int_t q = 0;  // Charge State
   Int_t m = 1;  // Multiplicity  NOTE: Due the limitation of the TGenPhaseSpace accepting only pointers/arrays the
                 // maximum multiplicity has been set to 10 particles.
   Double_t px = 0.000 / a; // X-Momentum / per nucleon!!!!!!
   Double_t py = 0.000 / a; // Y-Momentum / per nucleon!!!!!!
   Double_t pz = 0.241296 / a; // Z-Momentum / per nucleon. 4He at 1.95 MeV/u -> p = 241.296 MeV/c
   Double_t BExcEner = 0.0;
   Double_t Bmass = 4.00260325415;
   Double_t NomEnergy = 7.8; // MeV, only used for cross-section bookkeeping

   AtTPCIonGenerator *ionGen = new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy);
   ionGen->SetSpotRadius(0, -100, 0);
   // add the ion generator

   primGen->AddGenerator(ionGen);

   // primGen->SetBeam(1,1,0,0); //These parameters change the position of the vertex of every track
   // added to the Primary Generator
   // primGen->SetTarget(30,0);

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
   ResEner = 7.8; // MeV, 4He at 1.95 MeV/u

   // ---- Beam ----  4He, TRACKID=0
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(4.00260325415); // uma
   ExE.push_back(BExcEner);

   // ---- Target ----  4He in the He:CO2 gas
   Zp.push_back(2);
   Ap.push_back(4); //
   Qp.push_back(0); //
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(4.00260325415); // uma
   ExE.push_back(0.0);            // In MeV

   //--- Scattered -----
   Zp.push_back(2); //
   Ap.push_back(4); //
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(4.00260325415); // uma
   ExE.push_back(0.0);            // elastic: no excitation, the compound is unbound 8Be

   // ---- Recoil -----
   Zp.push_back(2); //
   Ap.push_back(4); //
   Qp.push_back(0); //
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(4.00260325415); // uma
   ExE.push_back(0.0);            // In MeV

   Double_t ThetaMinCMS = 2.0;  // wide open on purpose: the acceptance is the unknown
   Double_t ThetaMaxCMS = 180.0;

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
