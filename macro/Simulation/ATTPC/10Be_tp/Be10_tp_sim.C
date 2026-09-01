/// @file Be10_tp_sim.C
/// @brief 10Be(t,p)12Be simulation -- the TWO-NEUTRON-TRANSFER counterpart of the 14C(d,p) study.
///
/// WHY THIS EXISTS. The 14C(d,p)15C campaign (macro/Simulation/ATTPC/14C_dp/) established, over a
/// 3 field x 2 pad-pitch matrix, that the backward proton from a one-neutron transfer is already
/// measured to sigma(Ex) ~ 0.2 MeV by the AT-TPC as it exists, and that all the gain from a higher
/// field or finer pads is in the FORWARD hemisphere. 10Be(t,p)12Be asks the same question of a
/// two-neutron transfer whose interesting physics is a set of levels 100-600 keV apart:
///
///     0+  g.s.
///     2+  2.109
///     0+  2.251     <- the intruder 0+_2; only 142 keV above the 2+, and weakly populated
///     1-  2.715
///
/// S_n(12Be) = 3.171 MeV, so ALL FOUR ARE BOUND -- unlike 15C, where only two levels were.
/// The 2.109/2.251 pair is the hard case: a 142 keV separation against a resolution the (d,p)
/// study measured at 0.20-0.71 MeV depending on hemisphere and configuration.
///
/// CONDITIONS, deliberately identical to the (d,p) campaign except for the target gas and the
/// reaction:
///     B = 2.85 / 4 / 7 T,  300 torr,  drift length 1000 mm, same vessel and window
///     pad planes: the real AT-TPC 8x12 mm, and a square 2 mm plane
///     beam 10Be at 115 MeV = 11.5 MeV/u -- the SAME ENERGY PER NUCLEON as a1954's 14C at 161 MeV,
///     so the beam velocity, and hence the CM system and the ejectile kinematics, sit on the same
///     footing as the channel this is being compared against.
///
/// Q VALUE. Q(t,p) = S_2n(12Be) - S_2n(3H) = 3.673 - 8.482 = -4.81 MeV. Negative but far below
/// the 115 MeV available, so every level here is open at every angle.
///
/// GAS. ATTPC_T300torr.root, medium TargetT2_300, rho = 9.89854e-5 g/cm3 -- T2 at 300 torr and
/// 293.15 K, the ideal-gas value, and the digitisation par (make_tp_pars.sh) is set to exactly the
/// same number so transport and digitisation model the same gas. That is the a1954 lesson: a 7.4 %
/// mismatch between the two was a real defect once material effects went on.
/// T2 at 300 torr has the SAME molecular number density, hence the same ELECTRON density, as the
/// H2 and D2 fills of the other campaigns, so the beam energy loss per metre is the H2 one scaled
/// by Z_beam^2 alone: 14C lost ~12 MeV/m at 11.5 MeV/u, so 10Be loses ~12 x (16/36) = 5.3 MeV/m.
/// That is what maxELoss below is set from.
///
/// Beam momentum: the ion generator multiplies the pz argument by A, so pz is GeV/c per nucleon
/// and the TOTAL momentum sets the beam energy.
///     m(10Be) = 10.0135341 u = 9.3275465 GeV
///     KE = 115 MeV  ->  p = sqrt(115*(2m+115)) = 1.469204 GeV/c  ->  pz = 1.469204/10 per nucleon
/// NomEnergy does NOT set the beam energy -- only pz does. Pass maxELoss explicitly; it is the
/// argument that controls how often a reaction is generated at all (see the comment below).
///
///   root -l 'Be10_tp_sim.C(2000)'
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. That discrepancy is
/// repo-wide and is kept here so this sim stays comparable to its siblings -- which means the beam
/// energy must be VERIFIED from MC truth rather than assumed. check_vertex_beam_Be10.C does that,
/// and it also returns the mean beam energy AT THE REACTION VERTEX, which is the number the
/// two-body inversion in the acceptance and resolution macros needs.

/// CM ANGULAR RANGE: 2-178 by default, i.e. everything. As in (d,p), the proton does NOT follow
/// theta_lab = (180-theta_cm)/2; the region where a transfer angular distribution carries its
/// yield is small theta_cm, which maps to the BACKWARD lab hemisphere. Do not preselect.
void Be10_tp_sim(Int_t nEvents = 2000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                 TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                 Double_t resEx = 0.0, UInt_t seed = 0, TString geoFile = "ATTPC_T300torr.root")
{
   // RNG seed. seed = 0 keeps ROOT's default (time-based) behaviour; pass a distinct value per
   // parallel job or the jobs produce byte-identical events and the "added statistics" is a copy.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed. gRandom->GetSeed() returns TRandom3's internal state counter
   // (always 624 right after SetSeed), not the seed -- useless for verifying a parallel run.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar.root";

   TStopwatch timer;
   timer.Start();

   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);
   run->SetOutputFile(outFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   run->SetMaterials("media.geo");

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName(geoFile); // default ATTPC_T300torr.root, rho = 9.89854e-5 g/cm3
   std::cout << "\033[1;33m[Be10_tp_sim] transport gas geometry: " << geoFile << "\033[0m" << std::endl;
   run->AddModule(ATTPC);

   // -----   Magnetic field   -----------------------------------------------
   // SIGN: generate in the a1954 DATA convention (negative), as the (d,p) and (p,p') campaigns do.
   // AtSpyralPID::Direction() infers forward/backward from the sense of rotation alone, so
   // generating with the opposite sign makes its consistency check reject most of the sample.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 10Be at 115 MeV (11.5 MeV/u)   --------------------------
   Int_t z = 4;  // Atomic number
   Int_t a = 10; // Mass number
   Int_t q = 0;  // Charge state
   Int_t m = 1;  // Multiplicity
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = 1.469204 / a; // GeV/c per nucleon -> 115.0 MeV total KE
   Double_t BExcEner = 0.0;
   Double_t Bmass = 10.0135341; // amu (repo convention -- see the note above)
   Double_t NomEnergy = 115.0;  // MeV

   // maxELoss CONTROLS THE REACTION RATE -- it is not cosmetic.
   // AtTPCIonGenerator does  fMaxEnLoss = (eLoss < 0 ? ener : eLoss)  and then samples
   //     Er = gRandom->Uniform(0, fMaxEnLoss);  SetRndELoss(Er)
   // while AtTpc::reactionOccursHere() fires only once  fELossAcc > GetRndELoss().
   // Leaving it at the -1 default makes the threshold uniform in [0, 115], so it exceeds the
   // available loss ~95 % of the time and almost no reaction is generated. Set it to the actual
   // traversal loss and essentially every reaction event fires, with the vertex still uniform.
   // 10Be at 11.5 MeV/u in T2 at 300 torr: same electron density as the H2/D2 fills, so the
   // 14C figure of ~12 MeV/m scales by Z^2 -> 12 * 16/36 = 5.3 MeV over the metre.
   Double_t maxELoss = 5.5; // MeV

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 10Be(t,p)12Be   ------------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, scattered, recoil -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 10Be ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(10.0135341); // uma
   ExE.push_back(BExcEner);

   // ---- Target : t ----
   Zp.push_back(1);
   Ap.push_back(3);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(3.0160493); // uma
   ExE.push_back(0.0);

   // ---- Scattered : 12Be (the residual, and the one that carries the excitation) ----
   Zp.push_back(4);
   Ap.push_back(12);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(12.0269221); // uma, 12Be ground state (mass excess 25.077 MeV)
   // Excitation of the RESIDUAL 12Be. 0 = 0+ g.s.; 2.109 = 2+; 2.251 = the intruder 0+_2;
   // 2.715 = 1-. A transfer to an excited state changes the two-body kinematics, so
   // theta_lab(theta_cm) and the proton energy both shift -- which is why acceptance has to be
   // measured per level rather than assumed from the ground state.
   ExE.push_back(resEx);

   // ---- Recoil : p ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.0078250); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg" << std::endl;

   AtTPC2Body *TwoBody = new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE, ResEner,
                                        ThetaMinCMS, ThetaMaxCMS);
   primGen->AddGenerator(TwoBody);

   run->SetGenerator(primGen);
   run->SetStoreTraj(kTRUE);

   run->Init();

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);
   run->CreateGeometryFile("./data/geofile_Be10_tp_full.root");

   timer.Stop();
   cout << endl << endl;
   cout << "Macro finished succesfully." << endl;
   cout << "Output file is " << outFile << endl;
   cout << "Parameter file is " << parFile << endl;
   cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << endl << endl;
}
