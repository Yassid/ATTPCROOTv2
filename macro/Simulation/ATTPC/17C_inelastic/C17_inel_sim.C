/// @file C17_inel_sim.C
/// @brief 17C(p,p') and 17C(d,d') simulation for the M_n/M_p proposal (C17p_FRIB_Proposal.pdf,
///        this directory). One macro for both channels, selected by `channel`.
///
///   root -l 'C17_inel_sim.C(2000, "pp")'
///   root -l 'C17_inel_sim.C(2000, "dd", 0.217)'
///
/// WHY ONE MACRO FOR TWO CHANNELS. The proposal extracts M_n/M_p from the RATIO of the proton and
/// deuteron deformation lengths, so a systematic difference between the two simulations propagates
/// straight into the answer. Two sibling macros drift; one macro with a target switch cannot. The
/// only things `channel` changes are the target/ejectile mass and the transport gas -- beam, field,
/// drift, beam hole, angular range and generator settings are shared by construction.
///
/// It is otherwise C14_pp_sim.C with the beam and level scheme changed and nothing else, so that a
/// difference against the 14C(p,p') reference campaign (14C_pp/highfield/) is the beam and the
/// reaction rather than a setting.
///
/// CONDITIONS (from the proposal, ReA6):
///     beam    17C at 8.37 MeV/u = 142.29 MeV, 940 pps
///     gas     300 torr H2 for (p,p'), 300 torr D2 for (d,d'), 293 K -- one day of beam each
///     field   2.85 T is what Ref.[24] (Ayyad et al., Front. Phys. 13, 1539148 (2025)) ran at and
///             is the campaign's nominal; 4 T is SOLARIS's design field and is the second cell of
///             the campaign. See RESULTS.md -- the field turns out to be the only lever there is.
///     drift   1000 mm
///
/// LEVELS. 17C has S_n = 733 keV and exactly three bound states:
///     0      3/2+        ground state (elastic)
///     217    1/2+        proposal target
///     332    5/2+        proposal target  (ENSDF adopts 331; the proposal says 332)
/// The two excited states are 115 keV apart -- NOT the 130 keV the proposal text states -- and the
/// ground state is 217 keV below the 1/2+, so all three sit inside one 300 keV resolution width
/// with the elastic channel far stronger than either inelastic. Extracting the two inelastic
/// amplitudes from under the elastic is the measurement, and it is what this campaign quantifies.
///
/// ANGULAR DISTRIBUTION IS FLAT, ON PURPOSE. AtTPC2Body samples flat in cos(theta_cm)
/// (AtTPC2Body.cxx:169). The FRESCO distributions in fresco/ are applied as WEIGHTS in the
/// analysis, not here: with a flat generator the per-bin acceptance is a clean ratio and the Ex
/// resolution is measured per angular slice, so both are independent of the true shape and a
/// revised DWBA calculation never requires regenerating anything.
///
/// Beam momentum: the ion generator multiplies the pz argument by A, so pz is GeV/c per nucleon
/// and the TOTAL momentum is what sets the beam energy.
///     m(17C) = 17.02257865 u = 15.8564298 GeV
///     KE = 142.29 MeV  ->  p = 2.1290066 GeV/c  ->  pz = 2.1290066/17 per nucleon
/// This value is inherited from the (d,p) arm of the same proposal, where it was verified against
/// MC truth (17C_dp/RESULTS.md: momentum exact, KE 142.29 MeV = 8.370 MeV/u).
/// NomEnergy does NOT set the beam energy -- only pz does -- but it is not inert either: it is the
/// default for maxELoss, which controls how often a reaction is generated. maxELoss is passed
/// explicitly below.
///
/// @note Mass argument convention: every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. That discrepancy is
/// repo-wide and is kept here so this sim stays comparable to its siblings -- which is why the
/// beam MUST be verified from MC truth. Run check_beam_C17.C on the output; it is beam-agnostic.

void C17_inel_sim(Int_t nEvents = 2000, TString channel = "pp", Double_t resEx = 0.0, UInt_t seed = 0,
                  Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                  // CM ANGULAR RANGE. theta_cm is the PROJECTILE angle, so the light recoil comes
                  // out at theta_lab = (180 - theta_cm)/2 and the range below spans theta_lab
                  // 1-85 deg. It starts at 10 rather than 2 because below that the recoil carries
                  // under 220 keV and cannot make a track at all -- generating those events only
                  // dilutes the sample. Acceptance is a per-bin ratio, so a restricted generation
                  // range biases nothing as long as no bin outside it is quoted.
                  Double_t thetaMinCM = 10.0, Double_t thetaMaxCM = 178.0, TString mcEngine = "TGeant4")
{
   // ---- channel switch: the ONLY difference between (p,p') and (d,d') ----------------------
   Int_t tgtZ = 1, tgtA = 1;
   Double_t tgtMass = 1.0078250322; // amu
   TString geoFile = "ATTPC_H300torr_RT.root";
   if (channel == "dd") {
      tgtA = 2;
      tgtMass = 2.0141017778;
      // D2 at 300 torr, 293 K, rho = 6.5643e-5 g/cm3 -- matched EXACTLY to the digitisation par
      // (Density 0.065643 mg/cm3) so transport and digitisation model the same gas.
      geoFile = "ATTPC_D300torr_v2.root";
   } else if (channel != "pp") {
      std::cout << "\033[1;31m[C17_inel_sim] channel must be \"pp\" or \"dd\", got \"" << channel << "\"\033[0m"
                << std::endl;
      return;
   }
   std::cout << "\033[1;33m[C17_inel_sim] channel " << channel << " : target A = " << tgtA << ", gas " << geoFile
             << ", Ex(17C) = " << resEx << " MeV\033[0m" << std::endl;

   // RNG seed. seed = 0 keeps ROOT's default (time-based) behaviour; pass a distinct value per
   // parallel job, or byte-identical events get counted as added statistics.
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
   ATTPC->SetGeometryFileName(geoFile);
   run->AddModule(ATTPC);

   // -----   Magnetic field   -----------------------------------------------
   // SIGN: the a1954 DATA corresponds to -2.85 T, and AtSpyralPID::Direction() infers
   // forward/backward purely from the sense of rotation. Generating with the opposite sign made
   // its consistency check reject 82 % of sim tracks. Generate in the DATA convention so nothing
   // downstream needs a sign special-case; the 4 T cell of this campaign passes -40.0 for the
   // same reason.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 17C at 142.29 MeV (8.37 MeV/u)   ------------------------
   Int_t z = 6;  // Atomic number
   Int_t a = 17; // Mass number
   Int_t q = 0;  // Charge state
   Int_t m = 1;  // Multiplicity
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = 2.1290066 / a; // GeV/c per nucleon -> 142.29 MeV total KE
   Double_t BExcEner = 0.0;
   Double_t Bmass = 17.02257865; // amu (repo convention -- see the note above)
   Double_t NomEnergy = 142.29;  // MeV

   // maxELoss CONTROLS THE REACTION RATE -- it is not cosmetic, and it also controls whether the
   // VERTEX DISTRIBUTION IS UNIFORM, which every resolution number in this campaign depends on.
   // AtTPCIonGenerator does  fMaxEnLoss = (eLoss < 0 ? ener : eLoss)  and then samples
   //     Er = gRandom->Uniform(0, fMaxEnLoss);  SetRndELoss(Er)
   // while AtTpc::reactionOccursHere() fires only once  fELossAcc > GetRndELoss().
   // Set BELOW the traversal loss, every event reacts but the vertices pile up upstream; set far
   // above it, most events never react at all. Set at the traversal loss, essentially every event
   // fires with the vertex spread across the full metre.
   // The (d,p) arm measured 14.5 MeV of loss over the metre in D2 from its own truth, and H2 at
   // 300 torr has the same electron density to 0.8 %, so the same number holds for both channels.
   // 15.0 errs marginally high, which costs a few per cent of events and protects the uniformity.
   Double_t maxELoss = 15.0; // MeV, ~ the 17C energy loss across the full drift length

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0);
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 17C(p,p')17C  /  17C(d,d')17C   ---------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, scattered, recoil -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 17C ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(Bmass);
   ExE.push_back(BExcEner);

   // ---- Target : p or d ----
   Zp.push_back(tgtZ);
   Ap.push_back(tgtA);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(tgtMass);
   ExE.push_back(0.0);

   // ---- Scattered : 17C, left at resEx ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(Bmass);
   // Excitation of the RESIDUAL 17C. 0 = elastic; 0.217 = 1/2+; 0.332 = 5/2+. Inelastic changes
   // the two-body kinematics, so theta_lab(theta_cm) and the recoil energy both shift -- which is
   // exactly why acceptance has to be measured per level rather than assumed from the elastic.
   // All three levels are BOUND (S_n = 733 keV), so the residual does not break up and the heavy
   // track in the chamber really is 17C.
   ExE.push_back(resEx);

   // ---- Recoil : p or d (this is the track that gets fitted) ----
   Zp.push_back(tgtZ);
   Ap.push_back(tgtA);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(tgtMass);
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg"
             << "   (recoil theta_lab ~ " << 0.5 * (180 - ThetaMaxCMS) << " - " << 0.5 * (180 - ThetaMinCMS) << " deg)"
             << std::endl;

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
   run->CreateGeometryFile("./data/geofile_C17_inel_full.root");

   timer.Stop();
   std::cout << std::endl << std::endl;
   std::cout << "Macro finished succesfully." << std::endl;
   std::cout << "Output file is " << outFile << std::endl;
   std::cout << "Parameter file is " << parFile << std::endl;
   std::cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl
             << std::endl;
}
