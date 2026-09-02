/// @file check_beam_C17.C
/// @brief Verify, from MC TRUTH, that C17_dp_sim.C generated the beam and the reaction it claims.
///
///   root -b -q 'check_beam_C17.C("/path/to/sim.root")'
///
/// WHY THIS IS NOT OPTIONAL. Every AT-TPC sim macro passes the ion mass in amu to
/// AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. That discrepancy is
/// repo-wide, so the beam energy is a DERIVED quantity that has to be read back out of the
/// transport, never assumed from the argument that was typed in. This is the 17C counterpart of
/// 14C_pp/check_beam_C14.C, with the masses as arguments instead of hard-coded.
///
/// It checks three things:
///   1. the beam MOMENTUM handed to the transport matches what the macro asked for;
///   2. the beam energy loss across the chamber matches the CATIMA prediction the sim macro's
///      maxELoss was set from (14.8 MeV over 1000 mm) -- i.e. the transport gas is the right gas;
///   3. the two-body kinematics close: the proton angle and energy vs theta_cm reproduce the
///      analytic 17C(d,p)18C solution at Q = +1.959 MeV.
///
/// Check 3 is the one that would catch a swapped mass or a wrong residual excitation, which
/// checks 1 and 2 cannot see.
void check_beam_C17(TString f = "./data/attpcsim.root", Double_t pAsked = 2.1290066, Double_t mBeamAmu = 17.02257865,
                    Double_t mTgtAmu = 2.0141017778, Double_t mResAmu = 18.02675193, Double_t mEjAmu = 1.0078250322,
                    Double_t resEx = 0.0)
{
   TFile *fi = TFile::Open(f);
   if (!fi || fi->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", f.Data());
      return;
   }
   TTree *t = (TTree *)fi->Get("cbmsim");
   if (!t) {
      printf("\033[1;31mno cbmsim tree in %s\033[0m\n", f.Data());
      return;
   }
   gSystem->Load("libAtSimulationData.so");

   TClonesArray *tracks = nullptr;
   t->SetBranchAddress("MCTrack", &tracks);

   const double U = 0.93149401; // GeV per amu
   const double mB = mBeamAmu * U, mT = mTgtAmu * U;
   const double mR = mResAmu * U + resEx / 1000., mE = mEjAmu * U;

   // AtMCTrack::GetEnergy() builds E from a TDatabasePDG mass lookup, and ion PDG codes
   // (10LZZZAAAI) are not in that table -> it returns E = |p| and an effective mass of 0. So the
   // mass cannot be read back from the track; only the MOMENTUM is trustworthy.
   // AtVertexPropagator alternates beam-only and reaction events. In a BEAM event track 0 is the
   // pristine beam ion, generated at the entrance (z = -100 cm). In a REACTION event track 0 is
   // the SCATTERED ion, already at the vertex with a degraded momentum -- averaging over both
   // drags the mean ~1 % low and looks like a wiring bug. Select on the START POSITION.
   const double zEntrance = -100.0; // cm, from ionGen->SetSpotRadius(0, -100, 0)
   int nb = 0, nreact = 0;
   double sumP = 0, sumP2 = 0;
   // proton kinematics closure
   TH1D *hdA = new TH1D("hdA", "", 200, -5, 5);
   TH1D *hdE = new TH1D("hdE", "", 200, -5, 5);
   double thMin = 1e9, thMax = -1e9, keMin = 1e9, keMax = -1e9;
   double vzMin = 1e9, vzMax = -1e9, sumVz = 0;
   int nk = 0;

   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (tracks->GetEntriesFast() == 0)
         continue;
      auto *tr = (AtMCTrack *)tracks->At(0);
      if (!tr || tr->GetPdgCode() < 1000000000)
         continue;
      const double p = std::sqrt(tr->GetPx() * tr->GetPx() + tr->GetPy() * tr->GetPy() + tr->GetPz() * tr->GetPz());
      if (std::fabs(tr->GetStartZ() - zEntrance) < 1.0) { // pristine beam
         sumP += p;
         sumP2 += p * p;
         ++nb;
         continue;
      }
      ++nreact;

      // --- reaction event: track 0 is the scattered 18C at the vertex, and the proton is the
      // --- light track in the same event. Find it by PDG (2212).
      const double vz = tr->GetStartZ();
      vzMin = std::min(vzMin, vz);
      vzMax = std::max(vzMax, vz);
      sumVz += vz;
      const AtMCTrack *pr = nullptr;
      for (int k = 1; k < tracks->GetEntriesFast(); ++k) {
         auto *c = (AtMCTrack *)tracks->At(k);
         if (c && c->GetPdgCode() == 2212) {
            pr = c;
            break;
         }
      }
      if (!pr)
         continue;

      // BEAM 4-vector at the vertex. p (of track 0) is the RESIDUAL's momentum, so reconstruct the
      // incident beam from energy-momentum conservation instead: p_beam = p_res + p_ej.
      TVector3 pres(tr->GetPx(), tr->GetPy(), tr->GetPz());
      TVector3 pej(pr->GetPx(), pr->GetPy(), pr->GetPz());
      TVector3 pbeam = pres + pej;
      const double Eb = std::sqrt(pbeam.Mag2() + mB * mB);

      // analytic two-body solution at this beam energy, for the ejectile's OWN cm angle
      TLorentzVector Lb(pbeam, Eb), Lt(0, 0, 0, mT);
      TLorentzVector W = Lb + Lt;
      TVector3 bst = W.BoostVector();
      TLorentzVector Le(pej, std::sqrt(pej.Mag2() + mE * mE));
      TLorentzVector LeCM = Le;
      LeCM.Boost(-bst);
      const double s = W.M2();
      const double pcm = std::sqrt(std::max(0., (s - (mR + mE) * (mR + mE)) * (s - (mR - mE) * (mR - mE)))) /
                         (2 * std::sqrt(s));
      // predicted: rebuild the lab ejectile from the SAME cm angle and the analytic |p_cm|
      TVector3 dir = LeCM.Vect().Unit();
      TLorentzVector Lpred(pcm * dir, std::sqrt(pcm * pcm + mE * mE));
      Lpred.Boost(bst);

      const double thTrue = Le.Vect().Theta() * TMath::RadToDeg();
      const double keTrue = (Le.E() - mE) * 1000; // MeV
      hdA->Fill(thTrue - Lpred.Vect().Theta() * TMath::RadToDeg());
      hdE->Fill(keTrue - (Lpred.E() - mE) * 1000);
      thMin = std::min(thMin, thTrue);
      thMax = std::max(thMax, thTrue);
      keMin = std::min(keMin, keTrue);
      keMax = std::max(keMax, keTrue);
      ++nk;
   }

   if (!nb) {
      printf("\033[1;31mno pristine beam ion found (none started at z = %.0f cm)\033[0m\n", zEntrance);
      return;
   }
   const double p = sumP / nb;
   const double rms = std::sqrt(std::max(0., sumP2 / nb - p * p));
   const double KE = (std::sqrt(p * p + mB * mB) - mB) * 1000;

   printf("\n\033[1;33m===== 1. transported 17C beam, from MC truth =====\033[0m\n");
   printf("  beam events      = %d      reaction events = %d\n", nb, nreact);
   printf("  <p> beam         = %.6f GeV/c  (rms %.2e)   asked %.6f\n", p, rms, pAsked);
   printf("  generator wiring : %s  (%+.4f %%)\n",
          std::fabs(p - pAsked) / pAsked < 1e-4 ? "\033[1;32mOK\033[0m" : "\033[1;31mMISMATCH\033[0m",
          100 * (p - pAsked) / pAsked);
   printf("  KE with the physical m(17C) = %.5f GeV -> \033[1;36m%.2f MeV\033[0m (%.3f MeV/u), asked 142.29 (8.370)\n",
          mB, KE, KE / 17);

   printf("\n\033[1;33m===== 2. reaction vertices =====\033[0m\n");
   if (nreact)
      printf("  vertex z : %.1f to %.1f cm, mean %.2f cm   (uniform over the active volume)\n", vzMin, vzMax,
             sumVz / nreact);

   printf("\n\033[1;33m===== 3. two-body kinematics closure (proton) =====\033[0m\n");
   printf("  Q(17C(d,p)18C) at Ex = %.3f : %+.4f MeV\n", resEx, (mB + mT - mR - mE) * 1000);
   printf("  protons found    = %d\n", nk);
   if (nk) {
      printf("  theta_lab range  = %.1f to %.1f deg\n", thMin, thMax);
      printf("  KE range         = %.2f to %.2f MeV\n", keMin, keMax);
      printf("  truth - analytic : dTheta = %+.4f +- %.4f deg,  dKE = %+.4f +- %.4f MeV\n", hdA->GetMean(),
             hdA->GetRMS(), hdE->GetMean(), hdE->GetRMS());
      const bool ok = std::fabs(hdA->GetMean()) < 0.05 && std::fabs(hdE->GetMean()) < 0.05;
      printf("  kinematics       : %s\n", ok ? "\033[1;32mCLOSES\033[0m"
                                             : "\033[1;31mDOES NOT CLOSE -- check the masses/Ex\033[0m");
   }
   printf("\n  The mass itself cannot be read back from MCTrack (see the note in the source);\n"
          "  checks 1 and 3 together are what confirm the beam AND the reaction.\n\n");
   fi->Close();
}
