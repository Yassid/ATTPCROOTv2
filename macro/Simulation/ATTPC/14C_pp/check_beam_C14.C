/// @file check_beam_C14.C
/// @brief Read the beam kinetic energy back out of the MC truth of C14_pp_sim.C.
///
/// Every AT-TPC sim macro passes the ion mass in amu to AtTPCIonGenerator, whose FairIon
/// mass parameter is documented as GeV. If FairIon actually honours that number, the
/// transported 14C would have m = 14.003 GeV instead of the true 13.044 GeV, and the beam
/// energy implied by pz = 2.05574 GeV/c would come out at ~150 MeV rather than 161. This
/// macro does not argue about it -- it reads the primary track's four-momentum from the
/// simulated file and prints the energy that was actually transported.
///
///   root -b -q 'check_beam_C14.C("./data/attpcsim.root")'

void check_beam_C14(TString f = "./data/attpcsim.root")
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
   t->SetBranchAddress("MCTrack", &tracks); // FairRoot names the branch MCTrack

   const double U = 0.93149401;          // GeV per amu
   const double mTrue = 14.003242 * U;   // 13.04394 GeV -- the physical 14C mass
   const double mAmuAsGeV = 14.003242;   // what FairIon gets if the amu number is taken literally

   // AtMCTrack::GetEnergy() builds E from a TDatabasePDG mass lookup, and ion PDG codes
   // (10LZZZAAAI) are not in that table -> it returns E = |p| and an effective mass of 0.
   // So the mass cannot be read back from the track; only the MOMENTUM is trustworthy.
   // AtVertexPropagator alternates beam-only and reaction events. In a BEAM event track 0
   // is the pristine beam ion, generated at the entrance (z = -100 cm). In a REACTION
   // event track 0 is the SCATTERED ion, already sitting at the vertex with a degraded
   // momentum -- averaging over both drags the mean ~1 % low and looks like a wiring bug.
   // Select on the start position, which is what actually distinguishes them.
   const double zEntrance = -100.0; // cm, from ionGen->SetSpotRadius(0, -100, 0)
   int nb = 0, nreact = 0;
   double sumP = 0, sumP2 = 0, sumPr = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (tracks->GetEntriesFast() == 0)
         continue;
      auto *tr = (AtMCTrack *)tracks->At(0);
      if (!tr || tr->GetPdgCode() < 1000000000)
         continue;
      const double p =
         std::sqrt(tr->GetPx() * tr->GetPx() + tr->GetPy() * tr->GetPy() + tr->GetPz() * tr->GetPz());
      if (std::fabs(tr->GetStartZ() - zEntrance) < 1.0) { // pristine beam
         sumP += p;
         sumP2 += p * p;
         ++nb;
      } else { // scattered ion at the reaction vertex
         sumPr += p;
         ++nreact;
      }
   }
   if (!nb) {
      printf("\033[1;31mno pristine beam ion found (none started at z = %.0f cm)\033[0m\n", zEntrance);
      return;
   }
   const double p = sumP / nb, pAsked = 2.055740;
   const double rms = std::sqrt(std::max(0., sumP2 / nb - p * p));

   printf("\n\033[1;33m===== transported 14C beam, from MC truth =====\033[0m\n");
   printf("  beam events      = %d      reaction events = %d\n", nb, nreact);
   printf("  <p> beam         = %.6f GeV/c  (rms %.2e)   asked %.6f\n", p, rms, pAsked);
   printf("  generator wiring : %s  (%+.3f %%)\n",
          std::fabs(p - pAsked) / pAsked < 1e-4 ? "\033[1;32mOK\033[0m" : "\033[1;31mMISMATCH\033[0m",
          100 * (p - pAsked) / pAsked);
   if (nreact)
      printf("  <p> scattered 14C at the vertex = %.6f GeV/c  (degraded, as expected)\n", sumPr / nreact);
   printf("\n  KE implied by that momentum:\n");
   printf("    with the physical m(14C) = %.5f GeV  ->  \033[1;36m%.2f MeV\033[0m  (%.3f MeV/u)\n", mTrue,
          (std::sqrt(p * p + mTrue * mTrue) - mTrue) * 1000, (std::sqrt(p * p + mTrue * mTrue) - mTrue) * 1000 / 14);
   printf("    if the amu number were taken as GeV (%.5f) -> %.2f MeV\n", mAmuAsGeV,
          (std::sqrt(p * p + mAmuAsGeV * mAmuAsGeV) - mAmuAsGeV) * 1000);
   printf("\n  The mass itself cannot be read back from MCTrack (see note above); this only\n"
          "  confirms the momentum handed to the transport. Compare ranges/dE-dx against the\n"
          "  data to settle which mass Geant4 actually used.\n\n");
   fi->Close();
}
