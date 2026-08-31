/// @file pid_points_3Hed.C
/// @brief Rebuild the 46Ar(3He,d) PID plane at a CHOSEN fMinPoints, from the reco files. No truth,
/// no gate, no selection -- every track that AtSpyralPID accepts at that cut goes in.
///
/// WHY THIS EXISTS. The chain's own _pid.root files are written by AtPIDTask, which never calls
/// SetMinPoints, so they are fixed at the class default fMinPoints = 30 -- and at 30 the deuterons
/// between theta_lab 77 and 104 deg are ALL rejected (measured: 0 of 89 on the probe sample, every
/// one of them failCode 1, too few points after the z-tie collapse). That band is theta_cm 30-55
/// deg, the middle of the proposal's own 15-80 deg range. A gate drawn on the 30 plane is drawn on
/// a plane with a hole where the deuterons near 90 deg lab should be: they run nearly
/// perpendicular to the drift axis, span almost no z, and collapse to too few distinct-z knots.
///
/// This re-runs AtSpyralPID::Estimate() in memory over the pattern tracks in *_reco.root, so the
/// cut is a parameter. Nothing is re-reconstructed and nothing is re-fitted -- it reads the same
/// reco files the chain already wrote.
///
/// Sweep the cut before settling on one: on a1975 the knee was 15-20, below which the direction
/// checks (codes 8, 9) start firing instead. That number was measured for protons in H2 and does
/// NOT transfer here -- measure it on this reaction, which is what the printed table is for.
///
///   root -b -q 'pid_points_3Hed.C(15)'                    // one cut
///   root -b -q 'pid_points_3Hed.C(-1)'                    // class default (30), for comparison
///
/// Output: plots/pid_points_mp<N>.root, a TTree "pts" of (x = sqrt(dE/dx), y = Brho, nClusters,
/// nPoints, polar_deg, vertex_z, entry, trackID) -- everything a gate needs, plus the vertex the
/// pre-fit kinematics need, plus enough to go back to the track it came from.

void pid_points_3Hed(Int_t minPoints = 15, TString outFile = "",
                     TString simDir = "/mnt/f/ar46_3hed_OLD_2.85T_placeholder",
                     TString tags = "gs_s3001,gs_s3002,360_s3011,360_s3012,2020_s3021,2020_s3022",
                     Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0)
      spy.SetMinPoints(minPoints);
   const int mp = minPoints > 0 ? minPoints : 30;
   if (outFile.IsNull())
      outFile = TString::Format("plots/pid_points_mp%d.root", mp);
   printf("\n  fMinPoints = %d\n", mp);

   std::vector<float> vx, vy, vpol, vvz;
   std::vector<int> vnc, vnp, vent, vtid;
   long nTracks = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      if (tg.IsNull())
         continue;
      TString fr = simDir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fr)) {
         printf("  skip %-12s (missing)\n", tg.Data());
         continue;
      }
      TFile *Fr = TFile::Open(fr);
      TTree *tr = Fr ? (TTree *)Fr->Get("cbmsim") : nullptr;
      if (!tr) {
         if (Fr)
            Fr->Close();
         continue;
      }
      TClonesArray *pa = nullptr;
      tr->SetBranchAddress("AtPatternEvent", &pa);
      long nAcc = 0, nAll = 0;
      for (Long64_t i = 0; i < tr->GetEntries(); ++i) {
         tr->GetEntry(i);
         if (!pa || !pa->GetEntriesFast())
            continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe)
            continue;
         for (auto &track : pe->GetTrackCand()) {
            AtTrack &t2 = const_cast<AtTrack &>(track);
            ++nAll;
            auto res = spy.Estimate(t2);
            if (!res.valid)
               continue;
            ++nAcc;
            vx.push_back(res.sqrtdEdx);
            vy.push_back(res.brho);
            vpol.push_back(res.polar * TMath::RadToDeg());
            vvz.push_back(res.vertex.Z()); // mm, reconstructed frame -- see kinematics_3Hed.C on handedness
            vnc.push_back(res.nClusters);
            vnp.push_back(res.nPoints);
            vent.push_back((int)i);
            vtid.push_back(res.trackID);
         }
      }
      nTracks += nAll;
      printf("  %-12s %6ld of %6ld pattern tracks accepted (%.0f %%)\n", tg.Data(), nAcc, nAll,
             nAll ? 100.0 * nAcc / nAll : 0.);
      Fr->Close();
   }
   delete ta;

   if (vx.empty()) {
      printf("\n  nothing accepted -- check that the reco files exist\n");
      return;
   }

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   gSystem->mkdir(here + "/plots", kTRUE);
   TFile out(here + "/" + outFile, "RECREATE");
   TTree t("pts", "46Ar(3He,d) PID points");
   float x, y, pol, vz;
   int nc, np, ent, tid;
   t.Branch("x", &x);
   t.Branch("y", &y);
   t.Branch("polar_deg", &pol);
   t.Branch("vertex_z", &vz);
   t.Branch("nClusters", &nc);
   t.Branch("nPoints", &np);
   t.Branch("entry", &ent);
   t.Branch("trackID", &tid);
   for (size_t i = 0; i < vx.size(); ++i) {
      x = vx[i]; y = vy[i]; pol = vpol[i]; vz = vvz[i]; nc = vnc[i]; np = vnp[i]; ent = vent[i]; tid = vtid[i];
      t.Fill();
   }
   t.Write();
   out.Close();

   // Where the added tracks land is the whole point of lowering the cut. The dead band is the
   // MIDDLE of the rigidity range here, not an end of it -- 12 to 30 MeV of deuteron is
   // Brho 0.66 to 1.10 T*m -- which is why it is quoted as a window rather than a tail.
   long mid = 0;
   for (size_t i = 0; i < vy.size(); ++i)
      if (vy[i] > 0.66 && vy[i] < 1.10)
         ++mid;
   printf("\n  %zu points from %ld pattern tracks (%.1f %%)\n", vx.size(), nTracks,
          100.0 * vx.size() / std::max(1L, nTracks));
   printf("  in the mp30 dead band, Brho 0.66-1.10 T*m: %ld (%.1f %%)\n", mid,
          100.0 * mid / std::max(1UL, vx.size()));
   printf("  wrote %s\n\n", outFile.Data());
}
