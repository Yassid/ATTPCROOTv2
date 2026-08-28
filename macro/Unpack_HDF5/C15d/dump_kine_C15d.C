/// @file dump_kine_C15d.C
/// @brief Flatten a GENFIT fit output into a small kinematics ntuple, selected by a PID gate.
///
///   root -b -q 'dump_kine_C15d.C("run_0017","_dd","/path/to/fits/","pid/sel_deuteron_C15d.root")'
///
/// Output tree `kin`, one entry per fitted track that passes:
///   run event trackID | ke (MeV) theta (deg) phi (deg) | vx vy vz (mm) | chi2ndf ndf | dirFwd
///
/// ★ THE GATE IS APPLIED HERE, BY JOINING ON (run, event, trackID) -- NOT inside the fitter.
/// AtGenfitter::SetPIDGate evaluates its own internal AtSpyralPID on RAW dE/dx, while every gate
/// in this workspace is drawn on the GAIN-MATCHED plane. Those are different x scales (~29 % of
/// per-run drift), so an in-fit gate would select a different set of tracks than the one you drew
/// and would do it silently. Fit everything, select afterwards on the identity triple.
///
/// ★ theta IS AS THE FITTER REPORTS IT and is NOT folded into the forward hemisphere. For (d,p)
/// the proton goes backward in the lab, so theta > 90 deg is physical and must survive to the
/// kinematics. `dirFwd` records the sign of the fitted momentum's z component so a backward
/// population can be separated from a forward one that was mis-seeded.
///
/// Both GetKinematics() (at the first fitted point) and GetKinematicsXtr() (back-extrapolated to
/// the beam axis) are written: keXtr/thetaXtr is the one to use for physics, ke/theta is kept so
/// the size of the back-extrapolation correction stays visible instead of being folded in.

void dump_kine_C15d(TString fileName = "run_0017", TString suffix = "_dd",
                    TString inDir = "/home/yassid/C15d_fit/", TString selFile = "",
                    TString species = "d", TString outFile = "", Double_t chi2Cut = 1e9)
{
   gSystem->Load("libAtReconstruction.so");

   TString in = inDir + fileName + "_genfit_" + species + suffix + ".root";
   if (gSystem->AccessPathName(in)) {
      std::cout << "\033[1;31mERROR: " << in << " not found.\033[0m\n";
      return;
   }
   const Int_t runNo = AtGainMatchTask::RunNumberFromName(fileName);

   // ---- selection: (run, event, trackID) triples from apply_gate_C15d.C -------------------
   std::set<Long64_t> keep;
   auto key = [](int r, int e, int t) { return ((Long64_t)r << 44) | ((Long64_t)e << 12) | (Long64_t)(t & 0xFFF); };
   if (selFile.Length()) {
      if (gSystem->AccessPathName(selFile)) {
         std::cout << "\033[1;31mERROR: selection " << selFile << " not found.\033[0m\n";
         return;
      }
      TFile *fs = TFile::Open(selFile);
      TTree *ts = fs ? (TTree *)fs->Get("sel") : nullptr;
      if (!ts) {
         std::cout << "\033[1;31mERROR: no tree 'sel' in " << selFile << "\033[0m\n";
         return;
      }
      Int_t r, e, t;
      ts->SetBranchAddress("run", &r);
      ts->SetBranchAddress("event", &e);
      ts->SetBranchAddress("trackID", &t);
      for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
         ts->GetEntry(i);
         if (r == runNo)
            keep.insert(key(r, e, t));
      }
      fs->Close();
      if (keep.empty()) {
         std::cout << "\033[1;33mWARNING: the selection contains no tracks for run " << runNo
                   << " -- nothing will pass.\033[0m\n";
      }
   }

   TFile *f = TFile::Open(in);
   TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
   if (!t) {
      std::cout << "\033[1;31mERROR: no tree 'cbmsim' in " << in << "\033[0m\n";
      return;
   }
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);

   if (outFile.IsNull())
      outFile = inDir + fileName + "_kin_" + species + suffix + ".root";
   TFile fo(outFile, "RECREATE");
   TTree kin("kin", "fitted kinematics");
   Int_t o_run = runNo, o_event, o_track, o_ndf, o_dirFwd;
   Double_t o_ke, o_th, o_phi, o_keX, o_thX, o_phiX, o_vx, o_vy, o_vz, o_c2n;
   kin.Branch("run", &o_run, "run/I");
   kin.Branch("event", &o_event, "event/I");
   kin.Branch("trackID", &o_track, "trackID/I");
   kin.Branch("ke", &o_ke, "ke/D");
   kin.Branch("theta", &o_th, "theta/D");
   kin.Branch("phi", &o_phi, "phi/D");
   kin.Branch("keXtr", &o_keX, "keXtr/D");
   kin.Branch("thetaXtr", &o_thX, "thetaXtr/D");
   kin.Branch("phiXtr", &o_phiX, "phiXtr/D");
   kin.Branch("vx", &o_vx, "vx/D");
   kin.Branch("vy", &o_vy, "vy/D");
   kin.Branch("vz", &o_vz, "vz/D");
   kin.Branch("chi2ndf", &o_c2n, "chi2ndf/D");
   kin.Branch("ndf", &o_ndf, "ndf/I");
   kin.Branch("dirFwd", &o_dirFwd, "dirFwd/I");

   Long64_t nFit = 0, nSel = 0, nOut = 0, nBadNdf = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!te || te->GetEntries() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev)
         continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         ++nFit;
         o_event = (Int_t)i;
         o_track = ft->GetTrackID();
         if (!keep.empty() && keep.find(key(runNo, o_event, o_track)) == keep.end())
            continue;
         ++nSel;

         const auto &md = ft->GetTrackMetadata();
         const double ndf = md ? md->GetNdf() : -1, chi2 = md ? md->GetChi2() : -1;
         // ndf <= 0 is a COLLAPSED fit: genfit still writes kinematics for it, and they look
         // ordinary. Counting them separately keeps a failure population from quietly becoming a
         // physics population.
         if (!(ndf > 0)) {
            ++nBadNdf;
            continue;
         }
         o_ndf = (Int_t)ndf;
         o_c2n = chi2 / ndf;
         if (o_c2n > chi2Cut)
            continue;

         const auto &k = ft->GetKinematics();
         const auto &kx = ft->GetKinematicsXtr();
         o_ke = k.kineticEnergy;
         o_th = k.theta * TMath::RadToDeg();
         o_phi = k.phi * TMath::RadToDeg();
         o_keX = kx.kineticEnergy;
         o_thX = kx.theta * TMath::RadToDeg();
         o_phiX = kx.phi * TMath::RadToDeg();
         const auto &tp = ft->GetTrackPropertiesStruct();
         o_vx = tp.initialPositionXtr.X();
         o_vy = tp.initialPositionXtr.Y();
         o_vz = tp.initialPositionXtr.Z();
         // theta is NOT folded: > 90 deg is a genuinely backward track, which (d,p) needs.
         o_dirFwd = (o_thX >= 0 && o_thX < 90) ? 1 : 0;
         kin.Fill();
         ++nOut;
      }
   }
   fo.cd();
   kin.Write();
   fo.Close();
   f->Close();

   std::cout << "\033[1;32m" << fileName << "\033[0m: " << nFit << " fitted";
   if (!keep.empty())
      std::cout << ", " << nSel << " in gate";
   std::cout << ", " << nBadNdf << " collapsed (ndf<=0), " << nOut << " written -> " << outFile << "\n";
}
