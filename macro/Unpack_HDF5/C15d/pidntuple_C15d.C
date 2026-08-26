/// @file pidntuple_C15d.C
/// @brief Cache the gain-matched PID plane of one reco as a small flat tree.
///
/// Reads <run>_reco.root (which already carries a gain-matched AtPIDEvent, computed once by
/// unpackReco_C15d.C) and writes <run>_pid.root holding one entry PER TRACK.
///
///   root -b -q 'pidntuple_C15d.C("run_0017", "/home/yassid/C15d_reco/", "/home/yassid/C15d_ic/")'
///
/// Why cache at all: the plane is a few MB per run against a multi-GB reco, so once this exists
/// every re-plot, re-binning and gate iteration is instant instead of a minutes-long re-read.
///
/// ★ THE CACHE KEEPS `event` AND `trackID`. That is the whole point and it is easy to get wrong:
/// a cache of only (sqrtdEdx, brho, ...) lets you DRAW a gate but not APPLY one, because there is
/// nothing to join back to the reco on. With both kept, a gate selected here can be turned into a
/// per-event, per-track selection for the fitting stage without recomputing AtSpyralPID.
///
/// `valid` is stored rather than filtered on, so the rejected population stays countable -- if a
/// band looks thin, the first question is how many tracks never made it into the plane at all.

void pidntuple_C15d(TString fileName = "run_0017", TString inDir = "/home/yassid/C15d_reco/",
                    TString outDir = "", TString recoSuffix = "_reco", TString treeName = "cbmsim")
{
   gSystem->Load("libAtReconstruction.so");

   TString inputFile = inDir + fileName + recoSuffix + ".root";
   TString outBase = (outDir.Length() ? outDir : inDir);
   TString outputFile = outBase + fileName + "_pid.root";

   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }

   TFile *fin = TFile::Open(inputFile);
   if (fin == nullptr || fin->IsZombie()) {
      std::cout << "\033[1;31mERROR: cannot open " << inputFile << "\033[0m\n";
      return;
   }
   auto *tree = dynamic_cast<TTree *>(fin->Get(treeName));
   if (tree == nullptr) {
      std::cout << "\033[1;31mERROR: no tree '" << treeName << "' in " << inputFile << "\033[0m\n";
      return;
   }
   if (tree->GetBranch("AtPIDEvent") == nullptr) {
      std::cout << "\033[1;31mERROR: " << inputFile << " has no AtPIDEvent branch.\n"
                << "  Re-reco with doPID=kTRUE (unpackReco_C15d.C default).\033[0m\n";
      return;
   }

   TClonesArray *pidArray = nullptr;
   tree->SetBranchAddress("AtPIDEvent", &pidArray);

   const Int_t runNo = AtGainMatchTask::RunNumberFromName(fileName);

   TFile *fout = TFile::Open(outputFile, "RECREATE");
   TTree *out = new TTree("pid", "gain-matched PID plane, one entry per track");

   Int_t t_run = runNo, t_event = 0, t_trackID = -1, t_direction = -1;
   Int_t t_nPoints = 0, t_nClusters = 0, t_valid = 0;
   Double_t t_dEdx = 0, t_sqrtdEdx = 0, t_dE = 0, t_arclength = 0, t_brho = 0;
   Double_t t_polar = 0, t_azimuthal = 0, t_radius = 0;
   Double_t t_vtxX = 0, t_vtxY = 0, t_vtxZ = 0, t_vtxR = 0;

   out->Branch("run", &t_run, "run/I");
   out->Branch("event", &t_event, "event/I");
   out->Branch("trackID", &t_trackID, "trackID/I");
   out->Branch("valid", &t_valid, "valid/I");
   out->Branch("direction", &t_direction, "direction/I");
   out->Branch("nPoints", &t_nPoints, "nPoints/I");
   out->Branch("nClusters", &t_nClusters, "nClusters/I");
   out->Branch("dEdx", &t_dEdx, "dEdx/D");
   out->Branch("sqrtdEdx", &t_sqrtdEdx, "sqrtdEdx/D");
   out->Branch("dE", &t_dE, "dE/D");
   out->Branch("arclength", &t_arclength, "arclength/D");
   out->Branch("brho", &t_brho, "brho/D");
   out->Branch("polar", &t_polar, "polar/D");
   out->Branch("azimuthal", &t_azimuthal, "azimuthal/D");
   out->Branch("radius", &t_radius, "radius/D");
   out->Branch("vtxX", &t_vtxX, "vtxX/D");
   out->Branch("vtxY", &t_vtxY, "vtxY/D");
   out->Branch("vtxZ", &t_vtxZ, "vtxZ/D");
   out->Branch("vtxR", &t_vtxR, "vtxR/D");

   const Long64_t nEntries = tree->GetEntries();
   Long64_t nTracks = 0, nValid = 0;

   for (Long64_t i = 0; i < nEntries; ++i) {
      tree->GetEntry(i);
      if (pidArray == nullptr || pidArray->GetEntriesFast() == 0)
         continue;
      auto *ev = dynamic_cast<AtPIDEvent *>(pidArray->At(0));
      if (ev == nullptr)
         continue;
      t_event = static_cast<Int_t>(i);
      for (const auto &r : ev->GetSpyral()) {
         t_trackID = r.trackID;
         t_valid = r.valid ? 1 : 0;
         t_direction = r.direction;
         t_nPoints = r.nPoints;
         t_nClusters = r.nClusters;
         t_dEdx = r.dEdx;
         t_sqrtdEdx = r.sqrtdEdx;
         t_dE = r.dE;
         t_arclength = r.arclength;
         t_brho = r.brho;
         t_polar = r.polar;
         t_azimuthal = r.azimuthal;
         t_radius = r.radius;
         t_vtxX = r.vertex.X();
         t_vtxY = r.vertex.Y();
         t_vtxZ = r.vertex.Z();
         t_vtxR = std::sqrt(r.vertex.X() * r.vertex.X() + r.vertex.Y() * r.vertex.Y());
         out->Fill();
         ++nTracks;
         if (r.valid)
            ++nValid;
      }
   }

   fout->cd();
   out->Write();
   fout->Close();
   fin->Close();

   std::cout << "\033[1;32m" << fileName << "\033[0m: " << nEntries << " events, " << nTracks << " tracks, "
             << nValid << " valid (" << (nTracks ? 100.0 * nValid / nTracks : 0.0) << "%) -> " << outputFile
             << std::endl;
}
