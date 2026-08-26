/// @file pid2csv_C15d.C
/// @brief Export a cached PID plane to CSV, so it can be compared against Spyral's parquet.
///
///   root -b -q 'tools/pid2csv_C15d.C("run_0017", "/home/yassid/C15d_reco/", "/tmp/at_0017.csv")'
///
/// Exists only because there is no uproot in the Spyral venv and no parquet reader in ROOT,
/// so CSV is the one format both sides can read. Valid tracks only -- the comparison is about
/// where the bands sit, and Spyral's estimate table contains only tracks it accepted too.

void pid2csv_C15d(TString fileName = "run_0017", TString inDir = "/home/yassid/C15d_reco/",
                  TString outCsv = "")
{
   TString inputFile = inDir + fileName + "_pid.root";
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }
   if (outCsv.Length() == 0)
      outCsv = inDir + fileName + "_pid.csv";

   TFile *f = TFile::Open(inputFile);
   auto *t = dynamic_cast<TTree *>(f->Get("pid"));
   if (t == nullptr) {
      std::cout << "\033[1;31mERROR: no 'pid' tree in " << inputFile << "\033[0m\n";
      return;
   }

   Int_t run, event, trackID, valid, nPoints, nClusters, direction;
   Double_t dEdx, sqrtdEdx, dE, arclength, brho, polar, azimuthal, radius, vtxZ, vtxR;
   t->SetBranchAddress("run", &run);
   t->SetBranchAddress("event", &event);
   t->SetBranchAddress("trackID", &trackID);
   t->SetBranchAddress("valid", &valid);
   t->SetBranchAddress("nPoints", &nPoints);
   t->SetBranchAddress("nClusters", &nClusters);
   t->SetBranchAddress("direction", &direction);
   t->SetBranchAddress("dEdx", &dEdx);
   t->SetBranchAddress("sqrtdEdx", &sqrtdEdx);
   t->SetBranchAddress("dE", &dE);
   t->SetBranchAddress("arclength", &arclength);
   t->SetBranchAddress("brho", &brho);
   t->SetBranchAddress("polar", &polar);
   t->SetBranchAddress("azimuthal", &azimuthal);
   t->SetBranchAddress("radius", &radius);
   t->SetBranchAddress("vtxZ", &vtxZ);
   t->SetBranchAddress("vtxR", &vtxR);

   std::ofstream out(outCsv.Data());
   out << "run,event,trackID,direction,nPoints,nClusters,dEdx,sqrt_dEdx,dE,arclength,brho,polar,azimuthal,"
          "radius,vtxZ,vtxR\n";
   Long64_t n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (valid != 1)
         continue;
      out << run << "," << event << "," << trackID << "," << direction << "," << nPoints << "," << nClusters
          << "," << dEdx << "," << sqrtdEdx << "," << dE << "," << arclength << "," << brho << "," << polar
          << "," << azimuthal << "," << radius << "," << vtxZ << "," << vtxR << "\n";
      ++n;
   }
   out.close();
   f->Close();
   std::cout << "\033[1;32mwrote\033[0m " << outCsv << "  (" << n << " valid tracks)\n";
}
