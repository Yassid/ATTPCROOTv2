/// @file pid_plane_a1975.C
/// @brief Standalone, high-resolution PID plane (sqrt_dEdx vs brho) for a1975.
///
/// Computes the Spyral PID observables (AtTools::AtPIDEstimator) for every track
/// in <run>_reco.root and draws ONE large sqrt_dEdx-vs-brho histogram so the
/// proton / deuteron / triton bands are clearly visible. Also caches the
/// observables to <run>_pidobs.root (a TNtuple) so the plot can be re-made
/// instantly (replot=true) without re-reading the big reco file.
///
/// Run:
///   root -b -q 'pid_plane_a1975.C("run_0116")'                 // full, brho 0..1.5
///   root -b -q 'pid_plane_a1975.C("run_0116", 15, 1.5, true)'  // fast re-plot from cache

void pid_plane_a1975(TString fileName = "run_0116", Double_t dedxMax = 12000.0, Double_t brMax = 3.0,
                     Bool_t replot = false, Int_t minClusters = 0, Double_t bField = 2.85,
                     Double_t smallPadRadius = 152.0)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString obsFile = fileName + "_pidobs.root";

   // Spyral-style axes: dEdx (counts) on x, brho on y.
   TH2F *h = new TH2F("hpid", "Particle ID  (C++ Spyral port);dEdx  [counts];B#rho  [T m]", 400, 0, dedxMax, 400, 0,
                      brMax);

   if (replot && gSystem->AccessPathName(obsFile) == 0) {
      // Fast path: read cached observables.
      TFile *fo = TFile::Open(obsFile);
      TNtuple *nt = (TNtuple *)fo->Get("pidobs");
      float dedx, br, ncl;
      nt->SetBranchAddress("dedx", &dedx);
      nt->SetBranchAddress("brho", &br);
      nt->SetBranchAddress("ncl", &ncl);
      Long64_t kept = 0;
      for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
         nt->GetEntry(i);
         if (ncl < minClusters)
            continue;
         h->Fill(dedx, br);
         ++kept;
      }
      printf("Re-plotted %lld / %lld cached tracks (minClusters=%d)\n", kept, nt->GetEntries(), minClusters);
   } else {
      // Compute from the reco file and cache.
      TString inputFile = fileName + "_reco.root";
      if (gSystem->AccessPathName(inputFile)) {
         printf("\033[1;31mERROR: %s not found\033[0m\n", inputFile.Data());
         return;
      }
      TFile *f = TFile::Open(inputFile);
      TTree *t = (TTree *)f->Get("cbmsim");
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      TClonesArray *peArr = nullptr;
      t->SetBranchAddress("AtPatternEvent", &peArr);

      AtTools::AtPIDEstimator estimator(bField, smallPadRadius);

      TFile *fo = new TFile(obsFile, "RECREATE");
      TNtuple *nt = new TNtuple("pidobs", "PID observables", "sqrtdedx:brho:dedx:polar:dE:arclen:ncl");
      Long64_t n = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (peArr->GetEntries() == 0)
            continue;
         auto *pe = (AtPatternEvent *)peArr->At(0);
         if (!pe)
            continue;
         for (auto &track : pe->GetTrackCand()) {
            auto r = estimator.Estimate(const_cast<AtTrack &>(track));
            if (!r.valid)
               continue;
            nt->Fill(r.sqrtdEdx, r.brho, r.dEdx, r.polar * TMath::RadToDeg(), r.dE, r.arclength, r.nClusters);
            if (r.nClusters < minClusters)
               continue;
            h->Fill(r.dEdx, r.brho);
            ++n;
         }
      }
      nt->Write();
      fo->Close();
      f->Close();
      printf("Computed %lld tracks, cached to %s\n", n, obsFile.Data());
   }

   TCanvas *c = new TCanvas("c", "pid_plane", 1200, 1000);
   c->SetRightMargin(0.13);
   c->SetLogz();
   h->Draw("colz");
   TString png = "pid/plots/" + fileName + "_pid_plane.png";
   c->SaveAs(png);
   printf("Saved %s  (dedxMax=%.0f, brMax=%.2f)\n", png.Data(), dedxMax, brMax);
}
