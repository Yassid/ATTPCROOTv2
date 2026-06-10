/// @file pid_spyral_a1975.C
/// @brief PID plane from the FAITHFUL Spyral estimator (AtTools::AtSpyralPID),
///        which isolates the first arc, fits a dedicated circle, gets polar from
///        a rho-vs-z regression, and integrates arc length via smoothing splines.
///        Compares yield against the old AtPIDEstimator recipe.
///
///   root -b -q 'pid_spyral_a1975.C("run_0116")'
///   root -b -q 'pid_spyral_a1975.C("run_0116", 3000)'   // first 3000 events

void pid_spyral_a1975(TString fileName = "run_0116", Long64_t nEvents = -1, double bField = 2.85)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *f = TFile::Open(fileName + "_reco.root");
   TTree *t = (TTree *)f->Get("cbmsim");
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   TClonesArray *peArr = nullptr;
   t->SetBranchAddress("AtPatternEvent", &peArr);
   Long64_t nE = t->GetEntries();
   if (nEvents > 0 && nEvents < nE)
      nE = nEvents;

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);
   AtTools::AtPIDEstimator old(bField, 152.0);

   TH2F *hS = new TH2F("hS", "Spyral PID;#sqrt{dEdx};B#rho [T m]", 350, 0, 60, 350, 0, 0.5);
   TH2F *hSd = new TH2F("hSd", "Spyral PID;dEdx [counts];B#rho [T m]", 350, 0, 1500, 350, 0, 0.5);
   TH2F *hO = new TH2F("hO", "Old recipe;dEdx [counts];B#rho [T m]", 350, 0, 1500, 350, 0, 0.5);

   TFile *fo = new TFile(fileName + "_spyralpid.root", "RECREATE");
   TNtuple *nt =
      new TNtuple("spid", "spyral pid", "dedx:sqrtdedx:brho:polar:dE:arclen:npts:dir:vtxr:vtxz:radius");

   long nTracks = 0, nValidS = 0, nValidO = 0;
   for (Long64_t i = 0; i < nE; ++i) {
      t->GetEntry(i);
      if (peArr->GetEntries() == 0)
         continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe)
         continue;
      for (auto &trk : pe->GetTrackCand()) {
         ++nTracks;
         AtTrack &tr = const_cast<AtTrack &>(trk);
         auto r = spy.Estimate(tr);
         if (r.valid) {
            ++nValidS;
            double th = r.polar * TMath::RadToDeg();
            hS->Fill(r.sqrtdEdx, r.brho);
            hSd->Fill(r.dEdx, r.brho);
            double vtxr = TMath::Hypot(r.vertex.X(), r.vertex.Y());
            nt->Fill(r.dEdx, r.sqrtdEdx, r.brho, th, r.dE, r.arclength, r.nPoints, r.direction, vtxr, r.vertex.Z(),
                     r.radius);
         }
         auto ro = old.Estimate(tr);
         if (ro.valid) {
            ++nValidO;
            hO->Fill(ro.dEdx, ro.brho);
         }
      }
   }
   nt->Write();

   printf("\n===== SPYRAL PID  %s  (%lld events) =====\n", fileName.Data(), nE);
   printf("tracks=%ld   Spyral-valid=%ld (%.1f%%)   old-recipe-valid=%ld (%.1f%%)\n", nTracks, nValidS,
          nTracks ? 100.0 * nValidS / nTracks : 0, nValidO, nTracks ? 100.0 * nValidO / nTracks : 0);

   TCanvas *c = new TCanvas("c", "spyralpid", 1500, 500);
   c->Divide(3, 1);
   c->cd(1);
   hS->Draw("colz");
   c->cd(2);
   hSd->Draw("colz");
   c->cd(3);
   hO->Draw("colz");
   c->SaveAs("pid/plots/" + fileName + "_spyral_pid.png");
   printf("saved %s_spyral_pid.png  (left/mid: Spyral, right: old recipe)\n", fileName.Data());
   fo->Close();
}
