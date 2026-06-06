/// @file ukf_results_a1975.C
/// @brief Full-statistics UKF results for a1975 16C+p: KE spectrum, KE-vs-theta
/// kinematics, and chi2/ndf. Reads <run>_ukf.root (AtTrackingEvent).
///
/// Run: root -b -q 'ukf_results_a1975.C("run_0116")'

void ukf_results_a1975(TString fileName = "run_0116", Double_t chi2Cut = 2.0, Double_t keMax = 40.0)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString inputFile = fileName + "_ukf.root";
   if (gSystem->AccessPathName(inputFile)) {
      printf("\033[1;31mERROR: %s not found\033[0m\n", inputFile.Data());
      return;
   }
   TFile *f = TFile::Open(inputFile);
   TTree *t = (TTree *)f->Get("cbmsim");
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtTrackingEvent*", 1);
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);

   TH2F *hkin = new TH2F("hkin", "UKF kinematics (run_0116, all evts);#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 200, 0,
                         keMax);
   TH1F *hke = new TH1F("hke", "UKF kinetic energy;KE [MeV];tracks", 200, 0, keMax);
   TH1F *hc = new TH1F("hc", "#chi^{2}/ndf;#chi^{2}/ndf;tracks", 100, 0, 5);

   long nTot = 0, nGood = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (te->GetEntries() == 0)
         continue;
      auto *e = (AtTrackingEvent *)te->At(0);
      if (!e)
         continue;
      for (auto &ft : e->GetFittedTracks()) {
         if (!ft)
            continue;
         nTot++;
         auto &k = ft->GetKinematics();
         double ndf = ft->GetTrackMetadata()->GetNdf();
         double chi2 = ft->GetTrackMetadata()->GetChi2();
         double c2n = ndf > 0 ? chi2 / ndf : 1e9;
         double th = k.theta * TMath::RadToDeg();
         double ke = k.kineticEnergy;
         hc->Fill(c2n);
         if (ke > 0 && ke < 1000 && c2n < chi2Cut) {
            hkin->Fill(th, ke);
            hke->Fill(ke);
            nGood++;
         }
      }
   }

   TCanvas *c = new TCanvas("c", "ukf", 1500, 500);
   c->Divide(3, 1);
   c->cd(1);
   c->cd(1)->SetLogz();
   hkin->Draw("colz");
   c->cd(2);
   hke->SetFillColor(kAzure - 9);
   hke->Draw();
   c->cd(3);
   c->cd(3)->SetLogy();
   hc->SetFillColor(kOrange - 3);
   hc->Draw();
   TString png = fileName + "_ukf_results.png";
   c->SaveAs(png);
   printf("Fitted tracks: %ld total, %ld with KE>0 & chi2/ndf<%.1f\n", nTot, nGood, chi2Cut);
   printf("Saved %s\n", png.Data());
}
