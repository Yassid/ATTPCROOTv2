/// @file backward_fit_eff.C
/// @brief How many BACKWARD tracks does genfit lose between PRA and the fit?
///        Compares PRA pattern tracks (reco.root, AtPatternEvent::GetTrackCand,
///        GetGeoTheta) to genfit-fitted tracks (genfitter_p.root, AtTrackingEvent::
///        GetFittedTracks), matched by trackID per event, binned in PRA geometric
///        theta. Answers "how many backward protons are we missing from the fit".
///
///   root -l -b -q 'backward_fit_eff.C("run_0016")'
///
/// Genfit processes events in order, so genfitter entry i == reco entry i. The
/// genfitter file is truncated at ~18k/41k events (crash) -> we compare only the
/// events genfit actually reached (its GetEntries()).

void backward_fit_eff(TString run = "run_0016", TString dir = "/mnt/f/a1975/reco_d2/", double chi2Cut = 5.0)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TFile *fr = TFile::Open(dir + run + "_reco.root");
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pat = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pat);

   TFile *fg = TFile::Open(dir + run + "_genfitter_p.root");
   TTree *tg = (TTree *)fg->Get("cbmsim");
   TClonesArray *trk = nullptr;
   tg->SetBranchAddress("AtTrackingEvent", &trk);

   Long64_t Ng = tg->GetEntries(); // events genfit actually reached
   printf("comparing %lld events (genfit reached) of run %s\n", Ng, run.Data());

   TH1F *hPat = new TH1F("hPat", "tracks vs PRA #theta_{geo};#theta_{geo} [deg];tracks", 36, 0, 180);
   TH1F *hAny = new TH1F("hAny", "fitted (any chi2)", 36, 0, 180);
   TH1F *hGood = new TH1F("hGood", "fitted (chi2/ndf<cut)", 36, 0, 180);
   hPat->SetLineColor(kBlack);
   hPat->SetLineWidth(2);
   hAny->SetLineColor(kAzure + 1);
   hAny->SetLineWidth(2);
   hGood->SetLineColor(kGreen + 2);
   hGood->SetLineWidth(2);

   long patBack = 0, anyBack = 0, goodBack = 0, patFwd = 0, anyFwd = 0, goodFwd = 0;
   long nMatchTot = 0, nPatTot = 0;

   for (Long64_t i = 0; i < Ng; ++i) {
      tr->GetEntry(i);
      tg->GetEntry(i);
      std::map<int, double> fitC2; // trackID -> chi2/ndf
      if (trk->GetEntries() > 0) {
         auto *ev = (AtTrackingEvent *)trk->At(0);
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft || !ft->GetTrackMetadata())
               continue;
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            fitC2[ft->GetTrackID()] = ndf > 0 ? chi2 / ndf : 1e9;
         }
      }
      if (pat->GetEntries() == 0)
         continue;
      auto *pe = (AtPatternEvent *)pat->At(0);
      for (auto &tk : pe->GetTrackCand()) {
         double th = tk.GetGeoTheta() * TMath::RadToDeg();
         if (th <= 0 || th > 180)
            continue; // unset/garbage geo angle
         ++nPatTot;
         bool back = th >= 90.0;
         hPat->Fill(th);
         if (back)
            ++patBack;
         else
            ++patFwd;
         auto it = fitC2.find(tk.GetTrackID());
         if (it != fitC2.end()) {
            ++nMatchTot;
            hAny->Fill(th);
            if (back)
               ++anyBack;
            else
               ++anyFwd;
            if (it->second < chi2Cut) {
               hGood->Fill(th);
               if (back)
                  ++goodBack;
               else
                  ++goodFwd;
            }
         }
      }
   }

   printf("\n================= FIT EFFICIENCY (PRA -> genfit) =================\n");
   printf("PRA pattern tracks with valid geo-theta : %ld  (trackID-matched to a fit: %ld)\n", nPatTot, nMatchTot);
   printf("\n  FORWARD  (theta_geo < 90):  PRA=%ld  fit(any)=%ld (%.1f%%)  fit(good)=%ld (%.1f%%)\n", patFwd, anyFwd,
          100.0 * anyFwd / std::max(1L, patFwd), goodFwd, 100.0 * goodFwd / std::max(1L, patFwd));
   printf("  BACKWARD (theta_geo >=90):  PRA=%ld  fit(any)=%ld (%.1f%%)  fit(good)=%ld (%.1f%%)\n", patBack, anyBack,
          100.0 * anyBack / std::max(1L, patBack), goodBack, 100.0 * goodBack / std::max(1L, patBack));
   printf("\n  ==> BACKWARD tracks MISSING from the (good) fit: %ld of %ld  (%.1f%%)\n", patBack - goodBack, patBack,
          100.0 * (patBack - goodBack) / std::max(1L, patBack));
   printf("==================================================================\n");

   // efficiency curve
   TH1F *hEff = (TH1F *)hGood->Clone("hEff");
   hEff->Divide(hPat);
   hEff->SetTitle("genfit good-fit efficiency vs PRA #theta_{geo};#theta_{geo} [deg];good fits / PRA tracks");
   hEff->SetLineColor(kRed);
   hEff->SetLineWidth(2);

   TCanvas *c = new TCanvas("c", "backward_fit_eff", 1300, 560);
   c->Divide(2, 1);
   c->cd(1);
   hPat->Draw("hist");
   hAny->Draw("hist same");
   hGood->Draw("hist same");
   TLegend *l = new TLegend(0.42, 0.72, 0.88, 0.88);
   l->AddEntry(hPat, "PRA pattern tracks", "l");
   l->AddEntry(hAny, "genfit fitted (any)", "l");
   l->AddEntry(hGood, Form("genfit fitted (chi2/ndf<%.0f)", chi2Cut), "l");
   l->Draw();
   TLine *l90 = new TLine(90, 0, 90, hPat->GetMaximum());
   l90->SetLineStyle(2);
   l90->SetLineColor(kGray + 1);
   l90->Draw();
   c->cd(2);
   hEff->SetMinimum(0);
   hEff->SetMaximum(1.05);
   hEff->Draw("hist");
   TLine *e90 = new TLine(90, 0, 90, 1.05);
   e90->SetLineStyle(2);
   e90->SetLineColor(kGray + 1);
   e90->Draw();
   c->SaveAs("plots/backward_fit_eff.png");
   printf("saved plots/backward_fit_eff.png\n");
}
