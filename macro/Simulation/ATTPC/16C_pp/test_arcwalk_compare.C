/// @file test_arcwalk_compare.C
/// @brief Compare arc-walk vs Smooth3D PRA output.
///
/// Run after:
///   root -b -q 'test_arcwalk.C(500, false)'
///   root -b -q 'test_arcwalk.C(500, true)'
///
/// Then: root -b -q test_arcwalk_compare.C

void test_arcwalk_compare()
{
   TString inOutDir = "./data/";
   auto fSmooth = TFile::Open(inOutDir + "output_smooth3d.root");
   auto fArcWalk = TFile::Open(inOutDir + "output_arcwalk.root");

   if (!fSmooth || !fArcWalk) {
      std::cerr << "Run test_arcwalk.C with both modes first." << std::endl;
      return;
   }

   auto tSmooth = (TTree *)fSmooth->Get("cbmsim");
   auto tArcWalk = (TTree *)fArcWalk->Get("cbmsim");

   TClonesArray *patSmooth = nullptr, *patArcWalk = nullptr;
   tSmooth->SetBranchAddress("AtPatternEvent", &patSmooth);
   tArcWalk->SetBranchAddress("AtPatternEvent", &patArcWalk);

   int nEntries = std::min((int)tSmooth->GetEntries(), (int)tArcWalk->GetEntries());

   int nSmoothTracks = 0, nArcWalkTracks = 0;
   int nSmoothEvents = 0, nArcWalkEvents = 0;
   double sumCl_s = 0, sumCl_a = 0;
   double sumHits_s = 0, sumHits_a = 0;
   int nTr_s = 0, nTr_a = 0;

   TH1F *hClSmooth = new TH1F("hClSmooth", "Clusters/track;N clusters;Tracks", 60, 0, 60);
   TH1F *hClArcWalk = new TH1F("hClArcWalk", "Clusters/track;N clusters;Tracks", 60, 0, 60);
   TH1F *hHitsSmooth = new TH1F("hHitsSmooth", "Hits/track;N hits;Tracks", 100, 0, 500);
   TH1F *hHitsArcWalk = new TH1F("hHitsArcWalk", "Hits/track;N hits;Tracks", 100, 0, 500);

   for (int i = 0; i < nEntries; i++) {
      tSmooth->GetEntry(i);
      tArcWalk->GetEntry(i);

      if (patSmooth && patSmooth->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)patSmooth->At(0);
         auto &tracks = pe->GetTrackCand();
         if (!tracks.empty()) {
            nSmoothEvents++;
            for (auto &tr : tracks) {
               nSmoothTracks++;
               int ncl = tr.GetHitClusterArray()->size();
               int nhits = tr.GetHitArray().size();
               sumCl_s += ncl;
               sumHits_s += nhits;
               nTr_s++;
               hClSmooth->Fill(ncl);
               hHitsSmooth->Fill(nhits);
            }
         }
      }

      if (patArcWalk && patArcWalk->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)patArcWalk->At(0);
         auto &tracks = pe->GetTrackCand();
         if (!tracks.empty()) {
            nArcWalkEvents++;
            for (auto &tr : tracks) {
               nArcWalkTracks++;
               int ncl = tr.GetHitClusterArray()->size();
               int nhits = tr.GetHitArray().size();
               sumCl_a += ncl;
               sumHits_a += nhits;
               nTr_a++;
               hClArcWalk->Fill(ncl);
               hHitsArcWalk->Fill(nhits);
            }
         }
      }
   }

   std::cout << "\n========================================" << std::endl;
   std::cout << "      Smooth3D vs Arc-Walk Clustering" << std::endl;
   std::cout << "========================================" << std::endl;
   std::cout << "Events:              " << nEntries << std::endl;
   printf("                     %-12s%-12s\n", "Smooth3D", "Arc-Walk");
   printf("Events w/ tracks:    %-12d%-12d\n", nSmoothEvents, nArcWalkEvents);
   printf("Total tracks:        %-12d%-12d\n", nSmoothTracks, nArcWalkTracks);
   printf("Avg clusters/track:  %-12.1f%-12.1f\n", nTr_s > 0 ? sumCl_s / nTr_s : 0., nTr_a > 0 ? sumCl_a / nTr_a : 0.);
   printf("Avg hits/track:      %-12.1f%-12.1f\n", nTr_s > 0 ? sumHits_s / nTr_s : 0., nTr_a > 0 ? sumHits_a / nTr_a : 0.);
   std::cout << "========================================" << std::endl;

   TCanvas *c = new TCanvas("cComp", "Clustering Comparison", 1200, 500);
   c->Divide(2, 1);

   c->cd(1);
   hClSmooth->SetLineColor(kBlue);
   hClSmooth->SetLineWidth(2);
   hClArcWalk->SetLineColor(kRed);
   hClArcWalk->SetLineWidth(2);
   hClSmooth->Draw("hist");
   hClArcWalk->Draw("hist same");
   auto leg1 = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg1->AddEntry(hClSmooth, "Smooth3D", "l");
   leg1->AddEntry(hClArcWalk, "Arc-walk", "l");
   leg1->Draw();

   c->cd(2);
   hHitsSmooth->SetLineColor(kBlue);
   hHitsSmooth->SetLineWidth(2);
   hHitsArcWalk->SetLineColor(kRed);
   hHitsArcWalk->SetLineWidth(2);
   hHitsSmooth->Draw("hist");
   hHitsArcWalk->Draw("hist same");
   auto leg2 = new TLegend(0.55, 0.7, 0.88, 0.88);
   leg2->AddEntry(hHitsSmooth, "Smooth3D", "l");
   leg2->AddEntry(hHitsArcWalk, "Arc-walk", "l");
   leg2->Draw();

   c->SaveAs("data/arcwalk_comparison.png");
   std::cout << "Saved: data/arcwalk_comparison.png" << std::endl;

   fSmooth->Close();
   fArcWalk->Close();
}
