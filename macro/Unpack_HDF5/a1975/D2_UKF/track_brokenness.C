/// @file track_brokenness.C
/// @brief Test the "backward tracks are broken" hypothesis. For each PRA pattern
///        track matched to its genfit fit, record track-quality observables
///        (raw hits, clusters, fit Ndf, chi2/ndf) and compare FORWARD vs BACKWARD,
///        and ask: does good-fit fraction depend on cluster count? If backward
///        failures are driven by few-cluster/short ("broken") tracks, the good-fit
///        fraction should rise steeply with nClusters and backward tracks should be
///        cluster-starved relative to forward.
///
///   root -l -b -q 'track_brokenness.C("run_0016")'

void track_brokenness(TString run = "run_0016", TString dir = "/mnt/f/a1975/reco_d2/", double chi2Cut = 5.0)
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
   Long64_t Ng = tg->GetEntries();

   TH1F *hClF = new TH1F("hClF", "clusters/track;# clusters;tracks (norm)", 60, 0, 60);
   TH1F *hClB = new TH1F("hClB", "clusters/track", 60, 0, 60);
   TH1F *hHitF = new TH1F("hHitF", "raw hits/track;# hits;tracks (norm)", 100, 0, 400);
   TH1F *hHitB = new TH1F("hHitB", "raw hits/track", 100, 0, 400);
   hClF->SetLineColor(kAzure + 1);
   hClB->SetLineColor(kRed);
   hHitF->SetLineColor(kAzure + 1);
   hHitB->SetLineColor(kRed);
   for (auto h : {hClF, hClB, hHitF, hHitB})
      h->SetLineWidth(2);
   TProfile *pClus = new TProfile("pClus", "mean clusters vs #theta_{geo};#theta_{geo} [deg];<# clusters>", 36, 0, 180);
   TProfile *pChi2 = new TProfile("pChi2", "mean chi2/ndf vs #theta_{geo};#theta_{geo} [deg];<chi2/ndf>", 36, 0, 180);
   // good-fit fraction vs nClusters, forward & backward
   TH1F *hAllClF = new TH1F("hAllClF", "", 30, 0, 60), *hGoodClF = new TH1F("hGoodClF", "", 30, 0, 60);
   TH1F *hAllClB = new TH1F("hAllClB", "", 30, 0, 60), *hGoodClB = new TH1F("hGoodClB", "", 30, 0, 60);

   long nF = 0, nB = 0;
   double sClF = 0, sClB = 0, sHitF = 0, sHitB = 0;
   std::vector<int> clF, clB; // for medians

   for (Long64_t i = 0; i < Ng; ++i) {
      tr->GetEntry(i);
      tg->GetEntry(i);
      std::map<int, double> fitC2;
      if (trk->GetEntries() > 0) {
         auto *ev = (AtTrackingEvent *)trk->At(0);
         for (auto &ft : ev->GetFittedTracks())
            if (ft && ft->GetTrackMetadata()) {
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
            continue;
         int nhit = tk.GetHitArray().size();
         int nclus = tk.GetHitClusterArray() ? tk.GetHitClusterArray()->size() : 0;
         bool back = th >= 90.0;
         pClus->Fill(th, nclus);
         auto it = fitC2.find(tk.GetTrackID());
         bool matched = it != fitC2.end();
         bool good = matched && it->second < chi2Cut;
         if (matched)
            pChi2->Fill(th, std::min(it->second, 50.0));
         if (back) {
            hClB->Fill(nclus);
            hHitB->Fill(nhit);
            hAllClB->Fill(nclus);
            if (good)
               hGoodClB->Fill(nclus);
            nB++;
            sClB += nclus;
            sHitB += nhit;
            clB.push_back(nclus);
         } else {
            hClF->Fill(nclus);
            hHitF->Fill(nhit);
            hAllClF->Fill(nclus);
            if (good)
               hGoodClF->Fill(nclus);
            nF++;
            sClF += nclus;
            sHitF += nhit;
            clF.push_back(nclus);
         }
      }
   }
   std::sort(clF.begin(), clF.end());
   std::sort(clB.begin(), clB.end());
   auto med = [](std::vector<int> &v) { return v.empty() ? 0 : v[v.size() / 2]; };

   printf("\n=============== TRACK BROKENNESS: forward vs backward ===============\n");
   printf("                    N tracks   <clusters>  median clusters   <raw hits>\n");
   printf("  FORWARD (<90):   %8ld   %8.1f   %12d   %10.1f\n", nF, sClF / std::max(1L, nF), med(clF),
          sHitF / std::max(1L, nF));
   printf("  BACKWARD(>=90):  %8ld   %8.1f   %12d   %10.1f\n", nB, sClB / std::max(1L, nB), med(clB),
          sHitB / std::max(1L, nB));
   printf("\n  good-fit fraction vs cluster count:\n");
   printf("    nClus bin   FORWARD good%%   BACKWARD good%%\n");
   for (int b = 1; b <= 6; ++b) {
      double lo = (b - 1) * 10, hi = b * 10;
      double aF = hAllClF->Integral(hAllClF->FindBin(lo), hAllClF->FindBin(hi) - 1);
      double gF = hGoodClF->Integral(hGoodClF->FindBin(lo), hGoodClF->FindBin(hi) - 1);
      double aB = hAllClB->Integral(hAllClB->FindBin(lo), hAllClB->FindBin(hi) - 1);
      double gB = hGoodClB->Integral(hGoodClB->FindBin(lo), hGoodClB->FindBin(hi) - 1);
      printf("    %3.0f-%3.0f      %6.1f (%5.0f)    %6.1f (%5.0f)\n", lo, hi, aF > 0 ? 100 * gF / aF : 0, aF,
             aB > 0 ? 100 * gB / aB : 0, aB);
   }
   printf("====================================================================\n");

   TCanvas *c = new TCanvas("c", "brokenness", 1400, 900);
   c->Divide(2, 2);
   c->cd(1);
   if (hClB->Integral() > 0)
      hClB->Scale(1. / hClB->Integral());
   if (hClF->Integral() > 0)
      hClF->Scale(1. / hClF->Integral());
   hClB->SetMaximum(1.15 * std::max(hClF->GetMaximum(), hClB->GetMaximum()));
   hClB->Draw("hist");
   hClF->Draw("hist same");
   TLegend *l = new TLegend(0.5, 0.7, 0.88, 0.88);
   l->AddEntry(hClF, "forward (<90)", "l");
   l->AddEntry(hClB, "backward (>=90)", "l");
   l->Draw();
   c->cd(2);
   if (hHitB->Integral() > 0)
      hHitB->Scale(1. / hHitB->Integral());
   if (hHitF->Integral() > 0)
      hHitF->Scale(1. / hHitF->Integral());
   hHitB->SetMaximum(1.15 * std::max(hHitF->GetMaximum(), hHitB->GetMaximum()));
   hHitB->Draw("hist");
   hHitF->Draw("hist same");
   c->cd(3);
   pClus->SetLineColor(kBlack);
   pClus->SetLineWidth(2);
   pClus->Draw();
   TLine *m1 = new TLine(90, 0, 90, pClus->GetMaximum());
   m1->SetLineStyle(2);
   m1->Draw();
   c->cd(4);
   // good-fit fraction vs nClusters
   TH1F *effF = (TH1F *)hGoodClF->Clone("effF");
   effF->Divide(hAllClF);
   TH1F *effB = (TH1F *)hGoodClB->Clone("effB");
   effB->Divide(hAllClB);
   effF->SetTitle("good-fit fraction vs # clusters;# clusters;good fits / matched");
   effF->SetLineColor(kAzure + 1);
   effB->SetLineColor(kRed);
   effF->SetLineWidth(2);
   effB->SetLineWidth(2);
   effF->SetMaximum(1.05);
   effF->SetMinimum(0);
   effF->Draw("hist");
   effB->Draw("hist same");
   TLegend *l2 = new TLegend(0.5, 0.2, 0.88, 0.38);
   l2->AddEntry(effF, "forward", "l");
   l2->AddEntry(effB, "backward", "l");
   l2->Draw();
   c->SaveAs("plots/track_brokenness.png");
   printf("saved plots/track_brokenness.png\n");
}
