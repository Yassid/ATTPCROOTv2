/// @file perf_stage1.C
/// @brief Stage 1: PRA cluster quality -> GenFit fit outcome, and the key
/// population — GOOD clusterized tracks that GenFit drops/garbages (and tracks
/// UKF fits fine but GenFit drops). Reads pd/perf/perf_master.root.
///
/// GenFit outcome per PRA candidate:
///   DROP    : no genfit fitted track for this trackID
///   GARBAGE : fitted but KE<=0 or KE>1000 (bad momentum)
///   BADCHI2 : fitted, KE-ok, chi2/ndf > chi2Cut
///   GOOD    : fitted, KE-ok, chi2-ok
///
///   root -b -q 'pd/perf/perf_stage1.C'

void perf_stage1(double chi2Cut = 5.0, double goodClus = 15.0, TString in = "pd/perf/perf_master.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(in);
   TNtuple *m = (TNtuple *)f->Get("m");
   if (!m) { printf("no master ntuple\n"); return; }
   float nclus, gffit, gfke, gfchi2, ukfit, ukke, ukchi2, radius, gtheta, spyvalid;
   m->SetBranchAddress("nclus", &nclus);
   m->SetBranchAddress("gffit", &gffit);
   m->SetBranchAddress("gfke", &gfke);
   m->SetBranchAddress("gfchi2", &gfchi2);
   m->SetBranchAddress("ukfit", &ukfit);
   m->SetBranchAddress("ukke", &ukke);
   m->SetBranchAddress("radius", &radius);
   m->SetBranchAddress("gtheta", &gtheta);
   m->SetBranchAddress("spyvalid", &spyvalid);

   auto gfOutcome = [&]() -> int { // 0 drop 1 garbage 2 badchi2 3 good
      if (gffit < 0.5) return 0;
      if (gfke <= 0 || gfke > 1000) return 1;
      if (gfchi2 > chi2Cut) return 2;
      return 3;
   };
   auto ukGood = [&]() { return ukfit > 0.5 && ukke > 0 && ukke < 1000; };

   // outcome vs nclus
   const char *onm[4] = {"DROP", "GARBAGE", "BADCHI2", "GOOD"};
   long tot = 0, oc[4] = {0}, ocGood[4] = {0}; // ocGood = among good-clusterized tracks
   long ukg_gfdrop = 0, ukg_gfgarb = 0, ukg_gfgood = 0, ukg = 0; // UKF-good cross genfit
   long ukg_gfbad_goodclus = 0;                                   // UKF-good, genfit drop/garbage, AND good clusters
   TH1F *hAll = new TH1F("hAll", "nClusters;nClusters;tracks", 60, 0, 60);
   TH1F *hGood = new TH1F("hGood", "", 60, 0, 60);
   TH1F *hDrop = new TH1F("hDrop", "", 60, 0, 60);
   TH1F *hUkgGfbad = new TH1F("hUkgGfbad", "", 60, 0, 60);
   hAll->SetDirectory(nullptr); hGood->SetDirectory(nullptr); hDrop->SetDirectory(nullptr); hUkgGfbad->SetDirectory(nullptr);
   TProfile *pGoodVsClus = new TProfile("pGoodVsClus", "GenFit GOOD-fit fraction vs nClusters;nClusters;P(GOOD)", 30, 0, 60, 0, 1);
   pGoodVsClus->SetDirectory(nullptr);

   for (Long64_t i = 0; i < m->GetEntries(); ++i) {
      m->GetEntry(i);
      int o = gfOutcome();
      ++tot; ++oc[o];
      hAll->Fill(nclus);
      if (o == 3) hGood->Fill(nclus);
      if (o == 0) hDrop->Fill(nclus);
      pGoodVsClus->Fill(nclus, o == 3 ? 1.0 : 0.0);
      bool goodclus = nclus >= goodClus;
      if (goodclus) ++ocGood[o];
      if (ukGood()) {
         ++ukg;
         if (o == 0) ++ukg_gfdrop;
         else if (o == 1) ++ukg_gfgarb;
         else if (o == 3) ++ukg_gfgood;
         if ((o == 0 || o == 1) && goodclus) { ++ukg_gfbad_goodclus; hUkgGfbad->Fill(nclus); }
      }
   }

   printf("\n=== STAGE 1: GenFit fit outcome vs PRA cluster quality (%ld candidate tracks) ===\n", tot);
   printf("  ALL tracks:        ");
   for (int o = 0; o < 4; ++o) printf("%s %.1f%%  ", onm[o], 100.0 * oc[o] / tot);
   printf("\n  nClus>=%.0f (good):  ", goodClus);
   long totGood = ocGood[0] + ocGood[1] + ocGood[2] + ocGood[3];
   for (int o = 0; o < 4; ++o) printf("%s %.1f%%  ", onm[o], 100.0 * ocGood[o] / std::max(1L, totGood));
   printf("  (n=%ld)\n", totGood);
   printf("\n=== UKF vs GenFit (the key cross-check) ===\n");
   printf("  UKF-good tracks: %ld\n", ukg);
   printf("    of which GenFit GOOD   : %ld (%.1f%%)\n", ukg_gfgood, 100.0 * ukg_gfgood / std::max(1L, ukg));
   printf("    of which GenFit DROP   : %ld (%.1f%%)\n", ukg_gfdrop, 100.0 * ukg_gfdrop / std::max(1L, ukg));
   printf("    of which GenFit GARBAGE: %ld (%.1f%%)\n", ukg_gfgarb, 100.0 * ukg_gfgarb / std::max(1L, ukg));
   printf("  >> UKF-good + GOOD clusters (>=%.0f) but GenFit DROP/GARBAGE: %ld <<\n", goodClus, ukg_gfbad_goodclus);

   TCanvas *c = new TCanvas("c", "stage1", 1300, 520);
   c->Divide(2, 1);
   c->cd(1);
   pGoodVsClus->SetLineColor(kBlue + 1); pGoodVsClus->SetLineWidth(2); pGoodVsClus->SetMinimum(0); pGoodVsClus->SetMaximum(1);
   pGoodVsClus->Draw();
   c->cd(2);
   hAll->SetLineColor(kBlack); hAll->SetLineWidth(2); hAll->SetTitle("nClusters by genfit outcome;nClusters;tracks"); hAll->Draw("hist");
   hGood->SetLineColor(kGreen + 2); hGood->SetLineWidth(2); hGood->Draw("hist same");
   hDrop->SetLineColor(kRed + 1); hDrop->SetLineWidth(2); hDrop->Draw("hist same");
   hUkgGfbad->SetLineColor(kMagenta + 1); hUkgGfbad->SetLineWidth(2); hUkgGfbad->SetLineStyle(2); hUkgGfbad->Draw("hist same");
   TLegend *lg = new TLegend(0.55, 0.7, 0.88, 0.88);
   lg->AddEntry(hAll, "all candidates", "l"); lg->AddEntry(hGood, "genfit GOOD", "l");
   lg->AddEntry(hDrop, "genfit DROP", "l"); lg->AddEntry(hUkgGfbad, "UKF-good, genfit bad", "l");
   lg->Draw();
   c->SaveAs("pd/perf/perf_stage1.png");
   printf("saved pd/perf/perf_stage1.png\n");
}
