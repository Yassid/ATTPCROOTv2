/// @file cluster_quality_dp.C
/// @brief How good is the clustering for BACKWARD (d,p) tracks vs forward? Reads the
///        PRA tracks (AtPatternEvent in reco_d2/<run>_reco.root) — each carries its raw
///        hits AND its clusters — and measures, vs PRA GeoTheta: clusters/track,
///        raw-hits/track, hits-per-cluster, mean & max consecutive-cluster gap (3D mm),
///        track chord length, and clustered-coverage fraction. Backward = GeoTheta>=90
///        (the genfit-validated convention). No re-run needed.
///
///   root -l -b -q 'cluster_quality_dp.C("run_0016,run_0017,run_0018")'

void cluster_quality_dp(TString runsCSV = "run_0016", TString dir = "/mnt/f/a1975/reco_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   // profiles vs GeoTheta(deg)
   TProfile *pNcl = new TProfile("pNcl", "mean # clusters vs #theta_{geo};#theta_{geo} [deg];<# clusters>", 36, 0, 180);
   TProfile *pHpc = new TProfile("pHpc", "mean hits/cluster vs #theta_{geo};#theta_{geo} [deg];<hits/cluster>", 36, 0,
                                 180);
   TProfile *pGap = new TProfile("pGap", "mean cluster gap vs #theta_{geo};#theta_{geo} [deg];<gap> [mm]", 36, 0, 180);
   TProfile *pMaxGap =
      new TProfile("pMaxGap", "max cluster gap vs #theta_{geo};#theta_{geo} [deg];<max gap> [mm]", 36, 0, 180);
   // backward vs forward distributions of clusters/track
   TH1F *hClF = new TH1F("hClF", "clusters/track;# clusters;tracks (norm)", 60, 0, 60);
   TH1F *hClB = new TH1F("hClB", "", 60, 0, 60);
   hClF->SetLineColor(kAzure + 1);
   hClB->SetLineColor(kRed);
   hClF->SetLineWidth(2);
   hClB->SetLineWidth(2);

   long nF = 0, nB = 0;
   double sNclF = 0, sNclB = 0, sHpcF = 0, sHpcB = 0, sGapF = 0, sGapB = 0, sCovF = 0, sCovB = 0;

   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = dir + run + "_reco.root";
      if (gSystem->AccessPathName(fn)) {
         printf("skip (missing) %s\n", fn.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pat = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pat);
      Long64_t N = t->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (pat->GetEntries() == 0)
            continue;
         auto *pe = (AtPatternEvent *)pat->At(0);
         for (auto &tk : pe->GetTrackCand()) {
            double th = tk.GetGeoTheta() * TMath::RadToDeg();
            if (th <= 0 || th > 180)
               continue;
            auto *cl = tk.GetHitClusterArray();
            int ncl = cl ? cl->size() : 0;
            int nhit = tk.GetHitArray().size();
            if (ncl < 2)
               continue;
            // consecutive-cluster 3D gaps (clusters are PRA-ordered along the track)
            double sumGap = 0, maxGap = 0, chord = 0;
            for (int k = 1; k < ncl; ++k) {
               auto a = cl->at(k - 1).GetPosition();
               auto b = cl->at(k).GetPosition();
               double d = std::sqrt(std::pow(a.X() - b.X(), 2) + std::pow(a.Y() - b.Y(), 2) + std::pow(a.Z() - b.Z(), 2));
               sumGap += d;
               if (d > maxGap)
                  maxGap = d;
            }
            auto p0 = cl->front().GetPosition();
            auto pN = cl->back().GetPosition();
            chord = std::sqrt(std::pow(p0.X() - pN.X(), 2) + std::pow(p0.Y() - pN.Y(), 2) + std::pow(p0.Z() - pN.Z(), 2));
            double meanGap = sumGap / (ncl - 1);
            double hpc = nhit > 0 ? (double)nhit / ncl : 0;
            double coverage = chord > 0 ? sumGap / chord : 0; // path-along-clusters / chord (>=1; big = curl)

            pNcl->Fill(th, ncl);
            pHpc->Fill(th, hpc);
            pGap->Fill(th, meanGap);
            pMaxGap->Fill(th, maxGap);
            bool back = th >= 90;
            if (back) {
               hClB->Fill(ncl);
               nB++;
               sNclB += ncl;
               sHpcB += hpc;
               sGapB += meanGap;
               sCovB += coverage;
            } else {
               hClF->Fill(ncl);
               nF++;
               sNclF += ncl;
               sHpcF += hpc;
               sGapF += meanGap;
               sCovF += coverage;
            }
         }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }

   printf("\n=============== CLUSTERING QUALITY: forward vs backward (PRA tracks) ===============\n");
   printf("                  N tracks   <clusters>  <hits/clus>  <gap mm>  <pathlen/chord>\n");
   printf("  FORWARD (<90):  %8ld   %9.1f   %10.1f   %7.1f   %.2f\n", nF, sNclF / std::max(1L, nF),
          sHpcF / std::max(1L, nF), sGapF / std::max(1L, nF), sCovF / std::max(1L, nF));
   printf("  BACKWARD(>=90): %8ld   %9.1f   %10.1f   %7.1f   %.2f\n", nB, sNclB / std::max(1L, nB),
          sHpcB / std::max(1L, nB), sGapB / std::max(1L, nB), sCovB / std::max(1L, nB));
   printf("===================================================================================\n");

   TCanvas *c = new TCanvas("c", "cluster_quality", 1400, 900);
   c->Divide(2, 2);
   c->cd(1);
   if (hClB->Integral() > 0)
      hClB->Scale(1. / hClB->Integral());
   if (hClF->Integral() > 0)
      hClF->Scale(1. / hClF->Integral());
   hClB->SetTitle("clusters/track (norm);# clusters;");
   hClB->SetMaximum(1.2 * std::max(hClF->GetMaximum(), hClB->GetMaximum()));
   hClB->Draw("hist");
   hClF->Draw("hist same");
   TLegend *l = new TLegend(0.5, 0.72, 0.88, 0.88);
   l->AddEntry(hClF, "forward (<90)", "l");
   l->AddEntry(hClB, "backward (>=90)", "l");
   l->Draw();
   auto drawP = [&](TProfile *p, int pad) {
      c->cd(pad);
      p->SetLineColor(kBlack);
      p->SetLineWidth(2);
      p->Draw();
      TLine *m = new TLine(90, p->GetMinimum(), 90, p->GetMaximum());
      m->SetLineStyle(2);
      m->SetLineColor(kGray + 2);
      m->Draw();
   };
   drawP(pNcl, 2);
   drawP(pHpc, 3);
   drawP(pMaxGap, 4);
   c->SaveAs("plots/cluster_quality_dp.png");
   printf("saved plots/cluster_quality_dp.png\n");
}
