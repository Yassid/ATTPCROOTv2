/// @file pid_ic_a1975.C
/// @brief End-to-end IC-gated PID for a1975 16C+p, combining multiple runs.
///
/// For each run it index-matches the two DAQ streams 1:1 (events are
/// number-ordered):
///   <run>_reco.root  : AtPatternEvent  (TPC tracks, from AtHDFUnpacker 'get')
///   <run>_FRIB.root   : AtRawEvent      (8 generic traces, from frib/evt)
/// The ion chamber is generic trace[0]; ic_amplitude = baseline-subtracted peak
/// in the IC time window. Tracks are kept only if ic_amplitude is inside
/// [icMin, icMax] (the 16C beam peak), removing the contaminant-beam background.
///
/// Computes Spyral PID observables (AtTools::AtPIDEstimator) and draws the IC
/// spectrum + ungated vs IC-gated PID plane. Observables are cached to
/// <cacheTag>_pidobs.root for instant re-plotting / gate tuning (replot=true).
///
/// Run:
///   root -b -q 'pid_ic_a1975.C("run_0116,run_0117,run_0118")'
///   root -b -q 'pid_ic_a1975.C("","",950,1350,4000,1.0,15,true)'  // fast re-gate from cache

void pid_ic_a1975(TString runs = "run_0116,run_0117,run_0118", TString cacheTag = "combined", Double_t icMin = 950,
                  Double_t icMax = 1350, Double_t dedxMax = 4000, Double_t brMax = 1.0, Int_t minClusters = 15,
                  Bool_t replot = false, Int_t icTbLo = 1000, Int_t icTbHi = 1350, Double_t bField = 2.85,
                  Double_t maxVertexR = 40.0, Double_t polarMin = 10.0, Double_t polarMax = 170.0)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString cacheFile = cacheTag + "_pidobs.root";

   if (!replot) {
      AtTools::AtPIDEstimator estimator(bField, 152.0);
      TFile *fo = new TFile(cacheFile, "RECREATE");
      TNtuple *nt = new TNtuple("pidobs", "PID + IC", "dedx:brho:ncl:ic:polar:vtxr:vtxz:elm:elt");

      TObjArray *toks = runs.Tokenize(",");
      for (int ir = 0; ir < toks->GetEntries(); ++ir) {
         TString run = ((TObjString *)toks->At(ir))->GetString().Strip(TString::kBoth);
         TString recoF = run + "_reco.root", fribF = run + "_FRIB.root";
         if (gSystem->AccessPathName(recoF) || gSystem->AccessPathName(fribF)) {
            printf("\033[1;33mskip %s (missing reco or FRIB)\033[0m\n", run.Data());
            continue;
         }
         TFile *fr = TFile::Open(recoF);
         TTree *tr = (TTree *)fr->Get("cbmsim");
         tr->SetBranchStatus("*", 0);
         tr->SetBranchStatus("AtPatternEvent*", 1);
         TClonesArray *pe = nullptr;
         tr->SetBranchAddress("AtPatternEvent", &pe);

         TFile *ff = TFile::Open(fribF);
         TTree *tf = (TTree *)ff->Get("cbmsim");
         TClonesArray *re = nullptr;
         tf->SetBranchAddress("AtRawEvent", &re);

         Long64_t nR = tr->GetEntries(), nF = tf->GetEntries();
         Long64_t N = std::min(nR, nF);
         printf("%s: reco %lld evts, FRIB %lld evts -> matching %lld\n", run.Data(), nR, nF, N);
         long filled = 0;
         for (Long64_t i = 0; i < N; ++i) {
            // IC amplitude from the FRIB stream (trace 0, baseline-sub peak in window)
            tf->GetEntry(i);
            double ic = -1;
            if (re->GetEntries() > 0) {
               auto *raw = (AtRawEvent *)re->At(0);
               if (raw && !raw->GetGenTraces().empty()) {
                  auto &adc = raw->GetGenTraces()[0]->GetADC();
                  double mx = -1e9;
                  for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
                     mx = std::max(mx, adc[b]);
                  ic = mx;
               }
            }
            // TPC tracks
            tr->GetEntry(i);
            if (pe->GetEntries() == 0)
               continue;
            auto *patt = (AtPatternEvent *)pe->At(0);
            if (!patt)
               continue;
            for (auto &track : patt->GetTrackCand()) {
               auto r = estimator.Estimate(const_cast<AtTrack &>(track));
               if (!r.valid)
                  continue;
               double vtxr = std::sqrt(r.vertex.X() * r.vertex.X() + r.vertex.Y() * r.vertex.Y());
               float row[9] = {(float)r.dEdx,    (float)r.brho,      (float)r.nClusters,
                               (float)ic,        (float)(r.polar * TMath::RadToDeg()), (float)vtxr,
                               (float)r.vertex.Z(), (float)r.elossMean, (float)r.elossTrunc};
               nt->Fill(row);
               ++filled;
            }
         }
         printf("  filled %ld tracks\n", filled);
         fr->Close();
         ff->Close();
      }
      delete toks;
      fo->cd(); // opening/closing the reco+FRIB files moved gDirectory; restore the cache file
      nt->Write();
      fo->Close();
      printf("cached -> %s\n", cacheFile.Data());
   }

   // ---- plot from cache ----
   TFile *fo = TFile::Open(cacheFile);
   TNtuple *nt = (TNtuple *)fo->Get("pidobs");
   float dedx, brho, ncl, ic, polar, vtxr, vtxz;
   nt->SetBranchAddress("dedx", &dedx);
   nt->SetBranchAddress("brho", &brho);
   nt->SetBranchAddress("ncl", &ncl);
   nt->SetBranchAddress("ic", &ic);
   nt->SetBranchAddress("polar", &polar);
   nt->SetBranchAddress("vtxr", &vtxr);
   nt->SetBranchAddress("vtxz", &vtxz);

   TH1F *hic = new TH1F("hic", "IC amplitude (quality cut);IC amp [ADC];tracks", 250, 0, 2500);
   TH2F *hall = new TH2F("hall", "PID quality-cut, no IC gate;dEdx [counts];B#rho [T m]", 400, 0, dedxMax, 400, 0, brMax);
   TH2F *hgate = new TH2F("hgate", "PID quality + IC gate (16C);dEdx [counts];B#rho [T m]", 400, 0, dedxMax, 400, 0,
                          brMax);
   TH2F *hkin = new TH2F("hkin", "Kinematics: quality + IC gate;#theta [deg];B#rho [T m]", 180, 0, 180, 400, 0, brMax);
   long all = 0, gated = 0, raw = 0;
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      ++raw;
      // quality cuts: enough clusters, vertex near beam axis, polar window
      if (ncl < minClusters || vtxr > maxVertexR || polar < polarMin || polar > polarMax)
         continue;
      hic->Fill(ic);
      hall->Fill(dedx, brho);
      ++all;
      if (ic >= icMin && ic <= icMax) {
         hgate->Fill(dedx, brho);
         hkin->Fill(polar, brho);
         ++gated;
      }
   }
   printf("\033[1;32mtracks: %ld raw -> %ld after quality(ncl>=%d, vtxR<%.0f, polar[%.0f,%.0f]) -> %ld IC-gated "
          "[%.0f-%.0f]\033[0m\n",
          raw, all, minClusters, maxVertexR, polarMin, polarMax, gated, icMin, icMax);

   TCanvas *c = new TCanvas("c", "pid_ic", 1500, 1100);
   c->Divide(2, 2);
   c->cd(1);
   hic->SetFillColor(kAzure - 9);
   hic->Draw();
   double ymax = hic->GetMaximum();
   for (double x : {icMin, icMax}) {
      auto *l = new TLine(x, 0, x, ymax);
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
   }
   c->cd(2);
   gPad->SetLogz();
   hall->Draw("colz");
   c->cd(3);
   gPad->SetLogz();
   hgate->Draw("colz");
   c->cd(4);
   gPad->SetLogz();
   hkin->Draw("colz");
   TString png = cacheTag + "_pid_ic.png";
   c->SaveAs(png);
   printf("Saved %s\n", png.Data());
}
