/// @file pid_stats_a1975.C
/// @brief High-statistics PID plane across multiple runs, OLD AtPIDEstimator (the
///        charge/arclen recipe that made run_0116_pid_plane.png with two visible
///        bands) side-by-side with the faithful AtSpyralPID. Settles whether the
///        second band survives at higher statistics and under the faithful estimator.
///
///   root -b -q 'pid_stats_a1975.C("run_0116,run_0117,run_0118")'

#include <vector>

void pid_stats_a1975(TString runsCSV = "run_0116,run_0117,run_0118", double bField = 2.85)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   AtTools::AtPIDEstimator oldEst(bField, 152.0);
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH2F *hOld = new TH2F("hOld", "OLD recipe (charge/arclen);dEdx [counts];B#rho [T m]", 350, 0, 4000, 350, 0, 1.0);
   TH2F *hSpy = new TH2F("hSpy", "Faithful Spyral port;dEdx [counts];B#rho [T m]", 350, 0, 4000, 350, 0, 1.0);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nOld = 0, nSpy = 0, nTrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = run + "_reco.root";
      if (gSystem->AccessPathName(fn)) {
         printf("skip (missing) %s\n", fn.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      TClonesArray *peArr = nullptr;
      t->SetBranchAddress("AtPatternEvent", &peArr);
      Long64_t nE = t->GetEntries();
      for (Long64_t i = 0; i < nE; ++i) {
         t->GetEntry(i);
         if (peArr->GetEntries() == 0)
            continue;
         auto *pe = (AtPatternEvent *)peArr->At(0);
         if (!pe)
            continue;
         for (auto &trk : pe->GetTrackCand()) {
            ++nTrk;
            AtTrack &tr = const_cast<AtTrack &>(trk);
            auto ro = oldEst.Estimate(tr);
            if (ro.valid) {
               hOld->Fill(ro.dEdx, ro.brho);
               ++nOld;
            }
            auto rs = spy.Estimate(tr);
            if (rs.valid) {
               hSpy->Fill(rs.dEdx, rs.brho);
               ++nSpy;
            }
         }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nruns=%d  tracks=%ld  old-valid=%ld  spyral-valid=%ld\n", runs->GetEntries(), nTrk, nOld, nSpy);

   TCanvas *c = new TCanvas("c", "pidstats", 1500, 650);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz();
   hOld->Draw("colz");
   c->cd(2);
   gPad->SetLogz();
   hSpy->Draw("colz");
   c->SaveAs("pid_stats_compare.png");
   printf("saved pid_stats_compare.png\n");
}
