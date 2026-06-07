/// @file pid_compare_methods_a1975.C
/// @brief High-statistics side-by-side of the TWO PID methods produced by AtPIDTask.
///
/// Reads the AtPIDEvent branch from <run>_pid.root for several runs and draws the
/// classic AtPIDEstimator plane (charge/arclength) next to the Spyral-style
/// AtSpyralPID plane (first-arc + spline), both B#rho vs dEdx. Demonstrates, at full
/// statistics straight from the pipeline output, that the classic method shows the
/// spurious 2nd band while the Spyral-style one collapses to a single proton band.
///
/// Run: root -b -q 'pid_compare_methods_a1975.C("run_0116,run_0117,run_0118")'

void pid_compare_methods_a1975(TString runsCSV = "run_0116,run_0117,run_0118", TString inDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TH2F *hC =
      new TH2F("hC", "classic AtPIDEstimator (charge/arclen);#sqrt{dEdx};B#rho [T m]", 300, 0, 30, 400, 0, 2.0);
   TH2F *hS =
      new TH2F("hS", "Spyral-style AtSpyralPID (first-arc+spline);#sqrt{dEdx};B#rho [T m]", 300, 0, 40, 400, 0, 2.0);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nC = 0, nS = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = inDir + run + "_pid.root";
      if (gSystem->AccessPathName(fn)) {
         printf("skip (missing) %s\n", fn.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &r : ev->GetClassic())
            if (r.valid) {
               hC->Fill(r.sqrtdEdx, r.brho);
               ++nC;
            }
         for (auto &r : ev->GetSpyral())
            if (r.valid) {
               hS->Fill(r.sqrtdEdx, r.brho);
               ++nS;
            }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nclassic valid=%ld   spyral valid=%ld\n", nC, nS);

   TCanvas *c = new TCanvas("c", "pidmethods", 1500, 640);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz();
   hC->Draw("colz");
   c->cd(2);
   gPad->SetLogz();
   hS->Draw("colz");
   c->SaveAs("pid_methods_compare.png");
   printf("saved pid_methods_compare.png\n");
}
