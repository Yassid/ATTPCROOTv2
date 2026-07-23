/// @file ic_Be12.C
/// @brief ATTPCROOT Ion Chamber spectrum for a1954 12Be, from <run>_FRIB.root.
///        IC = generic trace[0]; ic_amplitude = baseline-subtracted peak in the IC
///        time window [icTbLo, icTbHi] (peak sits at TB~1132 for a1954). This is the
///        ATTPCROOT-side beam-ID used to gate the 12Be beam species.
///
///   root -b -q 'ic_Be12.C("run_0142")'
///   root -b -q 'ic_Be12.C("run_0147,run_0148,run_0149,run_0150")'
void ic_Be12(TString runsCSV = "run_0142", TString inDir = "/home/yassid/a1954_Be12_reco/", Int_t icTbLo = 1000,
             Int_t icTbHi = 1350, TString outTag = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954_Be12/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   TH1D *hAmp = new TH1D("hAmp", "ATTPCROOT IC amplitude (12Be beam ID);ic_amplitude [ADC];events", 500, 0, 4096);
   TH2D *hAmpTb = new TH2D("hAmpTb", "IC amplitude vs peak TB;peak TB;ic_amplitude", 200, 900, 1400, 300, 0, 4096);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nEvt = 0, nIC = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(ff)) {
         printf("skip %s (no %s)\n", run.Data(), ff.Data());
         continue;
      }
      TFile *f = TFile::Open(ff);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *ra = nullptr;
      t->SetBranchAddress("AtRawEvent", &ra);
      Long64_t N = t->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         ++nEvt;
         if (ra->GetEntries() == 0)
            continue;
         auto *raw = (AtRawEvent *)ra->At(0);
         if (!raw || raw->GetGenTraces().empty())
            continue;
         auto &adc = raw->GetGenTraces()[0]->GetADC();
         double mx = -1e9;
         int mb = 0;
         for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
            if (adc[b] > mx) {
               mx = adc[b];
               mb = b;
            }
         if (mx > 0) {
            hAmp->Fill(mx);
            hAmpTb->Fill(mb, mx);
            ++nIC;
         }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nATTPCROOT IC: events=%ld  with IC signal=%ld\n", nEvt, nIC);

   TCanvas *c = new TCanvas("c", "ic", 1200, 500);
   c->Divide(2, 1);
   c->cd(1);
   c->cd(1)->SetLogy();
   hAmp->Draw();
   c->cd(2);
   hAmpTb->Draw("colz");
   TString png = plotDir + "ic_Be12_attpc" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
