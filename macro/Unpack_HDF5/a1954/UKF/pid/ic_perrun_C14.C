/// @file ic_perrun_C14.C
/// @brief Per-run ATTPCROOT IC amplitude spectra for the a1954 14C runs.
///        The 14C runs span 0011-0069 and were taken with different beam tunes, so the
///        IC gate has to be checked run by run before it is frozen: this overlays the
///        normalised IC amplitude of every run and prints the peak position / width.
///        IC = generic trace[0]; amplitude = max ADC in [icTbLo,icTbHi] (pulse sits ~TB1150).
///
///   root -b -q 'pid/ic_perrun_C14.C("run_0033,run_0055,...")'
void ic_perrun_C14(TString runsCSV = "run_0033,run_0055,run_0056,run_0057,run_0058",
                   TString inDir = "/home/yassid/a1954_C14_reco_hdb_slim/", Int_t icTbLo = 1050,
                   Int_t icTbHi = 1250, TString outTag = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   TObjArray *runs = runsCSV.Tokenize(",");
   std::vector<TH1D *> hists;
   std::vector<TString> names;
   const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};

   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(ff)) {
         printf("skip %s (no FRIB)\n", run.Data());
         continue;
      }
      TFile *f = TFile::Open(ff);
      TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t || t->GetEntries() == 0) {
         printf("skip %s (empty)\n", run.Data());
         if (f) f->Close();
         continue;
      }
      TClonesArray *ra = nullptr;
      t->SetBranchAddress("AtRawEvent", &ra);
      auto *h = new TH1D(Form("h_%s", run.Data()), run.Data(), 210, 0, 2100);
      h->SetDirectory(nullptr); // survive f->Close()
      Long64_t N = t->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (ra->GetEntries() == 0) continue;
         auto *raw = (AtRawEvent *)ra->At(0);
         if (!raw || raw->GetGenTraces().empty()) continue;
         auto &adc = raw->GetGenTraces()[0]->GetADC();
         double mx = -1e9;
         for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b) mx = std::max(mx, (double)adc[b]);
         if (mx > 0) h->Fill(mx);
      }
      f->Close();
      // peak position above 300 ADC (skip the no-signal pile at 0)
      int bmax = 0; double ymax = 0;
      for (int b = h->FindBin(300); b <= h->FindBin(2000); ++b)
         if (h->GetBinContent(b) > ymax) { ymax = h->GetBinContent(b); bmax = b; }
      printf("%-9s  entries=%8.0f  peak=%6.0f ADC (%5.0f counts)  frac>300ADC=%.3f\n", run.Data(), h->GetEntries(),
             h->GetBinCenter(bmax), ymax, h->Integral(h->FindBin(300), h->GetNbinsX()) / std::max(1.0, h->GetEntries()));
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
      h->SetLineColor(colors[hists.size() % 8]);
      h->SetLineWidth(2);
      hists.push_back(h);
      names.push_back(run);
   }
   if (hists.empty()) { printf("nothing to draw\n"); return; }

   TCanvas *c = new TCanvas("c", "ic per run", 1100, 700);
   c->SetLogy();
   auto *leg = new TLegend(0.62, 0.55, 0.89, 0.89);
   double ymx = 0;
   for (auto *h : hists) ymx = std::max(ymx, h->GetMaximum());
   hists[0]->SetTitle("a1954 14C: IC amplitude per run (area-normalised);ic_amplitude [ADC];fraction/bin");
   hists[0]->GetYaxis()->SetRangeUser(1e-5, ymx * 2);
   for (size_t i = 0; i < hists.size(); ++i) {
      hists[i]->Draw(i ? "hist same" : "hist");
      leg->AddEntry(hists[i], names[i], "l");
   }
   leg->Draw();
   TString png = plotDir + "ic_perrun_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
