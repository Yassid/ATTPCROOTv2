/// @file ic_gate_C14.C
/// @brief a1954 14C beam gate: the ion-chamber amplitude spectrum ALONE, with the accepted
///        window marked. This is the left panel of the old two-panel ic_C14.C figure
///        (pid/plots/ic_C14_attpc_win.png, shipped as figures/ic_gate_C14.png in the
///        write-up) with vertical lines at icLo / icHi.
///
///        IC = generic trace[0]; amplitude = peak in the IC time window [icTbLo, icTbHi].
///        The defaults are the PRODUCTION gate of pipeline/gate_events_C14.C
///        (amplitude 950-1350, TB 1050-1250) -- note the old figure was made with the
///        wider TB search window 1000-1350, so the two spectra are not byte-identical.
///        The pile-up (exactly-one-pulse) condition of the production gate is NOT applied
///        here: this figure shows the amplitude cut, which is what the lines mark.
///
///   root -b -q 'ic_gate_C14.C'                       // 14 production runs, ~2 GB each
///   root -b -q 'ic_gate_C14.C("run_0055")'           // one run
///   root -b -q 'ic_gate_C14.C("run_0055","",950,1350,1050,1250,"_r55")'
///
/// The filled histogram is cached beside the png so the plot can be restyled without
/// re-reading the raw files. The cache is only reused when the runs AND the time window
/// match; anything else re-scans.
void ic_gate_C14(TString runsCSV = "run_0055,run_0056,run_0057,run_0058,run_0059,run_0060,run_0061,"
                                   "run_0062,run_0063,run_0064,run_0065,run_0066,run_0068,run_0069",
                 TString inDir = "/mnt/f/a1954_C14_reco_hdb_slim/", double icLo = 950, double icHi = 1350,
                 Int_t icTbLo = 1050, Int_t icTbHi = 1250, TString outTag = "", Bool_t useCache = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);
   TString cacheFile = plotDir + "ic_gate_C14" + outTag + ".root";
   TString stamp = TString::Format("%s|tb%d-%d", runsCSV.Data(), icTbLo, icTbHi);

   TH1D *hAmp = nullptr;
   long nEvt = 0, nIC = 0;

   // ---- cache: reuse only if the INPUTS are the same ones ----------------------------
   if (useCache && !gSystem->AccessPathName(cacheFile)) {
      TFile *fc = TFile::Open(cacheFile);
      auto *st = fc ? (TNamed *)fc->Get("stamp") : nullptr;
      auto *hc = fc ? (TH1D *)fc->Get("hAmp") : nullptr;
      if (st && hc && stamp == st->GetTitle()) {
         hAmp = (TH1D *)hc->Clone("hAmp");
         hAmp->SetDirectory(nullptr);
         nIC = (long)hAmp->GetEntries();
         printf("cache HIT  %s  (%ld IC entries)\n", cacheFile.Data(), nIC);
      } else {
         printf("cache MISS %s  (stamp differs, re-scanning)\n  want: %s\n  have: %s\n", cacheFile.Data(),
                stamp.Data(), st ? st->GetTitle() : "(none)");
      }
      if (fc) fc->Close();
   }

   // ---- scan the raw files ------------------------------------------------------------
   if (!hAmp) {
      hAmp = new TH1D("hAmp", ";IC amplitude [ADC];events", 500, 0, 4096);
      hAmp->SetDirectory(nullptr);
      TObjArray *runs = runsCSV.Tokenize(",");
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
            for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
               if (adc[b] > mx)
                  mx = adc[b];
            if (mx > 0) {
               hAmp->Fill(mx);
               ++nIC;
            }
         }
         f->Close();
         printf("processed %s (%lld events)\n", run.Data(), N);
      }
      TFile *fo = TFile::Open(cacheFile, "RECREATE");
      hAmp->Write("hAmp");
      TNamed st("stamp", stamp.Data());
      st.Write("stamp");
      fo->Close();
      printf("cached %s\n", cacheFile.Data());
   }

   // ---- what the window keeps ---------------------------------------------------------
   double all = hAmp->Integral(1, hAmp->GetNbinsX());
   double in = hAmp->Integral(hAmp->FindBin(icLo), hAmp->FindBin(icHi) - 1);
   printf("\nIC amplitude: entries=%.0f  inside [%.0f,%.0f]=%.0f (%.1f %%)  events read=%ld\n", all, icLo, icHi, in,
          all > 0 ? 100. * in / all : 0., nEvt);

   // ---- the plot ----------------------------------------------------------------------
   TCanvas *c = new TCanvas("cIcGate", "IC beam gate", 700, 550);
   c->SetLogy();
   c->SetTicks(1, 1);
   hAmp->SetLineColor(kBlack);
   hAmp->SetLineWidth(2);
   hAmp->GetYaxis()->SetTitleOffset(1.25);
   // frame on the data: the 0-4096 ADC range is mostly empty (saturation sits near 2050)
   int last = hAmp->FindLastBinAbove(0.);
   if (last > 0)
      hAmp->GetXaxis()->SetRangeUser(0., hAmp->GetXaxis()->GetBinUpEdge(last) * 1.05);
   hAmp->Draw("hist");
   gPad->Update();

   double y1 = TMath::Power(10, gPad->GetUymin());
   double y2 = TMath::Power(10, gPad->GetUymax());
   for (double x : {icLo, icHi}) {
      auto *l = new TLine(x, y1, x, y2);
      l->SetLineColor(kRed + 1);
      l->SetLineStyle(2);
      l->SetLineWidth(2);
      l->Draw();
   }
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.033);
   tx->SetTextColor(kRed + 1);
   tx->DrawLatex(0.57, 0.87, TString::Format("accepted window %.0f#minus%.0f ADC", icLo, icHi));
   tx->SetTextColor(kBlack);
   tx->DrawLatex(0.57, 0.81, TString::Format("%.1f %% of IC signals", all > 0 ? 100. * in / all : 0.));

   TString png = plotDir + "ic_gate_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
