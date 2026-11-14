/// @file ic_gate_scan_C14.C
/// @brief a1954 14C beam gate, step 1 of 2: scan the raw AT-TPC files and write the ion-chamber
///        amplitude spectrum that ic_gate_C14.C draws. NEEDS ATTPCROOT (AtRawEvent) and the raw
///        <run>_FRIB.root files; the plotting step does not -- it runs on plain ROOT from the
///        histogram this writes, which the a1954_C14_absnorm repo vendors as data/ic_gate_C14.root.
///
///        IC = generic trace[0]; amplitude = peak in the IC time window [icTbLo, icTbHi].
///        The defaults are the PRODUCTION gate of pipeline/gate_events_C14.C (TB 1050-1250)
///        over the 14 production runs.
///
///   root -b -q 'ic_gate_scan_C14.C'                          // 14 production runs, ~2 GB each
///   root -b -q 'ic_gate_scan_C14.C("run_0055","",1050,1250,"_r55")'   // one run
///
/// Output: plots/ic_gate_C14<outTag>.root beside this macro, holding hAmp and a "stamp" TNamed
/// ("<runs>|tb<lo>-<hi>") that ic_gate_C14.C checks before it draws. To refresh the vendored copy,
/// copy that file to a1954_C14_absnorm/data/.
void ic_gate_scan_C14(TString runsCSV = "run_0055,run_0056,run_0057,run_0058,run_0059,run_0060,run_0061,"
                                        "run_0062,run_0063,run_0064,run_0065,run_0066,run_0068,run_0069",
                      TString inDir = "/mnt/f/a1954_C14_reco_hdb_slim/", Int_t icTbLo = 1050, Int_t icTbHi = 1250,
                      TString outTag = "")
{
   gSystem->Load("libAtReconstruction.so");
   if (inDir == "")
      inDir = "/mnt/f/a1954_C14_reco_hdb_slim/";

   TString dir = gSystem->DirName(__FILE__);
   TString plotDir = dir + "/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);
   TString outFile = plotDir + "ic_gate_C14" + outTag + ".root";
   TString stamp = TString::Format("%s|tb%d-%d", runsCSV.Data(), icTbLo, icTbHi);

   TH1D *hAmp = new TH1D("hAmp", ";IC amplitude [ADC];events", 500, 0, 4096);
   hAmp->SetDirectory(nullptr);
   long nEvt = 0, nIC = 0;
   int nRead = 0;
   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(ff)) {
         printf("\033[1;31mmissing %s -- not writing a partial spectrum under the full stamp\033[0m\n", ff.Data());
         return;
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
      ++nRead;
      printf("processed %s (%lld events)\n", run.Data(), N);
   }

   TFile *fo = TFile::Open(outFile, "RECREATE");
   hAmp->Write("hAmp");
   TNamed st("stamp", stamp.Data());
   st.Write("stamp");
   fo->Close();
   printf("\n%d runs, %ld events, %ld IC signals -> %s\n", nRead, nEvt, nIC, outFile.Data());
}
