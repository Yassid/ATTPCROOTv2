/// @file icsum_C15p.C
/// @brief Reduce a big <run>_FRIB.root to a tiny <run>_ic.root IC summary.
///
/// The FRIB output persists 8 generic channels x 2048 samples per event -- ~36 kB/event,
/// i.e. ~54 GB over the a2091 runs. The analysis only ever uses TWO numbers per event:
///   icmax  = max of generic trace[0] over time buckets [icTbLo, icTbHi]   (the ion chamber)
///   npulse = pulse count over [pkTbLo, pkTbHi] above peakThr              (pile-up rejection)
/// exactly as gate_events_C15.C computes them. Storing those instead is ~25 MB for everything,
/// so the IC spectrum can be re-binned and the beam gate re-chosen instantly, and the 54 GB of
/// traces never has to exist. Defaults mirror gate_events_C15.C so the numbers are identical.
///
///   root -b -q 'pipeline/icsum_C15p.C("run_0017","/path/to/frib/","/home/yassid/a2091_C15_ic/")'

static int CountPulsesIC(const std::vector<Double_t> &adc, double thr, int tbLo, int tbHi)
{
   int n = 0, b = tbLo, N = (int)adc.size();
   if (tbHi > N) tbHi = N;
   while (b < tbHi) {
      if (adc[b] > thr) { ++n; while (b < tbHi && adc[b] > thr * 0.5) ++b; }
      else ++b;
   }
   return n;
}

void icsum_C15p(TString run, TString inDir = "/home/yassid/a2091_C15_ic/tmp/",
               TString outDir = "/home/yassid/a2091_C15_ic/", Int_t icTbLo = 1050, Int_t icTbHi = 1250,
               double peakThr = 200, Int_t pkTbLo = 800, Int_t pkTbHi = 1500)
{
   gSystem->Load("libAtReconstruction.so");
   TString fin = inDir + run + "_FRIB.root";
   if (gSystem->AccessPathName(fin)) { printf("MISSING %s\n", fin.Data()); return; }

   TFile *f = TFile::Open(fin);
   TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
   if (!t || t->GetEntries() == 0) { printf("EMPTY %s\n", fin.Data()); if (f) f->Close(); return; }
   TClonesArray *ra = nullptr;
   t->SetBranchAddress("AtRawEvent", &ra);

   gSystem->mkdir(outDir.Data(), kTRUE);
   TFile *fo = new TFile((outDir + run + "_ic.root").Data(), "RECREATE");
   Int_t entry, npulse, ictb, evtid;
   Float_t icmax;
   TTree *ot = new TTree("ic", "IC summary");
   // ★ TWO INDICES, AND THEY ARE NOT THE SAME THING.
   //   entry  = position in the FRIB tree
   //   evtid  = the TRUE event number, carried on the AtRawEvent
   // /get and /frib live in the SAME merged file and both name their data by event number, so
   // event N is the same trigger in both by construction. But the two DAQs do not record the same
   // set: run_0151 has 39993 GET against 38595 here, run_0140 31800 against 32295. From the first
   // gap onward the tree POSITION is off while the event NUMBER still matches, so joining on
   // `entry` silently mismatches everything after it. Joining on `evtid` does not.
   ot->Branch("entry", &entry, "entry/I");   // position in the FRIB tree
   ot->Branch("evtid", &evtid, "evtid/I");   // the true event number -- JOIN ON THIS
   ot->Branch("icmax", &icmax, "icmax/F");
   ot->Branch("npulse", &npulse, "npulse/I");
   ot->Branch("ictb", &ictb, "ictb/I");   // time bucket of the maximum

   Long64_t N = t->GetEntries(), nfill = 0;
   TH1F *h = new TH1F("hic", (run + ": IC max;IC max [ADC];events").Data(), 400, 0, 4000);
   for (Long64_t i = 0; i < N; i++) {
      t->GetEntry(i);
      entry = (Int_t)i; icmax = -1; npulse = 0; ictb = -1; evtid = -1;
      if (ra && ra->GetEntries() > 0) {
         auto *r = (AtRawEvent *)ra->At(0);
         if (r) evtid = (Int_t)r->GetEventID();
         if (r && !r->GetGenTraces().empty()) {
            auto &adc = r->GetGenTraces()[0]->GetADC();
            // ★ TAKE THE MAXIMUM OVER THE SAME RANGE THE PULSES ARE COUNTED OVER.
            // This took the max over the NARROW [icTbLo,icTbHi]=[1050,1250] while counting pulses
            // over [pkTbLo,pkTbHi]=[800,1500]. Pulses arrive across the whole counting range, so
            // for an event whose pulse falls outside the narrow window the maximum was the
            // BASELINE -- while npulse still said 1, because the counter did see it. Those events
            // then fail any beam window and vanish. Measured at 12.9 % on these a2091 proton runs,
            // showing up as a spurious peak near 35-45 ADC that is pure artefact.
            double mx = -1e9;
            int mxTb = -1;
            for (int b = pkTbLo; b < pkTbHi && b < (int)adc.size(); b++)
               if (adc[b] > mx) { mx = adc[b]; mxTb = b; }
            icmax = (Float_t)mx;
            ictb = mxTb;
            npulse = CountPulsesIC(adc, peakThr, pkTbLo, pkTbHi);
            h->Fill(icmax);
            nfill++;
         }
      }
      ot->Fill();
   }
   fo->cd(); ot->Write(); h->Write();
   printf("%-10s entries=%-8lld withIC=%-8lld  mean=%.1f  peak=%.0f ADC  1pulse=%.1f%%\n", run.Data(), N, nfill,
          h->GetMean(), h->GetBinCenter(h->GetMaximumBin()),
          100.0 * ot->GetEntries("npulse==1") / std::max((Long64_t)1, N));
   fo->Close(); f->Close();
}
