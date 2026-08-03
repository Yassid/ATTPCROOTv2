/// @file dv_from_raw.C
/// @brief Drift velocity from the RAW pad traces, with no PSA anywhere in the chain.
///
/// The first version of this measurement used AtEvent hits and got 1.305 cm/us from three runs,
/// in apparent agreement with the parameter file's own TBEntrance. It was wrong. AtPSAMultiFit
/// gates every peak with
///        if (amp < ampCut || peakTB < 20 || peakTB > 500)
/// and searches only i = 21 .. 499, so the hit distribution is clipped to 21-500 by construction.
/// The measured "span" was that window, which is why every run gave the same answer and why the
/// edges were vertical walls rather than physical roll-offs.
///
/// So this reads AtRawEvent straight from the HDF5 GET stream and histograms the time bucket of
/// the samples themselves. Nothing here can impose a TB window: the only cut is on amplitude
/// above baseline, and the macro prints the distribution near both ends so a software edge would
/// be visible as a cliff instead of a roll-off.
///
///     dv = driftLength / (dTB * TBtime),  TBtime = 160 ns (SamplingRate 6 -> 6.25 MHz)
///
/// The physical picture: charge released at the pad plane arrives immediately, charge from the
/// far end takes the full drift time, so the span of the distribution is the crossing time. The
/// beam has to actually traverse the full length for the upper edge to mean anything, which is
/// why the per-event maximum TB is printed too -- if the beam stops short, that edge is the
/// range of the beam and not the length of the detector.
///
///   root -b -q 'dv_from_raw.C("run_0016", 300)'

void dv_from_raw(TString run = "run_0016", Long64_t nEvents = 300, double thr = 25.0, double driftLen_cm = 100.0,
                 double tbTime_ns = 160.0, int tbIgnoreEnd = 4, TString plotOut = "plots/dv_from_raw.png")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gStyle->SetOptStat(0);

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = "/mnt/f/a1975/h5/" + run + ".h5";
   if (gSystem->AccessPathName(inputFile)) { printf("no %s\n", inputFile.Data()); return; }
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());

   FairRunAna *fRun = new FairRunAna();
   fRun->SetOutputFile("/tmp/dv_from_raw_dummy.root");
   fRun->SetGeomFile(dir + "/geometry/ATTPC_H1bar_geomanager.root");
   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/ATTPC.a1975_deuterium.par").Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");

   auto map = std::make_shared<AtTpcMap>();
   map->ParseXMLMap((dir + "/scripts/ANL2023.xml").Data());
   map->GeneratePadPlane();

   auto unpacker = std::make_unique<AtHDFUnpacker>(map);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(1);
   unpacker->SetBaseLineSubtraction(true);
   unpacker->Init();

   auto *hAll = new TH1D("hAll", "raw samples above threshold;time bucket;samples", 512, 0, 512);
   auto *hMax = new TH1D("hMax", "per-pad peak TB;time bucket;pads", 512, 0, 512);
   auto *hEvMax = new TH1D("hEvMax", "per-event MAX pad-peak TB;time bucket;events", 512, 0, 512);
   auto *hEvMin = new TH1D("hEvMin", "per-event MIN pad-peak TB;time bucket;events", 512, 0, 512);

   Long64_t N = unpacker->GetNumEvents();
   if (nEvents > 0 && nEvents < N) N = nEvents;
   printf("\n=== dv_from_raw: %s, %lld events, threshold %.0f ADC over baseline ===\n", run.Data(), N, thr);

   TClonesArray arr("AtRawEvent", 1);
   long nPad = 0;
   for (Long64_t i = 0; i < N; ++i) {
      arr.Clear("C");
      auto *rawp = dynamic_cast<AtRawEvent *>(arr.ConstructedAt(0));
      unpacker->FillRawEvent(*rawp);
      int evLo = 9999, evHi = -1;
      for (const auto &pad : rawp->GetPads()) {
         if (!pad) continue;
         const auto &adc = pad->GetADC();
         int best = -1;
         double bestA = thr;
         // the last few samples of every record carry an end-of-trace artifact -- every event
         // has a pad "peaking" at 512 -- so they are excluded rather than allowed to define the edge
         for (std::size_t tb = 0; tb + tbIgnoreEnd < adc.size(); ++tb) {
            if (adc[tb] > thr) hAll->Fill(tb);
            if (adc[tb] > bestA) { bestA = adc[tb]; best = tb; }
         }
         if (best >= 0) {
            hMax->Fill(best);
            ++nPad;
            evLo = std::min(evLo, best);
            evHi = std::max(evHi, best);
         }
      }
      if (evHi >= 0) { hEvMax->Fill(evHi); hEvMin->Fill(evLo); }
      if (unpacker->IsLastEvent()) break;
   }
   printf("%ld pads with a peak above threshold\n", nPad);
   if (nPad < 500) { printf("too few\n"); return; }

   auto edges = [&](TH1D *hh, const char *tag) {
      std::vector<double> v;
      for (int b = 1; b <= hh->GetNbinsX(); ++b) if (hh->GetBinContent(b) > 0) v.push_back(hh->GetBinContent(b));
      std::sort(v.begin(), v.end());
      double plateau = v.empty() ? 0 : v[v.size() / 2];
      int lo = 0, hi = 0, flo = 0, fhi = 0;
      for (int b = 1; b <= hh->GetNbinsX(); ++b) if (hh->GetBinContent(b) > 0.5 * plateau) { lo = b; break; }
      for (int b = hh->GetNbinsX(); b >= 1; --b) if (hh->GetBinContent(b) > 0.5 * plateau) { hi = b; break; }
      for (int b = 1; b <= hh->GetNbinsX(); ++b) if (hh->GetBinContent(b) > 0) { flo = b; break; }
      for (int b = hh->GetNbinsX(); b >= 1; --b) if (hh->GetBinContent(b) > 0) { fhi = b; break; }
      double dtb = hi - lo, us = dtb * tbTime_ns * 1e-3;
      printf("%-22s nonzero %3d-%3d | half-plateau %3d-%3d  dTB %5.1f  %6.2f us  dv = %.4f cm/us\n", tag, flo, fhi, lo,
             hi, dtb, us, us > 0 ? driftLen_cm / us : -1);
   };
   printf("\n%-22s %s\n", "", "(a software window would show as a cliff at a round number)");
   edges(hAll, "all samples > thr");
   edges(hMax, "per-pad peak TB");
   printf("\nper-event extremes: min-TB peaks at %.0f, max-TB peaks at %.0f\n",
          hEvMin->GetBinCenter(hEvMin->GetMaximumBin()), hEvMax->GetBinCenter(hEvMax->GetMaximumBin()));
   printf("par: DriftVelocity 1.15, TBEntrance 480, DriftLength 1000 mm, NumTbs 512\n");
   printf("     TBEntrance 480 over 100 cm would imply dv = %.4f\n", driftLen_cm / (480 * tbTime_ns * 1e-3));

   auto *cv = new TCanvas("cvr", "dv_from_raw", 1150, 900);
   cv->Divide(1, 3);
   cv->cd(1); gPad->SetLogy(); hAll->Draw("hist");
   cv->cd(2); gPad->SetLogy(); hMax->Draw("hist");
   cv->cd(3); hEvMax->SetLineColor(kRed); hEvMax->Draw("hist"); hEvMin->Draw("hist same");
   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
