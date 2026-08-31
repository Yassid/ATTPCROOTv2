/// @file drift_check_C15d.C
/// @brief Measure the drift working point FROM THE DATA, with no kinematics involved.
///
///   root -b -q 'drift_check_C15d.C()'
///   root -b -q 'drift_check_C15d.C("/home/yassid/C15d_reco/", 6)'   // more runs
///
/// ★ WHY THIS EXISTS. The par file's drift velocity was derived from Spyral's
/// `window_time_bucket` / `micromegas_time_bucket` rather than measured, and it was wrong by a
/// factor 0.574 -- which nothing downstream revealed, because the error propagates COHERENTLY.
/// A wrong z scale stretches theta, hence the fitted KE, hence Brho (AtSpyralPID divides by
/// sin(polar)), hence Ex. The whole picture stayed self-consistent at a bogus beam energy while
/// the elastic (d,d) peak sat 2.9 MeV from zero, and every kinematic fit I threw at it was happy.
///
/// What caught it was a DETECTOR consistency check: reconstructed hit z spanned -20 to 1197 mm
/// inside a 1000 mm chamber. Hits cannot be outside the chamber. Run this BEFORE trusting any
/// kinematics from a new par.
///
/// THE TWO ANCHORS, and why one is easy and one is not:
///   * the PAD PLANE end is a sharp lower edge in the hit-time distribution -- easy, stable.
///   * the WINDOW end is NOT the edge of the hit distribution, because the DAQ keeps recording
///     past it: with 512 samples the hits run to the end of the record whatever the drift is.
///     Take it from the LATEST-ARRIVING charge per track, which piles up at the window and stops.
///
/// The measurement constrains mm PER TIME BUCKET, i.e. the PRODUCT of drift velocity and sampling
/// period. It cannot separate them -- that needs the DAQ run log.

void drift_check_C15d(TString recoDir = "/home/yassid/C15d_reco/", Int_t nRuns = 6,
                      Double_t mmPerTB = 3.63636, Double_t zPadPlane = 1000.0, Double_t tbEntrance = 300.0,
                      Double_t driftLengthMM = 1000.0, Long64_t maxEvents = 8000)
{
   gSystem->Load("libAtReconstruction.so");
   std::cout << "\033[1;33m=== C15d drift working point, measured ===\033[0m\n"
             << "  assuming the par currently in use: " << mmPerTB << " mm/TB, ZPadPlane " << zPadPlane
             << ", TBEntrance " << tbEntrance << "\n\n";

   TSystemDirectory dir(recoDir, recoDir);
   TList *files = dir.GetListOfFiles();
   if (!files) { std::cout << "cannot list " << recoDir << "\n"; return; }
   std::vector<TString> use;
   TIter next(files);
   while (auto *o = dynamic_cast<TSystemFile *>(next())) {
      TString n = o->GetName();
      if (!o->IsDirectory() && n.EndsWith("_reco.root")) use.push_back(recoDir + n);
   }
   std::sort(use.begin(), use.end());
   if (use.empty()) { std::cout << "no recos in " << recoDir << "\n"; return; }

   printf("  %-14s %10s %9s %9s %9s\n", "run", "hits", "TB lo", "TB win", "span");
   double accLo = 0, accHi = 0;
   int nUsed = 0;
   const int step = std::max(1, (int)(use.size() / std::max(1, nRuns)));
   for (size_t k = 0; k < use.size() && nUsed < nRuns; k += step) {
      TFile f(use[k]);
      auto *t = (TTree *)f.Get("cbmsim");
      if (!t) continue;
      TClonesArray *pa = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pa);
      TH1D hall("hall", "", 540, -20, 520), hmax("hmax", "", 540, -20, 520);
      long nh = 0;
      const Long64_t n = std::min(maxEvents, t->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         t->GetEntry(i);
         if (!pa || pa->GetEntriesFast() == 0) continue;
         auto *p = (AtPatternEvent *)pa->At(0);
         if (!p) continue;
         for (auto &trk : p->GetTrackCand()) {
            double mx = -1e9;
            for (auto &h : trk.GetHitArray()) {
               const double tb = tbEntrance - (zPadPlane - h->GetPosition().Z()) / mmPerTB;
               hall.Fill(tb);
               if (tb > mx) mx = tb;
               ++nh;
            }
            if (mx > -1e8) hmax.Fill(mx);
         }
      }
      if (nh < 20000) continue;
      // pad plane: sharp lower edge of ALL hits
      double mxAll = hall.GetMaximum();
      int lo = -1;
      for (int b = 1; b <= hall.GetNbinsX(); ++b)
         if (hall.GetBinContent(b) > 0.02 * mxAll) { lo = b; break; }
      // window: the top end of the LATEST-arriving charge, which piles up there and stops
      int hi = -1;
      for (int b = hmax.GetNbinsX(); b >= 1; --b)
         if (hmax.GetBinContent(b) > 0) { hi = b; break; }
      const double tlo = hall.GetBinCenter(lo), thi = hmax.GetBinCenter(hi);
      printf("  %-14s %10ld %9.1f %9.1f %9.1f\n", gSystem->BaseName(use[k].Data()), nh, tlo, thi, thi - tlo);
      accLo += tlo; accHi += thi; ++nUsed;
   }
   if (!nUsed) { std::cout << "  no run had enough hits\n"; return; }

   const double tlo = accLo / nUsed, thi = accHi / nUsed, span = thi - tlo;
   const double mmtb = driftLengthMM / span;
   printf("\n  mean over %d runs: pad plane TB %.1f, window TB %.1f, span %.1f TB\n", nUsed, tlo, thi, span);
   printf("  \033[1;32mmm per time bucket = %.3f\033[0m  (par implies %.3f, factor %.3f)\n",
          mmtb, mmPerTB, mmtb / mmPerTB);
   printf("  TBEntrance should be ~%.0f (par says %.0f)\n\n", thi, tbEntrance);
   printf("  the measurement fixes the PRODUCT dv * TBTime, not the split:\n");
   for (double ns : {320.0, 160.0, 80.0})
      printf("     TBTime %5.0f ns  ->  dv = %.4f cm/us\n", ns, mmtb / (ns * 1e-9) / 1e7);
   printf("\n  Take the sampling frequency from the DAQ run log; the HDF5 /meta group carries only\n"
          "  GRAW filenames, no frequency.\n");
}
