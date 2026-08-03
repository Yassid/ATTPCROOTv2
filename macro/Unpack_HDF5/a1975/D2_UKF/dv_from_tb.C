/// @file dv_from_tb.C
/// @brief Measure the drift velocity from the raw time-bucket span of the hits, with no fit,
///        no kinematics and no beam energy anywhere in it.
///
/// Every hit carries the time bucket at which its charge arrived at the pad plane, and TB is
/// the one quantity in the chain that does NOT already assume a drift velocity -- AtHit's z is
/// computed FROM dv, so using z here would be circular. Charge released at the pad plane
/// arrives at once, charge released at the far end takes the full drift time, so the span of
/// the TB distribution is the time to cross the active length:
///
///     dv = driftLength / (dTB * TBtime)
///
/// TBtime is 160 ns exactly: AtDigiPar::GetTBTime() maps the parameter file's "SamplingRate 6"
/// to 160 ns, i.e. 6.25 MHz, not 6 MHz. With 512 buckets the whole window is only 81.92 us, so
/// crossing 100 cm inside it requires dv >= 1.2207 cm/us -- if the measured span turns out to
/// be the full 512 the acquisition was clipping the far end of the detector and the span is a
/// lower limit on dv rather than a measurement of it.
///
/// The edges are taken where the distribution crosses a fraction of its own plateau rather than
/// at the first and last non-empty bin, because single noise hits anywhere in the window would
/// otherwise set the answer. Both edge definitions are printed so the difference is visible.
///
///   root -b -q 'dv_from_tb.C("/mnt/f/a1975/reco_d2/run_0016_multifit_reco.root")'

void dv_from_tb(TString file = "/mnt/f/a1975/reco_d2/run_0016_multifit_reco.root", Long64_t nEvents = 4000,
                double driftLen_cm = 100.0, double tbTime_ns = 160.0, double frac = 0.5,
                TString plotOut = "plots/dv_from_tb.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", file.Data()); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("no tree cbmsim\n"); return; }
   TClonesArray *evArr = nullptr;
   t->SetBranchAddress("AtEventH", &evArr);

   auto *h = new TH1D("htb", "hit time buckets;time bucket;hits", 512, 0, 512);
   auto *hq = new TH1D("htbq", "charge-weighted;time bucket;charge", 512, 0, 512);
   Long64_t N = t->GetEntries();
   if (nEvents > 0 && nEvents < N) N = nEvents;
   long nHit = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!evArr || evArr->GetEntries() == 0) continue;
      auto *ev = dynamic_cast<AtEvent *>(evArr->At(0));
      if (!ev) continue;
      for (const auto &hit : ev->GetHits()) {
         if (!hit) continue;
         double tb = hit->GetTimeStamp();
         h->Fill(tb);
         hq->Fill(tb, hit->GetCharge());
         ++nHit;
      }
   }
   printf("\n=== dv_from_tb: %s ===\n", gSystem->BaseName(file));
   printf("%lld events, %ld hits, TB duration %.0f ns (%.2f MHz), drift length %.0f cm\n", N, nHit, tbTime_ns,
          1e3 / tbTime_ns, driftLen_cm);
   if (nHit < 1000) { printf("too few hits\n"); return; }

   // --- extremes: first and last non-empty bucket (sets an upper bound on the span) --------
   int fLo = 0, fHi = 0;
   for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinContent(b) > 0) { fLo = b; break; }
   for (int b = h->GetNbinsX(); b >= 1; --b) if (h->GetBinContent(b) > 0) { fHi = b; break; }

   // --- plateau edges: where the distribution crosses `frac` of its own median plateau -----
   std::vector<double> v;
   for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinContent(b) > 0) v.push_back(h->GetBinContent(b));
   std::sort(v.begin(), v.end());
   double plateau = v.empty() ? 0 : v[v.size() / 2];
   double thr = frac * plateau;
   int pLo = 0, pHi = 0;
   for (int b = 1; b <= h->GetNbinsX(); ++b) if (h->GetBinContent(b) > thr) { pLo = b; break; }
   for (int b = h->GetNbinsX(); b >= 1; --b) if (h->GetBinContent(b) > thr) { pHi = b; break; }

   auto report = [&](const char *tag, int lo, int hi) {
      double dtb = hi - lo;
      double us = dtb * tbTime_ns * 1e-3;
      printf("%-22s TB %3d -> %3d   dTB %5.1f   drift %6.2f us   dv = %.4f cm/us\n", tag, lo, hi, dtb, us,
             us > 0 ? driftLen_cm / us : -1);
   };
   printf("\nplateau (median of non-empty buckets) = %.0f hits/TB, edge threshold %.0f\n", plateau, thr);
   report("first/last non-empty", fLo, fHi);
   report(Form("%.0f%% of plateau", 100 * frac), pLo, pHi);
   printf("\nfor reference, dv if the FULL 512 TB window spanned %.0f cm: %.4f cm/us\n", driftLen_cm,
          driftLen_cm / (512 * tbTime_ns * 1e-3));
   printf("production value in ATTPC.a1975_deuterium.par: 1.15 cm/us\n");
   if (fHi >= 511 || fLo <= 1)
      printf("\nWARNING: the distribution reaches the edge of the 512-TB window, so the far end of the\n"
             "detector was clipped by the acquisition -- this span is a LOWER LIMIT on dv, not a value.\n");

   auto *cv = new TCanvas("cvtb", "dv_from_tb", 1100, 800);
   cv->Divide(1, 2);
   cv->cd(1);
   gPad->SetLogy();
   h->Draw("hist");
   for (int b : {pLo, pHi}) { auto *l = new TLine(b, 0, b, h->GetMaximum()); l->SetLineColor(kRed); l->Draw(); }
   cv->cd(2);
   gPad->SetLogy();
   hq->Draw("hist");
   gSystem->Exec("mkdir -p " + TString(gSystem->DirName(plotOut)));
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
