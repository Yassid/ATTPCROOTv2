/// @file icspec_C15.C
/// @brief Combined ion-chamber spectrum for a2091 15C from the per-run _ic.root summaries.
///
/// Answers the two questions needed before any beam gate can be chosen:
///   1. where are the IC peaks, and which window isolates ONE of them?
///      (the inherited 950-1350 placeholder cuts straight through the ~1380 peak)
///   2. is the IC gain STABLE across the 40 runs? A single window is only meaningful if the
///      peak does not drift run to run -- otherwise the gate has to be per-run.
///
///   root -b -q 'pid/icspec_C15.C()'
void icspec_C15(TString icDir = "/home/yassid/a2091_C15_ic/", double lo = 0, double hi = 2600)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString plotDir = TString(gSystem->DirName(__FILE__)) + "/plots/";

   void *dirp = gSystem->OpenDirectory(icDir);
   if (!dirp) { printf("ERROR: cannot open %s\n", icDir.Data()); return; }
   std::vector<TString> runs;
   const char *e;
   while ((e = gSystem->GetDirEntry(dirp))) {
      TString f(e);
      if (f.EndsWith("_ic.root")) runs.push_back(f(0, f.Length() - 8));
   }
   gSystem->FreeDirectory(dirp);
   std::sort(runs.begin(), runs.end());
   printf("found %zu run summaries\n\n", runs.size());

   const int NB = 520; // 5 ADC per bin
   TH1F *hAll = new TH1F("hAll", "a2091 15C: combined IC;IC max [ADC];events", NB, lo, hi);
   TH1F *hOne = new TH1F("hOne", "combined IC, single pulse;IC max [ADC];events", NB, lo, hi);
   TH2F *hRun = new TH2F("hRun", "IC vs run;run index;IC max [ADC]", runs.size(), 0, runs.size(), NB, lo, hi);

   printf("%-11s %10s %10s %9s %9s %9s\n", "run", "events", "withIC", "1pulse%", "peak", "mean");
   std::vector<double> peaks;
   long totEv = 0, totOne = 0;
   for (size_t i = 0; i < runs.size(); i++) {
      TFile *f = TFile::Open(icDir + runs[i] + "_ic.root");
      TTree *t = f ? (TTree *)f->Get("ic") : nullptr;
      if (!t) { if (f) f->Close(); continue; }
      Int_t entry, npulse; Float_t icmax;
      t->SetBranchAddress("entry", &entry); t->SetBranchAddress("icmax", &icmax);
      t->SetBranchAddress("npulse", &npulse);
      TH1F hr("hr", "", NB, lo, hi);
      long n = 0, n1 = 0;
      for (Long64_t k = 0; k < t->GetEntries(); k++) {
         t->GetEntry(k);
         if (icmax < 0) continue;
         n++;
         hAll->Fill(icmax); hr.Fill(icmax);
         hRun->Fill(i + 0.5, icmax);
         if (npulse == 1) { hOne->Fill(icmax); n1++; }
      }
      // peak position above 500 ADC, ignoring the low-amplitude junk spike
      int b0 = hr.FindBin(500), bm = b0; double best = -1;
      for (int b = b0; b <= NB; b++) if (hr.GetBinContent(b) > best) { best = hr.GetBinContent(b); bm = b; }
      double pk = hr.GetBinCenter(bm);
      peaks.push_back(pk);
      printf("%-11s %10ld %10ld %8.1f%% %9.0f %9.1f\n", runs[i].Data(), (long)t->GetEntries(), n,
             n ? 100.0 * n1 / n : 0, pk, hr.GetMean());
      totEv += n; totOne += n1;
      f->Close();
   }

   double pmin = *std::min_element(peaks.begin(), peaks.end());
   double pmax = *std::max_element(peaks.begin(), peaks.end());
   double psum = 0; for (double p : peaks) psum += p;
   double pmean = psum / peaks.size();
   double pvar = 0; for (double p : peaks) pvar += (p - pmean) * (p - pmean);
   printf("\n==== TOTAL: %ld events with IC, %ld single-pulse (%.1f%%) ====\n", totEv, totOne,
          100.0 * totOne / std::max(1L, totEv));
   printf("per-run peak: mean %.0f  spread %.0f (rms)  range %.0f - %.0f  -> %s\n", pmean,
          std::sqrt(pvar / peaks.size()), pmin, pmax,
          (pmax - pmin) < 100 ? "STABLE, one global window is fine" : "DRIFTING, consider a per-run gate");

   // find the peaks in the combined single-pulse spectrum
   printf("\n---- peak search, combined single-pulse spectrum ----\n");
   TSpectrum sp(12);
   int npk = sp.Search(hOne, 2, "nodraw", 0.05);
   double *xp = sp.GetPositionX(); double *yp = sp.GetPositionY();
   std::vector<std::pair<double, double>> pk;
   for (int i = 0; i < npk; i++) pk.push_back({xp[i], yp[i]});
   std::sort(pk.begin(), pk.end());
   for (auto &p : pk)
      if (p.first > 200) printf("  peak at %7.0f ADC   height %8.0f\n", p.first, p.second);

   printf("\n---- counts in candidate windows (single pulse) ----\n");
   auto frac = [&](double a, double b) {
      double c = hOne->Integral(hOne->FindBin(a), hOne->FindBin(b));
      printf("  %5.0f - %-5.0f : %9.0f  (%5.1f%% of single-pulse)\n", a, b, c, 100.0 * c / hOne->Integral());
   };
   frac(950, 1350);   // the inherited placeholder -- cuts through the second peak
   frac(1000, 1250);
   frac(1250, 1500);
   frac(1900, 2150);

   TCanvas *c = new TCanvas("c", "ic", 1700, 950);
   c->Divide(2, 2);
   c->cd(1); gPad->SetLogy(); hAll->Draw();
   c->cd(2); gPad->SetLogy(); hOne->Draw();
   c->cd(3); gPad->SetLogz(); hRun->Draw("colz");
   c->cd(4); hOne->Draw(); hOne->GetXaxis()->SetRangeUser(800, 1700);
   c->SaveAs(plotDir + "ic_spectrum_C15.png");
   printf("\nwrote %sic_spectrum_C15.png\n", plotDir.Data());
}
