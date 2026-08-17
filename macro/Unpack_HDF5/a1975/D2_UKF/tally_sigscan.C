/// @file tally_sigscan.C
/// @brief Collapse rate vs assumed per-cluster sigma, for a collapsing run and a clean one.
///
/// Called by scan_meassigma_matfx.sh. Counts ndf <= 0 (or chi2 <= 0) against total fitted
/// tracks in each (run, measSigma) output, and writes a two-curve plot.

#include <cstdio>
#include <sstream>
#include <vector>

void tally_sigscan(TString base, TString sigmas, TString runList)
{
   gStyle->SetOptStat(0);
   gSystem->Load("libAtReconstruction.so");

   std::vector<double> sig;
   {
      std::stringstream ss(sigmas.Data());
      double v;
      while (ss >> v) sig.push_back(v);
   }
   std::vector<TString> runs;
   {
      std::stringstream ss(runList.Data());
      std::string r;
      while (ss >> r) runs.push_back(r);
   }

   auto *c = new TCanvas("c", "", 1000, 700);
   auto *mg = new TMultiGraph();
   auto *leg = new TLegend(0.58, 0.74, 0.88, 0.88);
   leg->SetBorderSize(0); leg->SetFillStyle(0);
   const int cols[] = {kRed + 1, kAzure + 2};

   printf("\n%-10s %8s %9s %9s %9s\n", "run", "sigma", "fitted", "collapsed", "pct");
   for (size_t ir = 0; ir < runs.size(); ++ir) {
      auto *g = new TGraph();
      for (double s : sig) {
         // the shell made the directory with the sigma string verbatim ("4.0"), and %g
         // would render that as "4" -- try both spellings rather than guess.
         TString f = Form("%s/run_%s_s%.1f/run_%s_multifit_genfitter_t.root",
                          base.Data(), runs[ir].Data(), s, runs[ir].Data());
         if (gSystem->AccessPathName(f))
            f = Form("%s/run_%s_s%g/run_%s_multifit_genfitter_t.root",
                     base.Data(), runs[ir].Data(), s, runs[ir].Data());
         auto *fl = TFile::Open(f);
         if (!fl || fl->IsZombie()) { printf("%-10s %8.1f   MISSING\n", runs[ir].Data(), s); continue; }
         auto *t = (TTree *)fl->Get("cbmsim");
         TClonesArray *te = nullptr;
         t->SetBranchAddress("AtTrackingEvent", &te);
         long n = 0, bad = 0;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!te || te->GetEntries() == 0) continue;
            auto *ev = (AtTrackingEvent *)te->At(0);
            if (!ev) continue;
            for (auto &ft : ev->GetFittedTracks()) {
               if (!ft) continue;
               const auto &m = ft->GetTrackMetadata();
               if (!m) continue;
               ++n;
               if (!(m->GetNdf() > 0 && m->GetChi2() > 0)) ++bad;
            }
         }
         const double pct = n ? 100.0 * bad / n : 0.0;
         printf("%-10s %8.1f %9ld %9ld %8.1f%%\n", runs[ir].Data(), s, n, bad, pct);
         if (n) g->SetPoint(g->GetN(), s, pct);
         fl->Close();
      }
      g->SetMarkerStyle(20); g->SetMarkerSize(1.5);
      g->SetMarkerColor(cols[ir % 2]); g->SetLineColor(cols[ir % 2]); g->SetLineWidth(3);
      mg->Add(g, "LP");
      leg->AddEntry(g, Form("run %s", runs[ir].Data()), "lp");
   }

   mg->SetTitle("matFX fit collapse vs assumed cluster #sigma;measSigma  (mm);fits collapsed  (%)");
   mg->Draw("A");
   mg->GetYaxis()->SetRangeUser(-3, 103);
   leg->Draw();
   c->SaveAs("plots/sigscan_collapse.png");
   printf("\nwrote plots/sigscan_collapse.png\n");
   printf("A FLAT response in both runs kills the measurement-error hypothesis.\n\n");
}
