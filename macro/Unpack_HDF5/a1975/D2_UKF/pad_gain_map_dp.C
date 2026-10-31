/// @file pad_gain_map_dp.C
/// @brief Data-driven PAD-GAIN map for the a1975 D2 AT-TPC, from the tracks themselves.
///
/// Each good track (AtHelix3D valid, coherent, >= 30 hits) is its own reference: for every hit on the
/// small pads (rho <= 152 mm) with an UNSATURATED amplitude (< aMax: the GET ceiling clips pulses at
/// ~4200 ADC, measured 2026-10-31), r = amplitude / median amplitude of that track's qualifying hits.
/// The track's dE/dx and most of the crossing geometry cancel in r; averaged over many tracks, a pad's
/// mean ln r is its relative gain. Robust: |ln r| < 1.5 kept, pads with fewer than minHits entries get 1.
/// Output: text map `pad gain nHits` (gain normalised to median 1 over calibrated pads) + PNG of the map.
///   root -l -b -q 'pad_gain_map_dp.C("0016,0020,0026","pid/padgain_d2.txt","plots/padgain_d2.png")'
void pad_gain_map_dp(TString runsCSV, TString out, TString png, TString recoDir = "/mnt/f/a1975/reco_d2_dv1104/",
                     double aMax = 3500, int minHits = 50)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gErrorIgnoreLevel = kError;
   gStyle->SetOptStat(0);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   struct Acc { double s = 0, s2 = 0, x = 0, y = 0; long n = 0; };
   std::map<int, Acc> pad;
   long nTr = 0;
   for (auto *o : *runsCSV.Tokenize(",")) {
      TString r = ((TObjString *)o)->GetString();
      // a token starting with '/' is a full file path (e.g. a simulation reco), otherwise a run number
      TFile f(r.BeginsWith("/") ? r : recoDir + "run_" + r + "_multifit_reco.root");
      auto *t = (TTree *)f.Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      t->SetBranchAddress("AtPatternEvent", &pa);
      for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
         t->GetEntry(ev);
         if (!pa || !pa->GetEntries()) continue;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            auto &H = tr.GetHitArray();
            if (H.size() < 30) continue;
            auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
            auto hr = h3.Estimate(tr, sr.valid ? sr.direction : -1);
            if (!hr.valid || hr.spacing >= 6) continue;
            std::vector<double> q;
            for (auto &h : H) {
               auto p = h->GetPosition();
               if (std::hypot(p.X(), p.Y()) <= 152 && h->GetCharge() > 0 && h->GetCharge() < aMax) q.push_back(h->GetCharge());
            }
            if (q.size() < 15) continue;
            std::nth_element(q.begin(), q.begin() + q.size() / 2, q.end());
            const double med = q[q.size() / 2];
            ++nTr;
            for (auto &h : H) {
               auto p = h->GetPosition();
               double a = h->GetCharge();
               if (std::hypot(p.X(), p.Y()) > 152 || a <= 0 || a >= aMax) continue;
               double lr = std::log(a / med);
               if (std::abs(lr) > 1.5) continue;
               auto &A = pad[h->GetPadNum()];
               A.s += lr; A.s2 += lr * lr; A.x += p.X(); A.y += p.Y(); ++A.n;
            }
         }
      }
      printf("run %s done, %ld tracks so far, %zu pads seen\n", r.Data(), nTr, pad.size());
   }
   // gains, normalised to median 1 over calibrated pads
   std::vector<double> lg;
   for (auto &[p, A] : pad) if (A.n >= minHits) lg.push_back(A.s / A.n);
   std::sort(lg.begin(), lg.end());
   const double lmed = lg[lg.size() / 2];
   std::ofstream fo(out.Data());
   fo << "# a1975 D2 pad-gain map (pad_gain_map_dp.C): pad gain nHits x y gainLocal; gain = exp(mean ln(amp/track median)"
      << " - median); gainLocal = gain / median gain of pads at the same radius (5 mm rings): pad-to-pad only,\n"
      << "# without the smooth radial trend (part of which is PHYSICS: tracks start near the axis at low dE/dx)\n"
      << "# runs " << runsCSV.Data() << ", small pads only, amplitude < " << aMax << ", pads with < " << minHits
      << " hits set to 1\n";
   auto *g = new TGraph2D;
   std::vector<double> gs;
   int nCal = 0;
   // ring medians (5 mm in rho) of the calibrated gains, for the local (pad-to-pad) version
   std::map<int, std::vector<double>> ring;
   for (auto &[p, A] : pad)
      if (A.n >= minHits) ring[int(std::hypot(A.x / A.n, A.y / A.n) / 5)].push_back(std::exp(A.s / A.n - lmed));
   std::map<int, double> ringMed;
   for (auto &[k, v] : ring) { std::sort(v.begin(), v.end()); ringMed[k] = v[v.size() / 2]; }
   std::vector<double> gl;
   for (auto &[p, A] : pad) {
      double gain = A.n >= minHits ? std::exp(A.s / A.n - lmed) : 1.0;
      double x = A.n ? A.x / A.n : 0, y = A.n ? A.y / A.n : 0;
      double gLoc = A.n >= minHits ? gain / ringMed[int(std::hypot(x, y) / 5)] : 1.0;
      if (A.n >= minHits) gl.push_back(gLoc);
      fo << p << " " << Form("%.4f", gain) << " " << A.n << " " << Form("%.1f %.1f %.4f", x, y, gLoc) << "\n";
      if (A.n >= minHits) {
         ++nCal;
         gs.push_back(gain);
         g->SetPoint(g->GetN(), A.x / A.n, A.y / A.n, gain);
      }
   }
   std::sort(gs.begin(), gs.end());
   int n = gs.size();
   printf("%ld tracks; %zu pads seen, %d calibrated (>= %d hits)\n", nTr, pad.size(), nCal, minHits);
   std::sort(gl.begin(), gl.end());
   printf("LOCAL pad-to-pad spread: p2 %.3f p16 %.3f median %.3f p84 %.3f p98 %.3f (half 16-84 width %.3f)\n",
          gl[gl.size() * 2 / 100], gl[gl.size() * 16 / 100], gl[gl.size() / 2], gl[gl.size() * 84 / 100],
          gl[gl.size() * 98 / 100], 0.5 * (gl[gl.size() * 84 / 100] - gl[gl.size() * 16 / 100]));
   printf("pad gain spread: p2 %.3f p16 %.3f median %.3f p84 %.3f p98 %.3f  (half 16-84 width %.3f)\n", gs[n * 2 / 100],
          gs[n * 16 / 100], gs[n / 2], gs[n * 84 / 100], gs[n * 98 / 100], 0.5 * (gs[n * 84 / 100] - gs[n * 16 / 100]));
   auto *c = new TCanvas("c", "c", 900, 800);
   gPad->SetRightMargin(0.15);
   auto *h = new TH2D("hg", "relative pad gain (small pads);x [mm];y [mm]", 160, -160, 160, 160, -160, 160);
   for (int i = 0; i < g->GetN(); ++i) h->Fill(g->GetX()[i], g->GetY()[i], g->GetZ()[i]);
   auto *hn = new TH2D("hn", "", 160, -160, 160, 160, -160, 160);
   for (int i = 0; i < g->GetN(); ++i) hn->Fill(g->GetX()[i], g->GetY()[i]);
   h->Divide(hn);
   h->SetMinimum(0.7); h->SetMaximum(1.3);
   h->Draw("colz");
   c->SaveAs(png);
   printf("wrote %s\n", out.Data());
}
