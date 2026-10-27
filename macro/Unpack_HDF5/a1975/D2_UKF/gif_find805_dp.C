/// @file gif_find805_dp.C
/// @brief Animated GIF of the ADOPTED FIND fits (gf_dp_find805) against the full-length baseline
/// (gf_dp_cateloss) on the same 16C(d,p) tracks, one track per frame. ROOT-only stopgap; the
/// polished version is meant to be done in VTK.
///
/// Per frame: XY (the spiral face-on) and XZ (beam axis horizontal). Grey = every cluster of the
/// pattern track; dark = clusters within kNearMM of a find805 fit point, i.e. the prefix FIND kept.
/// Blue = full-length fit, red = find805 fit, open/filled star = their vertices.
///
/// FIND WRITES NO findPct, so the kept fraction shown is (stored fit points / clusters), the same
/// proxy as find_rung_dp.C -- a LOWER BOUND, since some fits store fewer points than they used.
///
/// Matching is on (entry, trackID) as in cmp_gf_dirs.C. Duplicate tree entries (22 %) are skipped
/// by requiring a valid Spyral record, as in dump_tracks_dp.C. Fit z is in the fitter frame and is
/// flipped (z = 1000 - z_fit) onto the cluster frame; vertices too.
///
///   root -l -b -q 'gif_find805_dp.C("0016", 24, "out.gif")'
#include <cmath>
#include <map>
#include <set>
#include <vector>
void gif_find805_dp(TString run = "0016", int nFrames = 24, TString outGif = "find805_fits.gif",
                    int minCl = 60, int delayCs = 250,
                    /// entries forced into the animation first (e11488 is the 14 % case)
                    TString forceCSV = "11488",
                    TString dirFind = "/mnt/f/a1975/gf_dp_find805/",
                    TString dirFull = "/mnt/f/a1975/gf_dp_cateloss/",
                    /// non-empty: also write every frame as <pngPrefix>_NN.png, for checking by eye
                    TString pngPrefix = "")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   const double kZPadFit = 1000.0, kNearMM = 12.0;

   std::set<Long64_t> force;
   { TObjArray *t = forceCSV.Tokenize(",");
     for (int i = 0; i < t->GetEntries(); ++i) force.insert(((TObjString*)t->At(i))->GetString().Atoll()); }

   auto *fF = TFile::Open(dirFind + "run_" + run + "_multifit_genfitter_p.root");
   auto *fB = TFile::Open(dirFull + "run_" + run + "_multifit_genfitter_p.root");
   if (!fF || fF->IsZombie() || !fB || fB->IsZombie()) { printf("  cannot open both files\n"); return; }
   auto *tF = (TTree*)fF->Get("cbmsim"); auto *tB = (TTree*)fB->Get("cbmsim");
   TClonesArray *aF = nullptr, *pF = nullptr, *aB = nullptr;
   tF->SetBranchAddress("AtTrackingEvent", &aF); tF->SetBranchAddress("AtPIDEvent", &pF);
   tB->SetBranchAddress("AtTrackingEvent", &aB);
   const Long64_t N = std::min(tF->GetEntries(), tB->GetEntries());

   auto c2n = [](AtFittedTrack *f) { double n = f->GetTrackMetadata()->GetNdf();
                                     return n > 0 ? f->GetTrackMetadata()->GetChi2() / n : -1.0; };

   // Frame list: forced entries first, then truncated tracks spread over the run.
   struct Pick { Long64_t entry; int tid; };
   std::vector<Pick> picks;
   std::vector<Pick> cand;
   for (Long64_t i = 0; i < N; ++i) {
      tF->GetEntry(i);
      if (!pF || pF->GetEntries() == 0 || !aF || aF->GetEntries() == 0) continue;
      auto *pid = (AtPIDEvent*)pF->At(0); auto *ev = (AtTrackingEvent*)aF->At(0);
      if (!pid || !ev) continue;
      std::map<int, int> ncl;
      for (auto &tr : ev->GetTrackArray()) ncl[tr.GetTrackID()] = tr.GetHitClusterArray()->size();
      std::map<int, AtFittedTrack*> fm;
      for (auto &ft : ev->GetFittedTracks()) if (ft) fm[ft->GetTrackID()] = ft.get();
      for (auto &sr : pid->GetSpyral()) {
         if (!sr.valid) continue;
         auto it = fm.find(sr.trackID); if (it == fm.end() || !ncl.count(sr.trackID)) continue;
         int n = ncl[sr.trackID], np = it->second->GetSmoothedPositions().size();
         if (force.count(i)) { picks.push_back({i, sr.trackID}); continue; }
         if (n < minCl || np < 3) continue;
         if ((double)np / n < 0.75) cand.push_back({i, sr.trackID});   // FIND really truncated it
      }
   }
   const int want = std::max(0, nFrames - (int)picks.size());
   const int stride = (want > 0 && (int)cand.size() > want) ? cand.size() / want : 1;
   for (size_t k = 0; k < cand.size() && (int)picks.size() < nFrames; k += stride) picks.push_back(cand[k]);
   printf("  run %s: %zu truncated candidates (>= %d clusters, kept < 75 %%), %zu frames\n",
          run.Data(), cand.size(), minCl, picks.size());

   gSystem->Unlink(outGif);
   auto *c = new TCanvas("c", "", 1400, 720);
   c->SetFillColor(kWhite);
   int frame = 0;
   for (auto &pk : picks) {
      tF->GetEntry(pk.entry); tB->GetEntry(pk.entry);
      auto *evF = (AtTrackingEvent*)aF->At(0); auto *evB = (AtTrackingEvent*)aB->At(0);
      if (!evF) continue;
      // GetTrackArray() returns BY VALUE: keep the copy alive, or trk dangles (garbage size and
      // stray clusters outside the chamber -- seen on the first version of this macro)
      auto trkArr = evF->GetTrackArray();
      AtTrack *trk = nullptr;
      for (auto &tr : trkArr) if (tr.GetTrackID() == pk.tid) { trk = &tr; break; }
      AtFittedTrack *ftF = nullptr, *ftB = nullptr;
      for (auto &ft : evF->GetFittedTracks()) if (ft && ft->GetTrackID() == pk.tid) ftF = ft.get();
      if (evB) for (auto &ft : evB->GetFittedTracks()) if (ft && ft->GetTrackID() == pk.tid) ftB = ft.get();
      if (!trk || !ftF) continue;
      auto *hc = trk->GetHitClusterArray();

      std::vector<XYZPoint> fitF, fitB;
      for (auto &p : ftF->GetSmoothedPositions()) fitF.emplace_back(p.X(), p.Y(), kZPadFit - p.Z());
      if (ftB) for (auto &p : ftB->GetSmoothedPositions()) fitB.emplace_back(p.X(), p.Y(), kZPadFit - p.Z());

      // clusters split by whether the find805 fit passes next to them
      std::vector<double> xi, yi, zi, xo, yo, zo;
      double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9, zmin = 1e9, zmax = -1e9;
      auto grow = [&](double x, double y, double z) {
         xmin = std::min(xmin, x); xmax = std::max(xmax, x); ymin = std::min(ymin, y);
         ymax = std::max(ymax, y); zmin = std::min(zmin, z); zmax = std::max(zmax, z); };
      for (auto &cl : *hc) {
         auto p = cl.GetPosition(); grow(p.X(), p.Y(), p.Z());
         double d2 = 1e18;
         for (auto &q : fitF) d2 = std::min(d2, (p - q).Mag2());
         if (d2 < kNearMM * kNearMM) { xi.push_back(p.X()); yi.push_back(p.Y()); zi.push_back(p.Z()); }
         else                        { xo.push_back(p.X()); yo.push_back(p.Y()); zo.push_back(p.Z()); }
      }
      for (auto &q : fitF) grow(q.X(), q.Y(), q.Z());
      for (auto &q : fitB) grow(q.X(), q.Y(), q.Z());
      for (auto *f : {ftF, ftB}) if (f) { auto v = f->GetVertex(); grow(v.X(), v.Y(), kZPadFit - v.Z()); }

      c->Clear();
      auto *top = new TPad("top", "", 0, 0.86, 1, 1); top->Draw();
      auto *pl = new TPad("pl", "", 0, 0, 0.42, 0.86); pl->Draw();
      auto *pr = new TPad("pr", "", 0.42, 0, 1, 0.86); pr->Draw();

      auto &kF = ftF->GetKinematicsXtr();
      const int n = hc->size(), np = fitF.size();
      top->cd();
      TLatex tx; tx.SetNDC(); tx.SetTextFont(42);
      tx.SetTextSize(0.24);
      tx.DrawLatex(0.01, 0.62, TString::Format("#bf{16C(d,p)  run_%s  entry %lld  track %d}   %d clusters   "
                   "kept #geq %.0f%% (%d fit pts)   [%d/%zu]",
                   run.Data(), pk.entry, pk.tid, n, 100.0 * np / n, np, frame + 1, picks.size()));
      tx.SetTextColor(kRed + 1);
      // find805 ran with c2Max 0.1: above it, no rung passed and the minimum-chi2 rung was kept
      tx.DrawLatex(0.01, 0.34, TString::Format("find805 (adopted):  KE %.3f MeV   #theta %.1f#circ   #chi^{2}/ndf %.3f%s",
                   kF.kineticEnergy, kF.theta * TMath::RadToDeg(), c2n(ftF),
                   c2n(ftF) > 0.1 ? "   #bf{no rung passed 0.1: min-#chi^{2} fallback}" : ""));
      tx.SetTextColor(kAzure + 2);
      if (ftB) { auto &kB = ftB->GetKinematicsXtr();
         tx.DrawLatex(0.01, 0.06, TString::Format("full length (cateloss):  KE %.3f MeV   #theta %.1f#circ   #chi^{2}/ndf %.3f",
                      kB.kineticEnergy, kB.theta * TMath::RadToDeg(), c2n(ftB))); }
      else tx.DrawLatex(0.01, 0.06, "full length (cateloss):  no fit for this track");

      auto draw = [&](TPad *p, int ia, int ib, double lo1, double hi1, double lo2, double hi2,
                      const char *t1, const char *t2) {
         p->cd(); p->SetLeftMargin(0.13); p->SetBottomMargin(0.12); p->SetGrid();
         double m1 = 0.06 * (hi1 - lo1) + 5, m2 = 0.06 * (hi2 - lo2) + 5;
         auto *fr = p->DrawFrame(lo1 - m1, lo2 - m2, hi1 + m1, hi2 + m2,
                                 TString::Format(";%s (mm);%s (mm)", t1, t2));
         fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
         auto pick = [&](const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, int a) -> const std::vector<double>& {
            return a == 0 ? x : a == 1 ? y : z; };
         auto gr = [&](const std::vector<double> &a, const std::vector<double> &b, int col, int st, double sz) {
            if (a.empty()) return; auto *g = new TGraph(a.size(), a.data(), b.data());
            g->SetMarkerStyle(st); g->SetMarkerColor(col); g->SetMarkerSize(sz); g->Draw("P"); };
         gr(pick(xo, yo, zo, ia), pick(xo, yo, zo, ib), kGray, 20, 0.55);
         gr(pick(xi, yi, zi, ia), pick(xi, yi, zi, ib), kGray + 3, 20, 0.55);
         auto line = [&](const std::vector<XYZPoint> &v, int col, int ls) {
            if (v.size() < 2) return; std::vector<double> a, b;
            for (auto &q : v) { double c3[3] = {q.X(), q.Y(), q.Z()}; a.push_back(c3[ia]); b.push_back(c3[ib]); }
            auto *g = new TGraph(a.size(), a.data(), b.data());
            g->SetLineColor(col); g->SetLineWidth(3); g->SetLineStyle(ls); g->Draw("L"); };
         line(fitB, kAzure + 2, 1);
         line(fitF, kRed + 1, 1);
         auto vtx = [&](AtFittedTrack *f, int col, int st) {
            if (!f) return; auto v = f->GetVertex(); double c3[3] = {v.X(), v.Y(), kZPadFit - v.Z()};
            auto *m = new TMarker(c3[ia], c3[ib], st); m->SetMarkerColor(col); m->SetMarkerSize(2.4); m->Draw(); };
         vtx(ftB, kAzure + 2, 30); vtx(ftF, kRed + 1, 29);
      };
      // XY with equal scales
      double cx = 0.5 * (xmin + xmax), cy = 0.5 * (ymin + ymax), h = 0.5 * std::max(xmax - xmin, ymax - ymin);
      draw(pl, 0, 1, cx - h, cx + h, cy - h, cy + h, "x", "y");
      draw(pr, 2, 0, zmin, zmax, xmin, xmax, "z", "x");
      top->cd();
      auto *lg = new TLegend(0.73, 0.02, 0.995, 0.98); lg->SetTextSize(0.20); lg->SetBorderSize(0);
      lg->SetNColumns(2);
      auto *e1 = new TGraph(); e1->SetLineColor(kRed + 1); e1->SetLineWidth(3);
      auto *e2 = new TGraph(); e2->SetLineColor(kAzure + 2); e2->SetLineWidth(3);
      auto *e3 = new TGraph(); e3->SetMarkerStyle(20); e3->SetMarkerColor(kGray + 3);
      auto *e4 = new TGraph(); e4->SetMarkerStyle(20); e4->SetMarkerColor(kGray);
      lg->AddEntry(e1, "find805 fit", "l"); lg->AddEntry(e2, "full-length fit", "l");
      lg->AddEntry(e3, "kept prefix", "p"); lg->AddEntry(e4, "dropped", "p");
      lg->Draw();

      c->Modified(); c->Update();
      c->Print(outGif + TString::Format("+%d", delayCs));
      if (pngPrefix != "") c->Print(TString::Format("%s_%02d.png", pngPrefix.Data(), frame));
      ++frame;
   }
   c->Print(outGif + "++");
   printf("  wrote %s : %d frames\n", outGif.Data(), frame);
}
