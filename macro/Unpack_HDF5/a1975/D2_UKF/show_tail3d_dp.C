/// @file show_tail3d_dp.C
/// @brief Draw 3D-only proton-band deuterons (tail3d_dp.C list): hits coloured along the MST order, the
///        vertex-end helix circle used for Brho (red), the circle Spyral's Brho implies through the
///        same vertex region (green dashed), the end windows. PNG only.
#include "../UKF/pid/helix3d.h"
void show_tail3d_dp(TString png, std::vector<std::pair<long, int>> L = {{124, 1}, {199, 0}, {260, 0}, {886, 0}, {1206, 1}, {1448, 0}},
                    TString reco = "/mnt/f/a1975/reco_d2_tau/run_0016_multifit_reco.root")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gErrorIgnoreLevel = kError;
   gStyle->SetOptStat(0);
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   auto *c = new TCanvas("c", "c", 1800, 1200);
   c->Divide(3, 2);
   for (size_t k = 0; k < L.size() && k < 6; ++k) {
      t->GetEntry(L[k].first);
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr.GetTrackID() != L[k].second) continue;
         std::vector<h3d::P> h;
         for (auto &x : tr.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto r = h3d::estimate(h, 2.85, -1, 1000);
         auto sr = spy.Estimate(tr);
         int spyEnd = -1;
         if (sr.valid) spyEnd = (sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1;
         auto V = h3d::majorityEnd(r, spyEnd);
         auto &E = r.e[V.end];
         double xl = 1e9, xh = -1e9, yl = 1e9, yh = -1e9;
         for (auto &p : h) { xl = std::min(xl, p.x); xh = std::max(xh, p.x); yl = std::min(yl, p.y); yh = std::max(yh, p.y); }
         double cx = 0.5 * (xl + xh), cy = 0.5 * (yl + yh), hw = 0.6 * std::max({xh - xl, yh - yl, 2.2 * E.R, 60.0});
         c->cd(k + 1);
         gPad->DrawFrame(cx - hw, cy - hw, cx + hw, cy + hw,
                         Form("ev %ld tid %d n=%zu (-%d blob): R_{3D} %.0f, B#rho 3D %.3f / Spy %.3f, dca %.0f, votes %d/%d%s;x [mm];y [mm]",
                              L[k].first, L[k].second, h.size(), r.nDropped, E.R, E.brho, sr.valid ? sr.brho : -1, E.dca,
                              V.nFor, V.nVotes, V.ambiguous ? " AMB" : ""));
         // hits dropped as blobs (not in r.order) in grey
         {
            std::vector<char> used(h.size(), 0);
            for (int i : r.order) used[i] = 1;
            for (size_t i = 0; i < h.size(); ++i)
               if (!used[i]) { auto *m = new TMarker(h[i].x, h[i].y, 20); m->SetMarkerSize(0.45); m->SetMarkerColor(kGray + 1); m->Draw(); }
         }
         const int n = r.order.size();
         for (int i = 0; i < n; ++i) {
            auto &p = h[r.order[i]];
            double fr = double(i) / n;
            auto *m = new TMarker(p.x, p.y, 20);
            m->SetMarkerSize(0.45);
            m->SetMarkerColor(TColor::GetColor(int(30 + 220 * fr), 40, int(250 - 220 * fr)));
            m->Draw();
         }
         auto *e = new TEllipse(E.xc, E.yc, E.R, E.R);
         e->SetFillStyle(0); e->SetLineColor(kRed); e->SetLineWidth(2); e->Draw();
         if (sr.valid) { // Spyral's circle: its own centre and the radius its Brho implies at its polar
            double Rs = sr.radius;
            auto *es = new TEllipse(sr.center.X(), sr.center.Y(), Rs, Rs);
            es->SetFillStyle(0); es->SetLineColor(kGreen + 2); es->SetLineStyle(2); es->SetLineWidth(2); es->Draw();
         }
         auto &cv = V.end == 0 ? r.cen.front() : r.cen.back();
         auto *mv = new TMarker(cv.x, cv.y, 29); mv->SetMarkerColor(kOrange + 7); mv->SetMarkerSize(2.5); mv->Draw();
         auto *ax = new TMarker(0, 0, 5); ax->SetMarkerSize(2); ax->Draw();
      }
   }
   c->SaveAs(png);
}
