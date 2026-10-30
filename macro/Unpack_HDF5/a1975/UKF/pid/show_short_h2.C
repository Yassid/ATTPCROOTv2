/// @file show_short_h2.C
/// @brief Draw short H2 data tracks from the dir3d bench (n < 150, curling, not beam-parallel), so we
///        can see what the tracks the forward-only referee fails on actually are. PNG only, no GUI.
///
/// Per track two pads: x-y and rho vs z_lab (z_lab = 1000 - z_digi, measured on the sim). The track
/// is coloured along its 3D backbone order (blue = end A, red = end B), the event's other tracks
/// are grey, green star = the end the forward-only referee calls the vertex (lower z_lab), the
/// dashed circles are the local end helices. Title: hits, GeoTheta, Spyral call, 3D sense call.
///
///   root -l -b -q 'pid/show_short_h2.C(1,"pid/plots/short_h2_atpra_wrong.png")'  // AtPRA disagrees
///   root -l -b -q 'pid/show_short_h2.C(0,"pid/plots/short_h2_atpra_right.png")'  // AtPRA agrees
#include "helix3d.h"

void show_short_h2(int wantWrong = 1, TString png = "pid/plots/short_h2.png", int nShow = 8, int skip = 0,
                   TString bench = "/mnt/f/a1975/dir3d/d3_h2_0106.root",
                   TString reco = "/mnt/f/a1975/reco/run_0106_reco.root")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kError;
   gStyle->SetOptStat(0);
   gStyle->SetTitleFontSize(0.045);
   TFile fb(bench);
   auto *d = (TTree *)fb.Get("d3");
   TString E = "ok&&rhoMax>60&&abs(ezd[0]-ezd[1])>30&&th[0]>15&&th[0]<165&&th[1]>15&&th[1]<165&&n<150";
   TString atpraWrong = "((geoTh>90?(ezd[0]<ezd[1]?0:1):(ezd[0]<ezd[1]?1:0))!=(ezd[0]>ezd[1]?0:1))";
   TString sel = E + "&&" + (wantWrong ? atpraWrong : "!" + atpraWrong);
   d->Draw("ev:tid", sel, "goff");
   int nSel = d->GetSelectedRows();
   std::vector<std::pair<Long64_t, int>> L;
   for (int i = 0; i < nSel; ++i) L.push_back({(Long64_t)d->GetV1()[i], (int)d->GetV2()[i]});
   printf("%d tracks selected (%s)\n", nSel, wantWrong ? "AtPRA vs referee DISAGREE" : "agree");
   int stride = std::max(1, (nSel - skip) / nShow);

   TFile fr(reco);
   auto *t = (TTree *)fr.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   auto *c = new TCanvas("c", "c", 1800, 450 * ((nShow + 1) / 2));
   c->Divide(4, (nShow + 1) / 2);
   for (int k = 0; k < nShow && skip + k * stride < nSel; ++k) {
      auto [ev, tid] = L[skip + k * stride];
      t->GetEntry(ev);
      auto *pe = (AtPatternEvent *)pa->At(0);
      auto oxy = new TGraph, orz = new TGraph;
      oxy->SetMarkerStyle(7); oxy->SetMarkerColor(kGray + 1);
      orz->SetMarkerStyle(7); orz->SetMarkerColor(kGray + 1);
      AtTrack *me = nullptr;
      for (auto &tr : pe->GetTrackCand()) {
         if (tr.GetTrackID() == tid) { me = &tr; continue; }
         for (auto &h : tr.GetHitArray()) {
            auto p = h->GetPosition();
            oxy->AddPoint(p.X(), p.Y());
            orz->AddPoint(1000 - p.Z(), std::hypot(p.X(), p.Y()));
         }
      }
      if (!me) continue;
      std::vector<h3d::P> h;
      for (auto &x : me->GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
      auto res = h3d::estimate(h, 2.85, +1, 0);
      auto sr = spy.Estimate(*me);
      const int n = res.order.size();
      // colour along the backbone: blue (end A) -> red (end B)
      std::vector<TGraph *> gxy(5), grz(5);
      for (int s = 0; s < 5; ++s) {
         gxy[s] = new TGraph; grz[s] = new TGraph;
         int col = TColor::GetColor(int(40 + 50 * s), 60, int(240 - 50 * s));
         gxy[s]->SetMarkerStyle(20); gxy[s]->SetMarkerSize(0.45); gxy[s]->SetMarkerColor(col);
         grz[s]->SetMarkerStyle(20); grz[s]->SetMarkerSize(0.45); grz[s]->SetMarkerColor(col);
      }
      double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9, zlo = 1e9, zhi = -1e9, rhi = 0;
      for (int i = 0; i < n; ++i) {
         auto &p = h[res.order[i]];
         int s = std::min(4, 5 * i / std::max(1, n));
         double zl = 1000 - p.z, r = std::hypot(p.x, p.y);
         gxy[s]->AddPoint(p.x, p.y); grz[s]->AddPoint(zl, r);
         xlo = std::min(xlo, p.x); xhi = std::max(xhi, p.x); ylo = std::min(ylo, p.y); yhi = std::max(yhi, p.y);
         zlo = std::min(zlo, zl); zhi = std::max(zhi, zl); rhi = std::max(rhi, r);
      }
      // referee end: lower z_lab = higher z_digi
      int vRef = h[res.order.front()].z > h[res.order.back()].z ? 0 : 1;
      auto &pv = h[vRef == 0 ? res.order.front() : res.order.back()];
      int vSense = res.e[0].sense == res.e[1].sense ? -1 : (res.e[0].sense == 1 ? 0 : 1);
      const char *spyTxt = !sr.valid ? Form("fail%d", sr.failCode) : (sr.direction == 0 ? "fwd" : "bwd");
      // frame on the track itself, square in x-y, with the axis included
      double cx = 0.5 * (xlo + xhi), cy = 0.5 * (ylo + yhi), half = 0.6 * std::max({xhi - xlo, yhi - ylo, 40.0});
      double fx0 = std::min(cx - half, -15.0), fx1 = std::max(cx + half, 15.0);
      double fy0 = std::min(cy - half, -15.0), fy1 = std::max(cy + half, 15.0);
      c->cd(2 * k + 1);
      gPad->DrawFrame(fx0, fy0, fx1, fy1,
                      Form("ev %lld tid %d  n=%d  GeoTh %.0f  Spyral %s  3D-sense end %s;x [mm];y [mm]", ev, tid,
                           (int)h.size(), me->GetGeoTheta() * TMath::RadToDeg(), spyTxt,
                           vSense < 0 ? "?" : (vSense == vRef ? "=ref" : "!=ref")));
      oxy->Draw("P");
      for (int kk = 0; kk < 2; ++kk) {
         auto &Ek = res.e[kk];
         if (!Ek.ok) continue;
         auto el = new TEllipse(Ek.xc, Ek.yc, Ek.R, Ek.R);
         el->SetFillStyle(0); el->SetLineStyle(2); el->SetLineColor(kk == 0 ? kBlue : kRed); el->Draw();
      }
      for (auto *g : gxy) g->Draw("P");
      auto ax = new TMarker(0, 0, 5); ax->SetMarkerSize(1.6); ax->Draw();
      auto mv = new TMarker(pv.x, pv.y, 29); mv->SetMarkerColor(kGreen + 2); mv->SetMarkerSize(2.2); mv->Draw();
      c->cd(2 * k + 2);
      gPad->DrawFrame(zlo - 30, 0, zhi + 30, rhi * 1.15 + 5,
                      Form("dca A %.0f / B %.0f mm;z_{lab} [mm];#rho [mm]", res.e[0].dca, res.e[1].dca));
      orz->Draw("P");
      for (auto *g : grz) g->Draw("P");
      auto mz = new TMarker(1000 - pv.z, std::hypot(pv.x, pv.y), 29);
      mz->SetMarkerColor(kGreen + 2); mz->SetMarkerSize(2.2); mz->Draw();
   }
   c->SaveAs(png);
}
