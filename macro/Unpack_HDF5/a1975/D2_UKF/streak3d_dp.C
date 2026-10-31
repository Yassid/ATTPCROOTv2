/// @file streak3d_dp.C
/// @brief The thin streak on the 3D PID plane at Brho 0.11-0.16, sqrt(dE/dx) > 42 (beyond the end of the
///        proton/deuteron band): which estimator is right on those tracks? Uses the fit-free radius
///        arbiter d_max (largest distance from the beam axis = the circle's DIAMETER for a track launched
///        from the axis; /0.93 for a stopping spiral, measured on (d,t)) against R_3D (AtHelix3D) and
///        R_Spyral (= Brho sin(theta)/B). Also draws the first `nShow` streak tracks.
///   root -l -b -q 'streak3d_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","plots/streak3d_dp.png")'
void streak3d_dp(TString reco, TString png, int nShow = 6, double brLo = 0.11, double brHi = 0.16, double sqMin = 42)
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
   AtTools::AtHelix3D h3;
   std::vector<double> r3S, rSS, r3B, rSB;
   int c3 = 0, cS = 0;
   std::vector<std::pair<Long64_t, int>> show;
   for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
      t->GetEntry(ev);
      if (!pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         if (H.size() < 30) continue;
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
         if (!sr.valid) continue;
         auto hr = h3.Estimate(tr, sr.direction);
         if (!hr.valid || hr.spacing >= 6 || hr.brho < brLo || hr.brho > brHi) continue;
         bool streak = hr.sqrtdEdx > sqMin, band = hr.sqrtdEdx > 28 && hr.sqrtdEdx < 38;
         if (!streak && !band) continue;
         std::vector<double> rho;
         for (auto &h : H) rho.push_back(std::hypot(h->GetPosition().X(), h->GetPosition().Y()));
         std::sort(rho.begin(), rho.end());
         double half = 0.5 * rho[rho.size() * 98 / 100] / 0.93;
         double r3 = hr.radius / half, rS = sr.brho * std::abs(std::sin(sr.polar)) / 2.85 * 1000 / half;
         if (streak) {
            r3S.push_back(r3); rSS.push_back(rS);
            (std::abs(std::log(r3)) < std::abs(std::log(rS)) ? c3 : cS)++;
            if ((int)show.size() < nShow) show.push_back({ev, tr.GetTrackID()});
         } else { r3B.push_back(r3); rSB.push_back(rS); }
      }
   }
   auto q = [](const char *l, std::vector<double> v) {
      std::sort(v.begin(), v.end());
      int n = v.size();
      if (n > 5) printf("  %-30s n=%5d p16 %.2f med %.2f p84 %.2f\n", l, n, v[n * 16 / 100], v[n / 2], v[n * 84 / 100]);
   };
   printf("radius / (d_max/2/0.93)   (1.00 = consistent with the track's own extent)\n");
   q("STREAK: R_3D", r3S); q("STREAK: R_Spyral", rSS); q("band (same Brho): R_3D", r3B); q("band (same Brho): R_Spyral", rSB);
   printf("STREAK tracks where 3D is closer to d_max: %d, Spyral closer: %d\n", c3, cS);
   // draw
   auto *c = new TCanvas("c", "c", 1800, 1200);
   c->Divide(3, 2);
   for (size_t k = 0; k < show.size(); ++k) {
      t->GetEntry(show[k].first);
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr.GetTrackID() != show[k].second) continue;
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
         auto hr = h3.Estimate(tr, sr.direction);
         double xl = 1e9, xh = -1e9, yl = 1e9, yh = -1e9;
         auto *g = new TGraph;
         for (auto &h : tr.GetHitArray()) {
            auto p = h->GetPosition();
            g->AddPoint(p.X(), p.Y());
            xl = std::min(xl, p.X()); xh = std::max(xh, p.X()); yl = std::min(yl, p.Y()); yh = std::max(yh, p.Y());
         }
         double cx = 0.5 * (xl + xh), cy = 0.5 * (yl + yh), hw = 0.6 * std::max({xh - xl, yh - yl, 60.0});
         c->cd(k + 1);
         gPad->DrawFrame(cx - hw, cy - hw, cx + hw, cy + hw,
                         Form("ev %lld tid %d n=%zu: B#rho 3D %.3f / Spy %.3f, #sqrt{dE/dx} 3D %.1f / Spy %.1f;x [mm];y [mm]",
                              show[k].first, show[k].second, tr.GetHitArray().size(), hr.brho, sr.brho, hr.sqrtdEdx,
                              sr.sqrtdEdx));
         g->SetMarkerStyle(20); g->SetMarkerSize(0.45); g->Draw("P");
         // 3D vertex-end circle (radius and closest approach known; centre reconstructed from the vertex)
         auto *e3 = new TEllipse(sr.center.X(), sr.center.Y(), sr.radius, sr.radius);
         e3->SetFillStyle(0); e3->SetLineColor(kGreen + 2); e3->SetLineStyle(2); e3->SetLineWidth(2); e3->Draw();
         auto *ax = new TMarker(0, 0, 5); ax->SetMarkerSize(2); ax->Draw();
         auto *vx = new TMarker(hr.vertex.X(), hr.vertex.Y(), 29); vx->SetMarkerColor(kRed); vx->SetMarkerSize(2.2); vx->Draw();
         auto *tl = new TLatex(0.13, 0.85, Form("R_{3D} %.0f mm (red star = 3D vertex), R_{Spy} %.0f mm (green)", hr.radius, sr.radius));
         tl->SetNDC(); tl->SetTextSize(0.035); tl->Draw();
      }
   }
   c->SaveAs(png);
}
