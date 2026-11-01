/// @file display_theta_dp.C
/// @brief Event displays of the data tracks where the 3D arc-length polar angle (z vs xy arc length round the
///        AtHelix3D hit-refit circle, first 50 % of the vertex-end window) and Spyral's polar angle disagree most
///        (2026-11-01). Per track, three panels: xy (hits, hit-refit circle, vertex-end hit, beam axis); z vs arc
///        length (window hits, the fitted line, the slope Spyral's angle implies); z vs xy distance from the
///        vertex-end hit (Spyral's view). Pass 1 ranks the tracks, pass 2 draws the nShow worst -- optionally only
///        tracks with score in [sLo, sHi] and |cos theta_Spyral| <= maxCosS (to skip Spyral's along-the-beam failures).
///   root -l -b -q 'display_theta_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0026_multifit_reco.root","plots/theta_worst.png")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
struct ThetaTrack {
   bool ok{false};
   double score{0}, thA{0}, thS{0}, brho{0}, slope{0}, icpt{0};
   std::vector<double> x, y, z, q, s, rho; // window hits in order from the vertex end (lab z), arc length, xy distance
   int nFit{0};
   double xc{0}, yc{0}, R{0};
   std::vector<double> ax, ay, az; // all hits of the track (lab z)
};
ThetaTrack analyse(AtTrack &tr, AtTools::AtSpyralPID &spy)
{
   namespace HX = AtTools::Helix3D;
   ThetaTrack T;
   auto sr = spy.Estimate(tr);
   if (!sr.valid) return T;
   std::vector<HX::P> h;
   for (auto &x : tr.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
   HX::EndOpt eo;
   eo.hitRefit = true;
   auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
   if (!r.ok || r.length / std::max(1, r.nBackbone) >= 6) return T;
   const double zA = h[r.order.front()].z, zB = h[r.order.back()].z;
   int spyEnd = (sr.direction == 0 || sr.direction == 1) ? ((sr.direction == 0) == (zA < zB) ? 0 : 1) : -1;
   auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
   const auto &E = r.e[V.end];
   if (!E.ok || E.brho <= 0) return T;
   const int m = r.cen.size();
   const double lim = (V.end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[m - E.nUsed]) + 2.0;
   std::vector<std::pair<double, int>> wv;
   for (size_t i = 0; i < r.order.size(); ++i) {
      double sf = V.end == 0 ? r.s[i] : r.length - r.s[i];
      if (sf <= lim) wv.push_back({sf, r.order[i]});
   }
   std::sort(wv.begin(), wv.end());
   if (wv.size() < 12) return T;
   T.xc = E.xc; T.yc = E.yc; T.R = E.R;
   double prev = std::atan2(h[wv[0].second].y - E.yc, h[wv[0].second].x - E.xc), ph = 0;
   const auto &h0 = h[wv[0].second];
   for (size_t i = 0; i < wv.size(); ++i) {
      const auto &hp = h[wv[i].second];
      double a = std::atan2(hp.y - E.yc, hp.x - E.xc), d = a - prev;
      while (d > M_PI) d -= 2 * M_PI;
      while (d < -M_PI) d += 2 * M_PI;
      if (i) ph += d;
      prev = a;
      T.x.push_back(hp.x); T.y.push_back(hp.y); T.z.push_back(1000.0 - hp.z); T.q.push_back(hp.q);
      T.s.push_back(E.R * std::abs(ph)); T.rho.push_back(std::hypot(hp.x - h0.x, hp.y - h0.y));
   }
   T.nFit = std::min<int>(wv.size(), std::max(6, (int)(0.5 * wv.size())));
   double sx = 0, sy = 0, sxx = 0, sxy = 0;
   for (int i = 0; i < T.nFit; ++i) { sx += T.s[i]; sy += T.z[i]; sxx += T.s[i] * T.s[i]; sxy += T.s[i] * T.z[i]; }
   double den = T.nFit * sxx - sx * sx;
   if (std::abs(den) < 1e-9) return T;
   T.slope = (T.nFit * sxy - sx * sy) / den;
   T.icpt = (sy - T.slope * sx) / T.nFit;
   T.thA = std::atan2(1.0, T.slope);
   T.thS = sr.polar;
   T.brho = E.brho;
   for (auto &x : h) { T.ax.push_back(x.x); T.ay.push_back(x.y); T.az.push_back(1000.0 - x.z); }
   T.score = std::abs(std::log(std::sin(T.thA) / std::abs(std::sin(T.thS))));
   T.ok = true;
   return T;
}
void display_theta_dp(TString in, TString png, int nShow = 6, double bLo = 0.45, double bHi = 1.0, Long64_t maxEv = 1000000000,
                      double sLo = 0, double sHi = 1e9, double maxCosS = 1.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   // pass 1: rank
   struct Cand { double score; Long64_t ev; int itr; };
   std::vector<Cand> cand;
   std::vector<double> allScores;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      int itr = 0;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr.GetHitArray().size() >= 30) {
            auto T = analyse(tr, spy);
            if (T.ok && T.brho > bLo && T.brho < bHi) {
               allScores.push_back(T.score);
               if (T.score >= sLo && T.score <= sHi && std::abs(std::cos(T.thS)) <= maxCosS) cand.push_back({T.score, ev, itr});
            }
         }
         ++itr;
      }
   }
   std::sort(cand.begin(), cand.end(), [](const Cand &a, const Cand &b) { return a.score > b.score; });
   std::sort(allScores.begin(), allScores.end());
   printf("%zu tracks in Brho %.2f-%.2f; |ln(sin thA/sin thS)| median %.4f, p84 %.4f, p98 %.4f\n", allScores.size(), bLo, bHi,
          allScores[allScores.size() / 2], allScores[allScores.size() * 84 / 100], allScores[allScores.size() * 98 / 100]);
   // pass 2: draw
   auto *c = new TCanvas("c", "c", 1500, 420 * nShow);
   c->Divide(3, nShow);
   for (int k = 0; k < nShow && k < (int)cand.size(); ++k) {
      t->GetEntry(cand[k].ev);
      auto &tr = ((AtPatternEvent *)pa->At(0))->GetTrackCand()[cand[k].itr];
      auto T = analyse(tr, spy);
      printf("  #%d ev %lld track %d: Brho %.3f  theta arc50 %.1f deg  theta Spyral %.1f deg (sin ratio %+.3f), %zu window hits\n", k,
             cand[k].ev, cand[k].itr, T.brho, T.thA * TMath::RadToDeg(), T.thS * TMath::RadToDeg(), T.score, T.s.size());
      // xy
      c->cd(3 * k + 1);
      auto *gA = new TGraph(T.ax.size(), T.ax.data(), T.ay.data());
      gA->SetTitle(Form("ev %lld trk %d, B#rho %.2f;x [mm];y [mm]", cand[k].ev, cand[k].itr, T.brho));
      gA->SetMarkerStyle(20); gA->SetMarkerSize(0.4); gA->SetMarkerColor(kGray + 1);
      gA->Draw("AP");
      auto *gW = new TGraph(T.x.size(), T.x.data(), T.y.data());
      gW->SetMarkerStyle(20); gW->SetMarkerSize(0.5); gW->SetMarkerColor(kBlue + 1); gW->Draw("P");
      auto *el = new TEllipse(T.xc, T.yc, T.R, T.R); el->SetFillStyle(0); el->SetLineColor(kRed); el->Draw();
      auto *v = new TMarker(T.x[0], T.y[0], 29); v->SetMarkerColor(kGreen + 2); v->SetMarkerSize(2.2); v->Draw();
      auto *ax0 = new TMarker(0, 0, 5); ax0->SetMarkerSize(2); ax0->Draw();
      // z vs arc length
      c->cd(3 * k + 2);
      auto *gS = new TGraph(T.s.size(), T.s.data(), T.z.data());
      gS->SetTitle(Form("#theta_{arc50} %.1f#circ vs #theta_{Spyral} %.1f#circ;xy arc length from vertex end [mm];z_{lab} [mm]",
                        T.thA * TMath::RadToDeg(), T.thS * TMath::RadToDeg()));
      gS->SetMarkerStyle(20); gS->SetMarkerSize(0.5); gS->SetMarkerColor(kBlue + 1);
      gS->Draw("AP");
      std::vector<double> sf(T.s.begin(), T.s.begin() + T.nFit), zf(T.z.begin(), T.z.begin() + T.nFit);
      auto *gF = new TGraph(sf.size(), sf.data(), zf.data());
      gF->SetMarkerStyle(24); gF->SetMarkerSize(0.8); gF->SetMarkerColor(kBlack); gF->Draw("P");
      double smax = *std::max_element(T.s.begin(), T.s.end());
      auto *lA = new TLine(0, T.icpt, smax, T.icpt + T.slope * smax); lA->SetLineColor(kRed); lA->SetLineWidth(2); lA->Draw();
      // Spyral slope |cot theta_S| with the arc fit's sign, through the mean of the fitted hits
      double cS = std::abs(std::cos(T.thS) / std::sin(T.thS)) * (T.slope >= 0 ? 1 : -1);
      double ms = 0, mz = 0;
      for (int i = 0; i < T.nFit; ++i) { ms += T.s[i]; mz += T.z[i]; }
      ms /= T.nFit; mz /= T.nFit;
      auto *lS = new TLine(0, mz - cS * ms, smax, mz + cS * (smax - ms)); lS->SetLineColor(kGreen + 2); lS->SetLineWidth(2);
      lS->SetLineStyle(2); lS->Draw();
      // z vs xy distance from the vertex-end hit (Spyral's view)
      c->cd(3 * k + 3);
      auto *gR = new TGraph(T.rho.size(), T.rho.data(), T.z.data());
      gR->SetTitle(";xy distance from vertex-end hit [mm];z_{lab} [mm]");
      gR->SetMarkerStyle(20); gR->SetMarkerSize(0.5); gR->SetMarkerColor(kMagenta + 2);
      gR->Draw("AP");
   }
   c->SaveAs(png);
}
