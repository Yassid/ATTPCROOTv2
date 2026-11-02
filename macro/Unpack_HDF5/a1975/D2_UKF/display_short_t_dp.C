/// @file display_short_t_dp.C
/// @brief Event displays of the data tracks that the circle dE/dx PID calls TRITON-LIKE but whose range is SHORTER than
///        any triton at their circle Brho (2026-11-01), next to triton-like tracks that sit ON the triton range curve.
///        Rows: pid_sbpid_dp.C merged tree + range_dp.C (friend), restricted to one run; the row -> (event, track) map
///        is rebuilt by looping that run's reco file with the same rule (every >= 30-hit track of every event, in order).
///        Per track: xy (hits, saturated > 3500 ADC in red, vertex end marked), z vs xy path along the track, hit
///        amplitude vs path from the vertex end (3500 line), with circle Brho, Spyral Brho, range, ln R/R_t, saturation.
///   root -l -b -q 'display_short_t_dp.C("0026","plots/short_t_tracks.png")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
void display_short_t_dp(TString run, TString png, int nShort = 4, int nGood = 2, double bLo = 0.20, double bHi = 0.30,
                        TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root", TString rr = "/mnt/f/a1975/pid_range/R_all.root")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kWarning;
   gStyle->SetOptStat(0);
   namespace HX = AtTools::Helix3D;
   // sim triton median range vs circle Brho (range_dp.C on the triton sims)
   std::vector<double> bx, rt;
   {
      TChain t("w"), tr("r");
      for (const char *s : {"gs", "ex2"}) { t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); tr.Add(Form("/mnt/f/a1975/pid_range/r_%s.root", s)); }
      t.AddFriend(&tr, "r");
      for (double b = 0.10; b < 0.48; b += 0.02) {
         long n = t.Draw("r.len", Form("va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<955&&bA>%g&&bA<%g", b, b + 0.02), "goff");
         if (n < 8) continue;
         std::vector<double> v(t.GetV1(), t.GetV1() + n);
         std::sort(v.begin(), v.end());
         bx.push_back(b + 0.01); rt.push_back(v[n / 2]);
      }
   }
   auto rtAt = [&](double b) { for (size_t i = 0; i + 1 < bx.size(); ++i) if (b >= bx[i] && b < bx[i + 1]) return rt[i] + (rt[i + 1] - rt[i]) * (b - bx[i]) / (bx[i + 1] - bx[i]); return -1.; };
   // run row range in the merged tree (merge order 0016, 0026, 0058, 0088)
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   Long64_t off = 0, nrun = 0;
   for (int i = 0; i < 4; ++i) {
      TFile fr(Form("/mnt/f/a1975/pid3d_sbpid/pv_%s.root", runs[i]));
      Long64_t n = ((TTree *)fr.Get("w"))->GetEntries();
      if (run == runs[i]) { nrun = n; break; }
      off += n;
   }
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("r=r", rr);
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<955&&r.len>20";
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   TH1D hh("hh", "", 600, -4, 2);
   w->Draw(Form("log(s1_g/%g*bA)>>hh", md), sel + Form("&&bA>%g&&bA<%g", bLo, bHi), "goff");
   hh.Smooth(2);
   const double dm = hh.GetBinCenter(hh.GetMaximumBin());
   w->SetEstimate(w->GetEntries() + 1);
   long n = w->Draw("Entry$:bA:r.len:br:log(s1_g/" + TString::Format("%g", md) + "*bA)",
                    sel + Form("&&bA>%g&&bA<%g&&Entry$>=%lld&&Entry$<%lld", bLo, bHi, off, off + nrun), "goff");
   struct Pick { Long64_t row; double bA, len, br, lnm, lrt; };
   std::vector<Pick> shortT, goodT;
   for (long i = 0; i < n; ++i) {
      double lnm = w->GetVal(4)[i] - dm;
      if (lnm < 0.25 || lnm > 0.55) continue;
      double r0 = rtAt(w->GetV2()[i]);
      if (r0 <= 0) continue;
      Pick p{(Long64_t)w->GetV1()[i], w->GetV2()[i], w->GetV3()[i], w->GetV4()[i], lnm, std::log(w->GetV3()[i] / r0)};
      if (p.lrt < -0.3) shortT.push_back(p);
      else if (std::abs(p.lrt) < 0.15) goodT.push_back(p);
   }
   printf("run %s: triton-like in Brho %.2f-%.2f: %zu short, %zu on the triton curve\n", run.Data(), bLo, bHi, shortT.size(), goodT.size());
   std::vector<Pick> show;
   for (int i = 0; i < nShort && i < (int)shortT.size(); ++i) show.push_back(shortT[i * std::max<size_t>(1, shortT.size() / nShort)]);
   size_t nS = show.size();
   for (int i = 0; i < nGood && i < (int)goodT.size(); ++i) show.push_back(goodT[i * std::max<size_t>(1, goodT.size() / nGood)]);
   // row -> (event, track) in the run's reco file
   std::map<Long64_t, std::pair<Long64_t, int>> where;
   std::set<Long64_t> want;
   for (auto &p : show) want.insert(p.row - off);
   TFile fr(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", run.Data()));
   auto *t = (TTree *)fr.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   Long64_t row = 0;
   for (Long64_t ev = 0; ev < t->GetEntries() && where.size() < want.size(); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      int itr = 0;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr.GetHitArray().size() >= 30) { if (want.count(row)) where[row] = {ev, itr}; ++row; }
         ++itr;
      }
   }
   auto *c = new TCanvas("c", "c", 1500, 380 * show.size());
   c->Divide(3, show.size());
   HX::EndOpt eo;
   eo.hitRefit = true; eo.hitRefitZ = true;
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   for (size_t k = 0; k < show.size(); ++k) {
      auto &p = show[k];
      auto wh = where[p.row - off];
      t->GetEntry(wh.first);
      auto &tr = ((AtPatternEvent *)pa->At(0))->GetTrackCand()[wh.second];
      std::vector<HX::P> h;
      for (auto &q : tr.GetHitArray()) { auto ps = q->GetPosition(); h.push_back({ps.X(), ps.Y(), ps.Z(), q->GetCharge()}); }
      auto sr = spy.Estimate(tr);
      auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
      int spyEnd = (sr.valid && (sr.direction == 0 || sr.direction == 1)) ? ((sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1) : -1;
      auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
      // hits in order from the vertex end, xy path along that order
      std::vector<double> X, Y, Z, Q, S;
      std::vector<std::pair<double, int>> ord;
      for (size_t i = 0; i < r.order.size(); ++i) ord.push_back({V.end == 0 ? r.s[i] : r.length - r.s[i], r.order[i]});
      std::sort(ord.begin(), ord.end());
      long nsat = 0;
      for (auto &o : ord) { auto &q = h[o.second]; X.push_back(q.x); Y.push_back(q.y); Z.push_back(1000 - q.z); Q.push_back(q.q); S.push_back(o.first); nsat += q.q > 3500; }
      const char *tag = k < nS ? "SHORT" : "on t curve";
      c->cd(3 * k + 1);
      auto *g = new TGraph(X.size(), X.data(), Y.data());
      g->SetTitle(Form("%s: run %s ev %lld trk %d;x [mm];y [mm]", tag, run.Data(), wh.first, wh.second));
      g->SetMarkerStyle(20); g->SetMarkerSize(0.45); g->SetMarkerColor(kGray + 2); g->Draw("AP");
      std::vector<double> xs, ys;
      for (size_t i = 0; i < X.size(); ++i) if (Q[i] > 3500) { xs.push_back(X[i]); ys.push_back(Y[i]); }
      if (!xs.empty()) { auto *gs = new TGraph(xs.size(), xs.data(), ys.data()); gs->SetMarkerStyle(20); gs->SetMarkerSize(0.55); gs->SetMarkerColor(kRed); gs->Draw("P"); }
      auto *mv = new TMarker(X[0], Y[0], 29); mv->SetMarkerColor(kGreen + 2); mv->SetMarkerSize(2.2); mv->Draw();
      auto *m0 = new TMarker(0, 0, 5); m0->SetMarkerSize(2); m0->Draw();
      c->cd(3 * k + 2);
      auto *gz = new TGraph(S.size(), S.data(), Z.data());
      gz->SetTitle(Form("B#rho circle %.3f, Spyral %.3f;path from vertex end [mm];z_{lab} [mm]", p.bA, p.br));
      gz->SetMarkerStyle(20); gz->SetMarkerSize(0.45); gz->SetMarkerColor(kBlue + 1); gz->Draw("AP");
      c->cd(3 * k + 3);
      auto *gq = new TGraph(S.size(), S.data(), Q.data());
      gq->SetTitle(Form("range %.0f mm (ln R/R_{t} %+.2f), ln-mass %+.2f, %ld/%zu hits > 3500;path from vertex end [mm];hit amplitude [ADC]",
                        p.len, p.lrt, p.lnm, nsat, Q.size()));
      gq->SetMarkerStyle(20); gq->SetMarkerSize(0.45); gq->SetMarkerColor(kBlack); gq->Draw("AP");
      auto *l = new TLine(0, 3500, S.back(), 3500); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw();
      printf("  %-10s ev %lld trk %d: bA %.3f br %.3f range %.0f lnR/Rt %+.2f lnm %+.2f, %ld/%zu saturated hits, z extent %.0f mm\n", tag, wh.first, wh.second,
             p.bA, p.br, p.len, p.lrt, p.lnm, nsat, Q.size(), *std::max_element(Z.begin(), Z.end()) - *std::min_element(Z.begin(), Z.end()));
   }
   c->SaveAs(png);
}
