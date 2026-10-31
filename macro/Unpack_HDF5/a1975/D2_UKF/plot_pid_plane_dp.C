/// @file plot_pid_plane_dp.C
/// @brief The a1975 D2 PID plane (sqrt(dE/dx) vs Brho), all runs merged, four ways: Spyral and the 3D
///        estimator (AtHelix3D), each before and after the per-run gain correction (AtPIDGainTable).
///        The production proton gate is drawn on the uncorrected Spyral plane, the only one it belongs
///        to. Axes are framed on the data (the 3D dE/dx scale differs from Spyral's).
///        Input: PID trees from pid_variants_dp.C, <prefix><run>.root (sqa/br Spyral, s0/b0 3D).
///   root -l -b -q 'plot_pid_plane_dp.C("/mnt/f/a1975/pid3d_gaincal/pv_","pid/gain_d2_dv1104.txt","plots/pid_plane_dp.png")'
void plot_pid_plane_dp(TString prefix, TString gainTable, TString png, TString gate = "pid/proton_d2_dv1104.json",
                       TString sel = "va&&v3&&sp3<6&&nh>=30")
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   AtTools::AtPIDGainTable gt;
   if (!gt.Load(gainTable.Data())) { printf("ERROR: cannot load %s\n", gainTable.Data()); return; }
   // runs present in the prefix directory
   TString dir = gSystem->GetDirName(prefix), base = gSystem->BaseName(prefix);
   std::vector<TString> runs;
   void *d = gSystem->OpenDirectory(dir);
   while (const char *e = gSystem->GetDirEntry(d)) {
      TString s(e);
      if (s.BeginsWith(base) && s.EndsWith(".root")) {
         TString r = s(base.Length(), s.Length() - base.Length() - 5);
         if (r.IsDigit()) runs.push_back(r);
      }
   }
   gSystem->FreeDirectory(d);
   std::sort(runs.begin(), runs.end());
   // pass 1: data-framed x ranges (99th percentile of each estimator's sqrt(dE/dx))
   std::vector<double> xs, x3;
   for (auto &r : runs) {
      TFile f(prefix + r + ".root");
      auto *w = (TTree *)f.Get("w");
      w->Draw("sqa:s0", sel, "goff");
      for (int i = 0; i < w->GetSelectedRows(); i += 5) { xs.push_back(w->GetV1()[i]); x3.push_back(w->GetV2()[i]); }
   }
   std::sort(xs.begin(), xs.end());
   std::sort(x3.begin(), x3.end());
   const double xmS = 1.15 * xs[xs.size() * 99 / 100], xm3 = 1.15 * x3[x3.size() * 99 / 100];
   const char *ttl[4] = {"Spyral (production), no gain correction", "Spyral, per-run gain corrected",
                         "3D (AtHelix3D), no gain correction", "3D (AtHelix3D), per-run gain corrected"};
   TH2D *h[4];
   for (int k = 0; k < 4; ++k)
      h[k] = new TH2D(Form("h%d", k), Form("%s;#sqrt{dE/dx} [arb.];B#rho [T m]", ttl[k]), 300, 0, k < 2 ? xmS : xm3, 300, 0, 1.2);
   long nTr = 0;
   for (auto &r : runs) {
      int run = r.Atoi();
      if (!gt.Has(run)) { printf("WARNING: run %d not in the gain table -- skipped\n", run); continue; }
      double gS = gt.GetSpyral(run), g3 = gt.Get3D(run);
      TFile f(prefix + r + ".root");
      auto *w = (TTree *)f.Get("w");
      w->Draw("sqa:br:s0:b0", sel, "goff para");
      for (int i = 0; i < w->GetSelectedRows(); ++i) {
         double sq = w->GetVal(0)[i], br = w->GetVal(1)[i], s3 = w->GetVal(2)[i], b3 = w->GetVal(3)[i];
         h[0]->Fill(sq, br); h[1]->Fill(sq * gS, br); h[2]->Fill(s3, b3); h[3]->Fill(s3 * g3, b3);
         ++nTr;
      }
   }
   printf("%zu runs, %ld tracks\n", runs.size(), nTr);
   auto *c = new TCanvas("c", "c", 1800, 1500);
   c->Divide(2, 2);
   const int pad[4] = {1, 2, 3, 4};
   for (int k = 0; k < 4; ++k) {
      c->cd(pad[k]);
      gPad->SetLogz();
      gPad->SetRightMargin(0.12);
      h[k]->Draw("colz");
   }
   // production gate on the uncorrected Spyral plane (x = sqrtdedx, y = brho)
   if (!gate.IsNull() && !gSystem->AccessPathName(gate)) {
      auto g = AtTools::AtParticleID::LoadJSON(gate.Data());
      auto *gr = new TGraph;
      for (auto &[x, y] : g.GetCut().GetVertices()) gr->AddPoint(x, y);
      gr->SetLineColor(kRed);
      gr->SetLineWidth(2);
      c->cd(1);
      gr->Draw("L");
      auto *tx = new TLatex(0.45, 0.85, Form("production gate %s", gSystem->BaseName(gate)));
      tx->SetNDC(); tx->SetTextColor(kRed); tx->SetTextSize(0.03); tx->Draw();
   }
   c->SaveAs(png);
}
