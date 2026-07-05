/// @file resolution_test.C
/// @brief Resolution residual distributions for UKF vs GENFIT on a test channel
///        (branch 8 pi or branch 10 K). Truth is analytic (fixed |p|, theta=90 deg,
///        vertex (0,0,75) mm), so no sim file is needed. Four panels — Dp/p, Dtheta,
///        Dz, Dr(xy) — with both fitters overlaid and per-fitter sigma_IQR.
/// Run: root -b -q resolution_test.C                       // pions, output_digi_pi.root
///      root -b -q 'resolution_test.C("K",0.777,"./data/output_digi_K.root")'
double iqr(std::vector<double> v)
{
   if (v.size() < 4) return 0;
   std::sort(v.begin(), v.end());
   return (v[3 * v.size() / 4] - v[v.size() / 4]) / 1.349;
}
double med(std::vector<double> v)
{
   if (v.empty()) return 0;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}

void resolution_test(TString species = "pi", Double_t testEnergy = -1,
                     TString digiFile = "./data/output_digi_pi.root", TString outPng = "./data/resolution_test.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const bool isK = (species == "K" || species == "kaon");
   const double mass = isK ? 493.677 : 139.57039;
   const double E0 = (testEnergy > 0) ? testEnergy : (isK ? 0.777 : 0.4);
   const double p0 = std::sqrt(E0 * E0 - (mass / 1000) * (mass / 1000)) * 1000; // MeV/c
   const double vz0 = 75.0;

   TFile f(digiFile); TTree *t = (TTree *)f.Get("cbmsim");
   if (!t) { printf("missing %s\n", digiFile.Data()); return; }
   TClonesArray *ukf = new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventUKF", &ukf);
   TClonesArray *gf = new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventGenfit", &gf);

   // residual vectors: [fitter][metric] 0=dp/p,1=dtheta,2=dz,3=dr
   std::vector<double> R[2][4];
   TClonesArray *arr[2] = {ukf, gf};
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      for (int fi = 0; fi < 2; ++fi) {
         if (!arr[fi]->GetEntries()) continue;
         for (const auto &ft : ((AtTrackingEvent *)arr[fi]->At(0))->GetFittedTracks()) {
            const auto &kin = ft->GetKinematics(0);
            double KE = kin.kineticEnergy;
            if (!(KE > 0)) continue;
            double p = std::sqrt(KE * KE + 2 * KE * mass);
            R[fi][0].push_back((p - p0) / p0);
            R[fi][1].push_back(kin.theta * 180.0 / M_PI - 90.0);
            const auto &pr = ft->GetTrackPropertiesStruct();
            double vx = pr.initialPositionXtr.X(), vy = pr.initialPositionXtr.Y(), vz = pr.initialPositionXtr.Z();
            if (std::abs(vx) < 1e-9 && std::abs(vy) < 1e-9 && std::abs(vz) < 1e-9) { const auto &v = ft->GetVertex(0); vx = v.X(); vy = v.Y(); vz = v.Z(); }
            R[fi][2].push_back(vz - vz0);
            R[fi][3].push_back(std::hypot(vx, vy));
         }
      }
   }

   // Dp/p is filled in PERCENT (scale 100), so the axis range is in percent too;
   // [-50,50]% covers ~94% of tracks. The heavy tail (near-straight tracks fit as
   // near-infinite radius) overflows and is reported as the >50% fraction below.
   struct Cfg { const char *title; double lo, hi; };
   Cfg cfg[4] = {{"#Deltap/p [%]", -50, 50}, {"#Delta#theta [deg]", -12, 12}, {"#Deltaz [mm]", -15, 15}, {"#Deltar(xy) [mm]", 0, 25}};
   const char *unit[4] = {"", " deg", " mm", " mm"};
   double scale[4] = {100, 1, 1, 1}; // dp/p shown in %

   auto *c = new TCanvas("res", "resolution UKF vs genfit", 1300, 950);
   c->Divide(2, 2);
   int col[2] = {kAzure + 2, kRed + 1};
   const char *fn[2] = {"UKF", "genfit"};
   for (int m = 0; m < 4; ++m) {
      c->cd(m + 1); gPad->SetGrid();
      TH1F *h[2];
      double ymax = 0;
      for (int fi = 0; fi < 2; ++fi) {
         h[fi] = new TH1F(Form("h%d_%d", m, fi), Form("%s;%s%s;tracks (norm.)", cfg[m].title, cfg[m].title, unit[m]), 80, cfg[m].lo, cfg[m].hi);
         for (double x : R[fi][m]) h[fi]->Fill((m == 0 ? scale[0] * x : x));
         if (h[fi]->Integral() > 0) h[fi]->Scale(1.0 / h[fi]->Integral());
         h[fi]->SetLineColor(col[fi]); h[fi]->SetLineWidth(2);
         ymax = std::max(ymax, h[fi]->GetMaximum());
      }
      h[0]->SetMaximum(ymax * 1.25); h[0]->Draw("HIST"); h[1]->Draw("HIST SAME");
      auto *leg = new TLegend(0.58, 0.72, 0.98, 0.90); leg->SetTextSize(0.032);
      for (int fi = 0; fi < 2; ++fi) {
         double s = iqr(R[fi][m]) * (m == 0 ? 100 : 1), md = med(R[fi][m]) * (m == 0 ? 100 : 1);
         leg->AddEntry(h[fi], Form("%s  med %+.2f  #sigma %.2f", fn[fi], md, s), "l");
      }
      leg->Draw();
   }
   c->SaveAs(outPng);
   printf("resolution test (%s, |p|=%.1f MeV/c) -> %s\n", species.Data(), p0, outPng.Data());
   for (int fi = 0; fi < 2; ++fi) {
      size_t n = R[fi][0].size(), tail = 0, core = 0;
      for (double x : R[fi][0]) { if (std::abs(x) > 0.5) tail++; if (std::abs(x) < 0.2) core++; }
      printf("  %-7s: dp/p sig(IQR) %.1f%%  [core |dp/p|<20%%: %.0f%%, tail >50%%: %.1f%%]  "
             "dtheta %.2f  dz %.2f mm  dr %.2f mm  (n=%zu)\n",
             fn[fi], 100 * iqr(R[fi][0]), 100.0 * core / n, 100.0 * tail / n,
             iqr(R[fi][1]), iqr(R[fi][2]), iqr(R[fi][3]), n);
   }
}
