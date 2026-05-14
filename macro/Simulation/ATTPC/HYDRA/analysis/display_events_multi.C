/// @file display_events_multi.C
/// @brief 2x3 grid of top-down (xy) event displays with the FULL MC helix
/// drawn from the production vertex through the chamber. Useful for
/// visualising the GLAD-style wide-field setup where the pion curves
/// the whole way from target to detector.
///
/// The MC helix outside drift_volume is reconstructed analytically from
/// the production vertex position and momentum (for a uniform Bz field):
///   center = (x0 ± R sin φ_0, y0 ∓ R cos φ_0)
/// with R = p_MeV / (0.3·B_T) mm; sign chosen by charge (pi- curves CCW
/// in +Bz, pi+ curves CW).
///
/// Usage: root -b -q 'analysis/display_events_multi.C(800, 0)'
///   shows events [start, start+6) at p = P_MeV.

void display_events_multi(int P_MeV = 800, int startEvent = 0,
                          double B_T = 2.0,
                          const char *runDir = "data",
                          const char *outPng = "data/event_multi_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.055, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.135);
   gStyle->SetPadBottomMargin(0.135);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.10);

   TString simF = Form("%s/HYDRAsim_p%d.root", runDir, P_MeV);
   TString digiF = Form("%s/output_digi_p%d.root", runDir, P_MeV);
   TString ukfF  = Form("%s/output_ukf_HYDRA_p%d.root", runDir, P_MeV);

   TFile fS(simF), fD(digiF), fU(ukfF);
   auto *tS = (TTree *)fS.Get("cbmsim");
   auto *tD = (TTree *)fD.Get("cbmsim");
   auto *tU = (TTree *)fU.Get("cbmsim");
   auto *trks = new TClonesArray("AtMCTrack");
   auto *pts  = new TClonesArray("AtMCPoint");
   auto *ev   = new TClonesArray("AtEvent");
   auto *te   = new TClonesArray("AtTrackingEvent");
   tS->SetBranchAddress("MCTrack", &trks);
   tS->SetBranchAddress("AtTpcPoint", &pts);
   tD->SetBranchAddress("AtEventH", &ev);
   tU->SetBranchAddress("AtTrackingEvent", &te);

   auto *c = new TCanvas("c",
                         Form("HYDRA p = %d MeV/c, events %d..%d (xy)", P_MeV, startEvent,
                              startEvent + 5),
                         1700, 1100);
   c->Divide(3, 2, 0.012, 0.020);

   for (int panel = 0; panel < 6; ++panel) {
      int eventIdx = startEvent + panel;
      tS->GetEntry(eventIdx);
      tD->GetEntry(eventIdx);
      tU->GetEntry(eventIdx);
      if (trks->GetEntries() == 0) continue;

      auto *mc = (AtMCTrack *)trks->At(0);
      double mcVx = mc->GetStartX() * 10., mcVy = mc->GetStartY() * 10.;
      double mcPx = mc->GetPx(), mcPy = mc->GetPy(), mcPz = mc->GetPz();
      double pMC = std::sqrt(mcPx * mcPx + mcPy * mcPy + mcPz * mcPz) * 1000.; // MeV/c
      double phi_0 = std::atan2(mcPy, mcPx);

      // Collect MC points in drift_volume (primary)
      std::vector<double> mcX, mcY;
      double xFirst = 1e9, yFirst = 0, phMCfirst = 0;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (p->GetTrackID() != 0) continue;
         double x_mm = p->GetX() * 10., y_mm = p->GetY() * 10.;
         mcX.push_back(x_mm);
         mcY.push_back(y_mm);
         if (x_mm < xFirst) {
            xFirst = x_mm; yFirst = y_mm;
            phMCfirst = std::atan2(p->GetPy(), p->GetPx());
         }
      }

      // Analytic helix from MC vertex to first MC point inside drift_volume.
      // pi-: center 90° CCW from velocity, dψ/ds = +1/R.
      // pi+: center 90° CW  from velocity, dψ/ds = −1/R.
      int charge_sign = (mc->GetPdgCode() < 0) ? -1 : +1; // pi- = -211
      double R = (B_T > 0) ? pMC / (0.3 * B_T) : 0; // mm
      double cx, cy;
      if (charge_sign < 0) {
         cx = mcVx - R * std::sin(phi_0);
         cy = mcVy + R * std::cos(phi_0);
      } else {
         cx = mcVx + R * std::sin(phi_0);
         cy = mcVy - R * std::cos(phi_0);
      }
      double psi_start = std::atan2(mcVy - cy, mcVx - cx);
      double psi_end   = (xFirst < 1e8) ? std::atan2(yFirst - cy, xFirst - cx) : psi_start;
      if (charge_sign < 0)
         while (psi_end < psi_start) psi_end += 2 * M_PI;
      else
         while (psi_end > psi_start) psi_end -= 2 * M_PI;

      // Pad hits
      std::vector<double> hx, hy;
      if (ev->GetEntries() > 0) {
         for (auto &h : ((AtEvent *)ev->At(0))->GetHits()) {
            auto pos = h->GetPosition();
            hx.push_back(pos.X()); hy.push_back(pos.Y());
         }
      }

      // UKF
      std::vector<double> uX, uY;
      double ukfVx = 0, ukfVy = 0, phUKF = 0;
      bool haveUKF = false;
      if (te->GetEntries() > 0) {
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (best) {
            haveUKF = true;
            auto v = best->GetVertex();
            ukfVx = v.X(); ukfVy = v.Y();
            phUKF = best->GetKinematics().phi;
            for (auto &pp : best->GetSmoothedPositions()) {
               uX.push_back(pp.X()); uY.push_back(pp.Y());
            }
         }
      }

      double dPhi_mrad =
         haveUKF ? (phUKF - phMCfirst) * 1000. : 0;

      // Draw panel
      c->cd(panel + 1);
      gPad->SetGrid();
      auto *fr = new TH2F(Form("fr_%d", panel),
                          Form("evt %d:  #varphi_{0}=%.1f#circ  #Delta#varphi_{UKF-MC@first}=%.1f mrad;x (mm, beam);y (mm, transverse)",
                               eventIdx, phi_0 * 180. / M_PI, dPhi_mrad),
                          10, -450, 280, 10, -30, 100);
      fr->Draw();

      // Chamber + field-region boxes
      auto *cham = new TBox(0, 0, 256, 88);
      cham->SetFillStyle(0); cham->SetLineColor(kBlue + 1); cham->SetLineWidth(2);
      cham->Draw();
      auto *fld = new TBox(-500, -150, 300, 250);
      fld->SetFillStyle(0); fld->SetLineColor(kGray + 2); fld->SetLineStyle(2);
      fld->Draw();

      // MC helix outside chamber: analytical sample
      if (R > 0 && xFirst < 1e8) {
         auto *gPre = new TGraph();
         const int Npre = 80;
         for (int i = 0; i <= Npre; ++i) {
            double psi = psi_start + i * (psi_end - psi_start) / Npre;
            gPre->SetPoint(i, cx + R * std::cos(psi), cy + R * std::sin(psi));
         }
         gPre->SetLineColor(kGreen + 3); gPre->SetLineStyle(2); gPre->SetLineWidth(2);
         gPre->Draw("L");
      }

      // Pad hits
      if (!hx.empty()) {
         auto *gHits = new TGraph(hx.size(), hx.data(), hy.data());
         gHits->SetMarkerStyle(20); gHits->SetMarkerSize(0.4); gHits->SetMarkerColor(kBlack);
         gHits->Draw("P");
      }
      // MC inside chamber
      if (!mcX.empty()) {
         auto *gMC = new TGraph(mcX.size(), mcX.data(), mcY.data());
         gMC->SetLineColor(kGreen + 3); gMC->SetLineWidth(2);
         gMC->Draw("L");
      }
      // UKF
      if (!uX.empty()) {
         auto *gU = new TGraph(uX.size(), uX.data(), uY.data());
         gU->SetMarkerStyle(24); gU->SetMarkerSize(0.7); gU->SetMarkerColor(kRed + 1);
         gU->SetLineColor(kRed + 1); gU->SetLineWidth(2);
         gU->Draw("PL");
      }
      // Vertex markers
      auto *mMC = new TMarker(mcVx, mcVy, 29);
      mMC->SetMarkerColor(kGreen + 3); mMC->SetMarkerSize(1.5); mMC->Draw();
      if (haveUKF) {
         auto *mU = new TMarker(ukfVx, ukfVy, 29);
         mU->SetMarkerColor(kRed + 1); mU->SetMarkerSize(1.5); mU->Draw();
      }
   }

   // Shared legend
   c->cd(0);
   auto *leg = new TLegend(0.36, 0.97, 1.00, 1.00);
   leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
   leg->SetNColumns(4); leg->SetTextSize(0.02);
   auto *gMClg = new TGraph();  gMClg->SetLineColor(kGreen + 3);  gMClg->SetLineWidth(2);
   auto *gMClgD = new TGraph(); gMClgD->SetLineColor(kGreen + 3); gMClgD->SetLineWidth(2); gMClgD->SetLineStyle(2);
   auto *gUKFlg = new TGraph(); gUKFlg->SetMarkerStyle(24); gUKFlg->SetMarkerColor(kRed + 1);
   gUKFlg->SetLineColor(kRed + 1); gUKFlg->SetLineWidth(2);
   leg->AddEntry((TObject *)0, "TPC (blue) + field (gray)", "");
   leg->AddEntry(gMClgD, "MC helix (analytic, pre-chamber)", "l");
   leg->AddEntry(gMClg,  "MC (in chamber)", "l");
   leg->AddEntry(gUKFlg, "UKF smoothed", "pl");
   leg->Draw();

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
