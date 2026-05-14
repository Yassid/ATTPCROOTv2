/// @file display_event3d.C
/// @brief Single-event explanatory figure: xy projection (top-down) and
/// xz projection (side, drift direction). Overlays MC trajectory, pad
/// hits, UKF smoothed track, MC vertex and UKF back-extrap vertex.
/// Annotates the relevant angles to help diagnose why the back-extrap
/// reconstruction is short of the production vertex.
///
/// Default inputs are the most recent fixed-φ narrow-field scan files.
/// Pass eventIdx to choose another event.
///
/// Usage:
///   root -b -q 'analysis/display_event3d.C(800, 100)'

void display_event3d(int P_MeV = 800, int eventIdx = 100,
                     const char *runDir = "data/runs/scan_phi0_narrowfield_2k",
                     const char *outPng = "data/event3d_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.05, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.04, "XY");
   gStyle->SetPadLeftMargin(0.12);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.08);

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

   tS->GetEntry(eventIdx);
   tD->GetEntry(eventIdx);
   tU->GetEntry(eventIdx);

   auto *mc = (AtMCTrack *)trks->At(0);
   double mcVx = mc->GetStartX() * 10., mcVy = mc->GetStartY() * 10., mcVz = mc->GetStartZ() * 10.;
   double mcPx = mc->GetPx(), mcPy = mc->GetPy(), mcPz = mc->GetPz();
   double phMC_vtx = std::atan2(mcPy, mcPx) * 180. / M_PI; // deg

   // Collect MC points in drift_volume (primary pion)
   std::vector<double> mcX, mcY, mcZ;
   double phMC_first = 0;
   double xFirst = 1e9, yFirst = 0;
   for (int j = 0; j < pts->GetEntries(); ++j) {
      auto *p = (AtMCPoint *)pts->At(j);
      if (p->GetVolName() != TString("drift_volume")) continue;
      if (p->GetTrackID() != 0) continue;
      double x_mm = p->GetX() * 10., y_mm = p->GetY() * 10., z_mm = p->GetZ() * 10.;
      mcX.push_back(x_mm);
      mcY.push_back(y_mm);
      mcZ.push_back(z_mm);
      if (x_mm < xFirst) {
         xFirst = x_mm;
         yFirst = y_mm;
         phMC_first = std::atan2(p->GetPy(), p->GetPx()) * 180. / M_PI;
      }
   }

   // Pad hits from AtEvent
   std::vector<double> hx, hy, hz, hQ;
   if (ev->GetEntries() > 0) {
      for (auto &h : ((AtEvent *)ev->At(0))->GetHits()) {
         auto pos = h->GetPosition();
         hx.push_back(pos.X());
         hy.push_back(pos.Y());
         hz.push_back(pos.Z());
         hQ.push_back(h->GetCharge());
      }
   }

   // UKF: best track + smoothed positions + back-extrap vertex
   std::vector<double> uX, uY, uZ;
   double ukfVx = 0, ukfVy = 0, ukfVz = 0;
   double phUKF = 0;
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
         ukfVx = v.X(); ukfVy = v.Y(); ukfVz = v.Z();
         phUKF = best->GetKinematics().phi * 180. / M_PI;
         for (auto &pp : best->GetSmoothedPositions()) {
            uX.push_back(pp.X());
            uY.push_back(pp.Y());
            uZ.push_back(pp.Z());
         }
      }
   }

   // ------------------------------------------------------------------
   // Two panels: top-down xy (left), side xz (right).
   auto *c = new TCanvas("c", Form("Event %d, p = %d MeV/c", eventIdx, P_MeV), 1600, 800);
   c->Divide(2, 1, 0.01, 0.02);

   // ----- xy projection -----
   c->cd(1); gPad->SetGrid();
   auto *fr1 = new TH2F("fr1",
                        Form("evt %d, p = %d MeV/c — top-down (xy plane);x (mm, beam);y (mm, transverse)",
                             eventIdx, P_MeV),
                        10, -450, 280, 10, -30, 100);
   fr1->Draw();

   // Chamber boundary
   auto *box = new TBox(0, 0, 256, 88);
   box->SetFillStyle(0);
   box->SetLineColor(kBlue + 1); box->SetLineWidth(2);
   box->Draw();
   // Field-region boundary
   auto *fieldBox = new TBox(-20, -20, 300, 120);
   fieldBox->SetFillStyle(0);
   fieldBox->SetLineColor(kGray + 2); fieldBox->SetLineStyle(2);
   fieldBox->Draw();
   auto *labelField = new TLatex(-20, 95, " field region");
   labelField->SetTextSize(0.028); labelField->SetTextColor(kGray + 2);
   labelField->Draw();
   auto *labelChamber = new TLatex(256, -10, "TPC ");
   labelChamber->SetTextSize(0.030); labelChamber->SetTextColor(kBlue + 1);
   labelChamber->SetTextAlign(33); labelChamber->Draw();

   // Pad hits
   if (!hx.empty()) {
      auto *gHits = new TGraph(hx.size(), hx.data(), hy.data());
      gHits->SetMarkerStyle(20); gHits->SetMarkerSize(0.4); gHits->SetMarkerColor(kBlack);
      gHits->Draw("P");
   }

   // MC trajectory
   if (!mcX.empty()) {
      auto *gMC = new TGraph(mcX.size(), mcX.data(), mcY.data());
      gMC->SetLineColor(kGreen + 3); gMC->SetLineWidth(2);
      gMC->Draw("L");
      // Straight line from MC vertex to first drift-volume MC point (field-free portion)
      double tlX[2] = {mcVx, mcX.front()};
      double tlY[2] = {mcVy, mcY.front()};
      auto *gPre = new TGraph(2, tlX, tlY);
      gPre->SetLineColor(kGreen + 3); gPre->SetLineStyle(2); gPre->SetLineWidth(2);
      gPre->Draw("L");
   }

   // UKF smoothed positions
   if (!uX.empty()) {
      auto *gUKF = new TGraph(uX.size(), uX.data(), uY.data());
      gUKF->SetMarkerStyle(24); gUKF->SetMarkerSize(0.9); gUKF->SetMarkerColor(kRed + 1);
      gUKF->SetLineColor(kRed + 1); gUKF->SetLineWidth(2);
      gUKF->Draw("PL");
   }

   // MC vertex marker
   auto *mMC = new TMarker(mcVx, mcVy, 29);
   mMC->SetMarkerColor(kGreen + 3); mMC->SetMarkerSize(2.0); mMC->Draw();
   auto *labMC = new TLatex(mcVx + 5, mcVy - 8, "MC vertex");
   labMC->SetTextColor(kGreen + 3); labMC->SetTextSize(0.030); labMC->Draw();

   // UKF back-extrap vertex marker
   if (haveUKF) {
      auto *mU = new TMarker(ukfVx, ukfVy, 29);
      mU->SetMarkerColor(kRed + 1); mU->SetMarkerSize(2.0); mU->Draw();
      auto *labU = new TLatex(ukfVx + 5, ukfVy + 8, "UKF POCA");
      labU->SetTextColor(kRed + 1); labU->SetTextSize(0.030); labU->Draw();
   }

   // Angle arrows: short tangent lines at each angle reference point.
   auto drawArrow = [&](double x0, double y0, double phiDeg, double len, int color) {
      double phiRad = phiDeg * M_PI / 180.;
      auto *a = new TArrow(x0, y0, x0 + len * std::cos(phiRad), y0 + len * std::sin(phiRad), 0.02, "|>");
      a->SetLineColor(color); a->SetFillColor(color); a->SetLineWidth(2);
      a->Draw();
   };
   drawArrow(mcVx, mcVy, phMC_vtx, 80, kGreen + 3);  // MC at production
   if (!mcX.empty())
      drawArrow(xFirst, yFirst, phMC_first, 60, kGreen + 1); // MC at first MC point
   if (haveUKF)
      drawArrow(ukfVx, ukfVy, phUKF, 60, kRed + 1);  // UKF at POCA

   // Annotation box
   auto *box2 = new TPaveText(0.45, 0.78, 0.93, 0.92, "NDC");
   box2->SetFillColor(0); box2->SetBorderSize(1); box2->SetTextAlign(12); box2->SetTextSize(0.025);
   box2->AddText(Form("#varphi^{MC}_{prod vtx}        = %.2f#circ", phMC_vtx));
   box2->AddText(Form("#varphi^{MC}_{first MC pt}    = %.2f#circ  (#Delta = %.1f mrad)",
                       phMC_first, (phMC_first - phMC_vtx) * M_PI / 180. * 1000.));
   if (haveUKF) {
      box2->AddText(Form("#varphi^{UKF}_{POCA}           = %.2f#circ  (#Delta vs MC@first = %.1f mrad)",
                          phUKF, (phUKF - phMC_first) * M_PI / 180. * 1000.));
      box2->AddText(Form("UKF POCA       (%.1f, %.1f) mm", ukfVx, ukfVy));
      box2->AddText(Form("first MC point  (%.1f, %.1f) mm", xFirst, yFirst));
   }
   box2->Draw();

   // Legend
   auto *leg = new TLegend(0.16, 0.62, 0.45, 0.88);
   leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.028);
   leg->AddEntry((TObject *)0, "TPC active   (blue box)", "");
   leg->AddEntry((TObject *)0, "field region (gray box)", "");
   auto *gMClg = new TGraph(); gMClg->SetLineColor(kGreen + 3); gMClg->SetLineWidth(2);
   leg->AddEntry(gMClg, "MC track (in chamber)", "l");
   auto *gPrelg = new TGraph(); gPrelg->SetLineColor(kGreen + 3); gPrelg->SetLineStyle(2); gPrelg->SetLineWidth(2);
   leg->AddEntry(gPrelg, "MC pre-chamber (field-free)", "l");
   auto *gUKFlg = new TGraph(); gUKFlg->SetMarkerStyle(24); gUKFlg->SetMarkerColor(kRed + 1);
   gUKFlg->SetLineColor(kRed + 1); gUKFlg->SetLineWidth(2);
   leg->AddEntry(gUKFlg, "UKF smoothed track", "pl");
   leg->Draw();

   // ----- xz projection -----
   c->cd(2); gPad->SetGrid();
   auto *fr2 = new TH2F("fr2",
                        Form("evt %d, p = %d MeV/c — side (xz);x (mm, beam);z (mm, drift)",
                             eventIdx, P_MeV),
                        10, -450, 280, 10, -20, 320);
   fr2->Draw();
   // Chamber boundary
   auto *box3 = new TBox(0, 0, 256, 294);
   box3->SetFillStyle(0); box3->SetLineColor(kBlue + 1); box3->SetLineWidth(2);
   box3->Draw();

   if (!hx.empty()) {
      auto *gHitsXZ = new TGraph(hx.size(), hx.data(), hz.data());
      gHitsXZ->SetMarkerStyle(20); gHitsXZ->SetMarkerSize(0.4); gHitsXZ->SetMarkerColor(kBlack);
      gHitsXZ->Draw("P");
   }
   if (!mcX.empty()) {
      auto *gMCxz = new TGraph(mcX.size(), mcX.data(), mcZ.data());
      gMCxz->SetLineColor(kGreen + 3); gMCxz->SetLineWidth(2);
      gMCxz->Draw("L");
      double tlX[2] = {mcVx, mcX.front()};
      double tlZ[2] = {mcVz, mcZ.front()};
      auto *gPreXZ = new TGraph(2, tlX, tlZ);
      gPreXZ->SetLineColor(kGreen + 3); gPreXZ->SetLineStyle(2); gPreXZ->SetLineWidth(2);
      gPreXZ->Draw("L");
   }
   if (!uX.empty()) {
      auto *gUxz = new TGraph(uX.size(), uX.data(), uZ.data());
      gUxz->SetMarkerStyle(24); gUxz->SetMarkerSize(0.9); gUxz->SetMarkerColor(kRed + 1);
      gUxz->SetLineColor(kRed + 1); gUxz->SetLineWidth(2);
      gUxz->Draw("PL");
   }
   auto *mMCz = new TMarker(mcVx, mcVz, 29);
   mMCz->SetMarkerColor(kGreen + 3); mMCz->SetMarkerSize(2.0); mMCz->Draw();
   if (haveUKF) {
      auto *mUz = new TMarker(ukfVx, ukfVz, 29);
      mUz->SetMarkerColor(kRed + 1); mUz->SetMarkerSize(2.0); mUz->Draw();
   }

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
   std::cout << "MC vertex: (" << mcVx << ", " << mcVy << ", " << mcVz << ") mm  φ_MC = " << phMC_vtx << "°\n";
   std::cout << "First MC pt in chamber: (" << xFirst << ", " << yFirst << ") mm  φ_MC@first = " << phMC_first << "°\n";
   if (haveUKF) {
      std::cout << "UKF POCA: (" << ukfVx << ", " << ukfVy << ", " << ukfVz << ") mm  φ_UKF = " << phUKF << "°\n";
   }
}
