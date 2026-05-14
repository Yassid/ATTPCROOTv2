/// @file make_performance.C
/// @brief 2x3 performance summary for the HYDRA Prototype momentum scan:
///   (a) σ_p/p vs p (Gaussian core)
///   (b) σ_θ, σ_φ vs p (Gaussian core, mrad)
///   (c) bias_p vs p
///   (d) bias_θ, bias_φ vs p (mrad)
///   (e) vertex Δx-Δy scatter (all scan points, color-coded by p)
///   (f) σ_vtx_x, σ_vtx_y vs p; fit efficiency overlay on second y-axis
///
/// Reads ./data/output_ukf_HYDRA_p{P}.root and ./data/HYDRAsim_p{P}.root for
/// each P in {200, 400, 600, 800, 1000, 1200} (MeV/c) — matches scan_p.sh.
///
/// MC θ = arccos(Pz/p), φ = atan2(Py, Px) — evaluated at the FIRST AtMCPoint
/// inside drift_volume (the spatial reference where the UKF anchors its
/// fitted direction), not at the production vertex. Comparing to the vertex
/// would fold in the helix rotation between vertex and first cluster as an
/// apparent bias of ~30-150 mrad.
///
/// UKF: AtFittedTrack::Kinematics{theta,phi} (rad).
/// MC vertex (for residuals in (e), (f)): AtMCTrack::GetStart{X,Y,Z} (cm).
/// UKF vertex: AtFittedTrack::GetVertex() (mm).
///
/// Usage: root -b -q analysis/make_performance.C

void make_performance(const char *outPng = "data/perf_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.055, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.14);
   gStyle->SetPadBottomMargin(0.135);
   gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadTopMargin(0.08);

   const std::vector<int> plist = {200, 400, 600, 800, 1000, 1200};
   const int nP = plist.size();
   const double mass_pi = 139.57039; // MeV/c²
   const double radToDeg = 180. / M_PI;
   const double radToMrad = 1000.;

   std::vector<double> pV(nP), sigP(nP), sigPerr(nP), biasP(nP), effV(nP);
   std::vector<double> sigTh(nP), sigPhi(nP), biasTh(nP), biasPhi(nP);
   std::vector<double> sigVx(nP), sigVy(nP);

   // 2D vertex residual scatter for panel (e) — all events from all scan points.
   auto *hVxy = new TH2F("hVxy", "(e) vertex (UKF #minus MC) projection on pad plane;#Deltax (mm, beam);#Deltay (mm, transverse)",
                          120, -120., +500., 120, -60., +60.);
   hVxy->SetDirectory(nullptr);

   const int trkColors[] = {kBlue + 1, kRed + 1, kGreen + 3, kMagenta + 1, kOrange + 7, kCyan + 2};

   for (int k = 0; k < nP; ++k) {
      int P = plist[k];
      pV[k] = P;
      TString simF = Form("data/HYDRAsim_p%d.root", P);
      TString ukfF = Form("data/output_ukf_HYDRA_p%d.root", P);
      TFile fS(simF), fU(ukfF);
      auto *tS = (TTree *)fS.Get("cbmsim");
      auto *tU = (TTree *)fU.Get("cbmsim");
      if (!tS || !tU) { std::cerr << "missing " << simF << " or " << ukfF << "\n"; continue; }

      auto *trks = new TClonesArray("AtMCTrack");
      auto *pts  = new TClonesArray("AtMCPoint");
      auto *te = new TClonesArray("AtTrackingEvent");
      tS->SetBranchAddress("MCTrack", &trks);
      tS->SetBranchAddress("AtTpcPoint", &pts);
      tU->SetBranchAddress("AtTrackingEvent", &te);

      TH1F hDP("hDP", "", 100, -0.3, 0.3);
      TH1F hDTh("hDTh", "", 120, -120., 120.); // mrad
      TH1F hDPhi("hDPhi", "", 120, -120., 120.);
      TH1F hVx("hVx", "", 100, -100., +500.);
      TH1F hVy("hVy", "", 100, -60., +60.);
      hDP.SetDirectory(nullptr); hDTh.SetDirectory(nullptr); hDPhi.SetDirectory(nullptr);
      hVx.SetDirectory(nullptr); hVy.SetDirectory(nullptr);

      int nThr = 0, nFit = 0;
      Long64_t n = std::min(tS->GetEntries(), tU->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         tS->GetEntry(i);
         tU->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;
         double Px = mc->GetPx(), Py = mc->GetPy(), Pz = mc->GetPz();
         double pmcGeV = std::sqrt(Px * Px + Py * Py + Pz * Pz);
         double pmc = pmcGeV * 1000.;
         double xV_mc = mc->GetStartX() * 10., yV_mc = mc->GetStartY() * 10.;
         // Angles at the production vertex. UKF now back-rotates φ through
         // the helix AND straight-tails to fBackExtrapTargetX in field-free
         // region, so this is the apples-to-apples reference.
         double thMC = std::acos(Pz / pmcGeV);
         double phMC = std::atan2(Py, Px);
         (void)pts;
         ++nThr;
         if (te->GetEntries() == 0) continue;
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double c = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (c < bestChi) { bestChi = c; best = t.get(); }
         }
         if (!best) continue;
         ++nFit;
         auto kin = best->GetKinematics();
         double KEfit = kin.kineticEnergy;
         double Efit = KEfit + mass_pi;
         double pfit = std::sqrt(Efit * Efit - mass_pi * mass_pi);
         hDP.Fill((pfit - pmc) / pmc);
         double dTh = (kin.theta - thMC) * radToMrad;
         double dPhi = (kin.phi - phMC);
         // Wrap dPhi to [-π, π] before scaling to mrad
         while (dPhi > M_PI) dPhi -= 2 * M_PI;
         while (dPhi < -M_PI) dPhi += 2 * M_PI;
         hDTh.Fill(dTh);
         hDPhi.Fill(dPhi * radToMrad);

         auto vrec = best->GetVertex();
         hVx.Fill(vrec.X() - xV_mc);
         hVy.Fill(vrec.Y() - yV_mc);
         hVxy->Fill(vrec.X() - xV_mc, vrec.Y() - yV_mc);
      }

      auto fitCore = [](TH1F &h, double &mu, double &sig, double &sigE) {
         double rms = h.GetRMS();
         TFitResultPtr fr = h.Fit("gaus", "SQR0", "", -3 * rms, 3 * rms);
         mu = fr.Get() ? fr->Parameter(1) : h.GetMean();
         sig = fr.Get() ? fr->Parameter(2) : rms;
         sigE = fr.Get() ? fr->ParError(2) : sig / std::sqrt(2.0 * std::max(1, (int)h.GetEntries()));
      };
      double mu, sig, sigE;
      fitCore(hDP, mu, sig, sigE);  biasP[k] = mu * 100.; sigP[k] = sig * 100.; sigPerr[k] = sigE * 100.;
      fitCore(hDTh, mu, sig, sigE); biasTh[k] = mu; sigTh[k] = sig;
      fitCore(hDPhi, mu, sig, sigE); biasPhi[k] = mu; sigPhi[k] = sig;
      fitCore(hVx, mu, sig, sigE); sigVx[k] = sig;
      fitCore(hVy, mu, sig, sigE); sigVy[k] = sig;
      effV[k] = nThr > 0 ? 100. * nFit / nThr : 0.;

      std::cout << "p=" << P << "  Nfit=" << nFit << "/" << nThr
                << "  σ_p=" << sigP[k] << "% bias_p=" << biasP[k] << "%"
                << "  σ_θ=" << sigTh[k] << " mrad  σ_φ=" << sigPhi[k] << " mrad"
                << "  σ_vtx (x,y)=(" << sigVx[k] << ", " << sigVy[k] << ") mm"
                << "  eff=" << effV[k] << "%\n";
   }

   auto *c = new TCanvas("c", "HYDRA performance", 1800, 1000);
   c->Divide(3, 2, 0.012, 0.015);

   // (a) σ_p/p
   c->cd(1); gPad->SetGrid();
   auto *gSig = new TGraphErrors(nP, pV.data(), sigP.data(), nullptr, sigPerr.data());
   gSig->SetTitle("(a) #sigma_{p}/p vs p (Gaussian core);p_{MC} (MeV/c);#sigma_{p}/p (%)");
   gSig->SetMarkerStyle(20); gSig->SetMarkerSize(1.5);
   gSig->SetMarkerColor(kBlue + 1); gSig->SetLineColor(kBlue + 1); gSig->SetLineWidth(2);
   gSig->Draw("APL");
   gSig->GetYaxis()->SetRangeUser(0., 10.);
   gSig->GetXaxis()->SetLimits(0., 1400.);

   // (b) σ_θ, σ_φ
   c->cd(2); gPad->SetGrid();
   auto *gTh = new TGraph(nP, pV.data(), sigTh.data());
   auto *gPhi = new TGraph(nP, pV.data(), sigPhi.data());
   for (auto *g : {gTh, gPhi}) { g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetLineWidth(2); }
   gTh->SetMarkerColor(kBlue + 1); gTh->SetLineColor(kBlue + 1);
   gPhi->SetMarkerColor(kRed + 1);  gPhi->SetLineColor(kRed + 1);
   double amax = std::max(*std::max_element(sigTh.begin(), sigTh.end()),
                           *std::max_element(sigPhi.begin(), sigPhi.end()));
   gTh->SetTitle("(b) #sigma_{#theta}, #sigma_{#varphi} vs p (Gaussian core);p_{MC} (MeV/c);#sigma (mrad)");
   gTh->Draw("APL");
   gTh->GetYaxis()->SetRangeUser(0., 1.3 * amax);
   gTh->GetXaxis()->SetLimits(0., 1400.);
   gPhi->Draw("PL same");
   auto *legB = new TLegend(0.60, 0.70, 0.93, 0.90);
   legB->SetBorderSize(0); legB->SetFillStyle(0); legB->SetTextSize(0.05);
   legB->AddEntry(gTh, "#sigma_{#theta}", "lp");
   legB->AddEntry(gPhi, "#sigma_{#varphi}", "lp");
   legB->Draw();

   // (c) bias_p
   c->cd(3); gPad->SetGrid();
   auto *gBP = new TGraph(nP, pV.data(), biasP.data());
   gBP->SetTitle("(c) bias_p vs p;p_{MC} (MeV/c);<#Deltap/p_{MC}> (%)");
   gBP->SetMarkerStyle(20); gBP->SetMarkerSize(1.5);
   gBP->SetMarkerColor(kRed + 1); gBP->SetLineColor(kRed + 1); gBP->SetLineWidth(2);
   gBP->Draw("APL");
   gBP->GetYaxis()->SetRangeUser(-10., 2.);
   gBP->GetXaxis()->SetLimits(0., 1400.);
   auto *zL = new TLine(0., 0., 1400., 0.);
   zL->SetLineStyle(2); zL->SetLineColor(kGray + 2);
   zL->Draw();

   // (d) bias_θ, bias_φ
   c->cd(4); gPad->SetGrid();
   auto *gBTh = new TGraph(nP, pV.data(), biasTh.data());
   auto *gBPhi = new TGraph(nP, pV.data(), biasPhi.data());
   for (auto *g : {gBTh, gBPhi}) { g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetLineWidth(2); }
   gBTh->SetMarkerColor(kBlue + 1); gBTh->SetLineColor(kBlue + 1);
   gBPhi->SetMarkerColor(kRed + 1);  gBPhi->SetLineColor(kRed + 1);
   double bmin = 0, bmax = 0;
   for (int k = 0; k < nP; ++k) {
      bmin = std::min({bmin, biasTh[k], biasPhi[k]});
      bmax = std::max({bmax, biasTh[k], biasPhi[k]});
   }
   gBTh->SetTitle("(d) bias_{#theta}, bias_{#varphi} vs p;p_{MC} (MeV/c);<#Delta> (mrad)");
   gBTh->Draw("APL");
   gBTh->GetYaxis()->SetRangeUser(1.3 * bmin - 5, 1.3 * bmax + 5);
   gBTh->GetXaxis()->SetLimits(0., 1400.);
   gBPhi->Draw("PL same");
   auto *zL2 = new TLine(0., 0., 1400., 0.);
   zL2->SetLineStyle(2); zL2->SetLineColor(kGray + 2);
   zL2->Draw();
   auto *legD = new TLegend(0.60, 0.70, 0.93, 0.90);
   legD->SetBorderSize(0); legD->SetFillStyle(0); legD->SetTextSize(0.05);
   legD->AddEntry(gBTh, "bias_{#theta}", "lp");
   legD->AddEntry(gBPhi, "bias_{#varphi}", "lp");
   legD->Draw();

   // (e) vertex Δx-Δy 2D scatter
   c->cd(5); gPad->SetGrid();
   hVxy->Draw("COLZ");
   // Mark MC vertex residual = 0
   auto *m0 = new TMarker(0, 0, 29);
   m0->SetMarkerColor(kRed + 1); m0->SetMarkerSize(2.0);
   m0->Draw();
   auto *txt = new TLatex(0, 0, "  MC");
   txt->SetTextSize(0.04); txt->SetTextColor(kRed + 1); txt->Draw();

   // (f) σ_vtx vs p + efficiency
   c->cd(6); gPad->SetGrid();
   auto *gVx = new TGraph(nP, pV.data(), sigVx.data());
   auto *gVy = new TGraph(nP, pV.data(), sigVy.data());
   for (auto *g : {gVx, gVy}) { g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetLineWidth(2); }
   gVx->SetMarkerColor(kBlue + 1); gVx->SetLineColor(kBlue + 1);
   gVy->SetMarkerColor(kGreen + 3); gVy->SetLineColor(kGreen + 3);
   double vmax = std::max(*std::max_element(sigVx.begin(), sigVx.end()),
                           *std::max_element(sigVy.begin(), sigVy.end()));
   gVx->SetTitle("(f) #sigma_{vtx,x}, #sigma_{vtx,y} vs p;p_{MC} (MeV/c);#sigma_{vtx} (mm)");
   gVx->Draw("APL");
   gVx->GetYaxis()->SetRangeUser(0., 1.4 * vmax);
   gVx->GetXaxis()->SetLimits(0., 1400.);
   gVy->Draw("PL same");
   auto *legF = new TLegend(0.60, 0.55, 0.93, 0.78);
   legF->SetBorderSize(0); legF->SetFillStyle(0); legF->SetTextSize(0.045);
   legF->AddEntry(gVx, "#sigma_{x} (beam)", "lp");
   legF->AddEntry(gVy, "#sigma_{y} (transverse)", "lp");
   legF->Draw();
   // Efficiency annotation
   auto *txtEff = new TLatex();
   txtEff->SetTextSize(0.05);
   txtEff->SetTextColor(kBlack);
   txtEff->DrawLatexNDC(0.18, 0.86,
                        Form("fit eff: %.1f-%.1f%%",
                             *std::min_element(effV.begin(), effV.end()),
                             *std::max_element(effV.begin(), effV.end())));

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
