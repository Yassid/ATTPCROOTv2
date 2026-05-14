/// @file make_distributions.C
/// @brief Per-scan-point residual distributions, one curve per p:
///   (a) Δp / p_MC
///   (b) Δθ  (mrad)
///   (c) Δφ  (mrad)
///   (d) Δy_vtx (mm) — transverse vertex residual
///
/// Reads ./data/output_ukf_HYDRA_p{P}.root and ./data/HYDRAsim_p{P}.root for
/// each P in {200, 400, 600, 800, 1000, 1200} (MeV/c).
///
/// Usage: root -b -q analysis/make_distributions.C

void make_distributions(const char *outPng = "data/dist_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.055, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.14);
   gStyle->SetPadBottomMargin(0.14);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.08);

   const std::vector<int> plist = {200, 400, 600, 800, 1000, 1200};
   const int nP = plist.size();
   const double mass_pi = 139.57039;
   const double radToMrad = 1000.;

   const int colors[] = {kBlue + 1, kRed + 1, kGreen + 3, kMagenta + 1, kOrange + 7, kCyan + 2};

   std::vector<TH1F *> hDP(nP), hDTh(nP), hDPhi(nP), hDVy(nP);

   for (int k = 0; k < nP; ++k) {
      int P = plist[k];
      TFile fS(Form("data/HYDRAsim_p%d.root", P));
      TFile fU(Form("data/output_ukf_HYDRA_p%d.root", P));
      auto *tS = (TTree *)fS.Get("cbmsim");
      auto *tU = (TTree *)fU.Get("cbmsim");
      if (!tS || !tU) continue;

      auto *trks = new TClonesArray("AtMCTrack");
      auto *pts  = new TClonesArray("AtMCPoint");
      auto *te = new TClonesArray("AtTrackingEvent");
      tS->SetBranchAddress("MCTrack", &trks);
      tS->SetBranchAddress("AtTpcPoint", &pts);
      tU->SetBranchAddress("AtTrackingEvent", &te);

      hDP[k]   = new TH1F(Form("hDP_%d", P),  "", 100, -0.3, 0.3);
      hDTh[k]  = new TH1F(Form("hDTh_%d", P), "", 100, -150., 150.);
      hDPhi[k] = new TH1F(Form("hDPhi_%d", P), "", 100, -300., 300.);
      hDVy[k]  = new TH1F(Form("hDVy_%d", P), "", 100, -60., 60.);
      for (auto *h : {hDP[k], hDTh[k], hDPhi[k], hDVy[k]}) {
         h->SetDirectory(nullptr);
         h->SetLineColor(colors[k % 6]);
         h->SetLineWidth(2);
      }

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
         double yV_mc = mc->GetStartY() * 10.;
         // Production-vertex angles — UKF back-rotates φ via helix and then
         // straight-tails to fBackExtrapTargetX upstream of the field region.
         double thMC = std::acos(Pz / pmcGeV);
         double phMC = std::atan2(Py, Px);
         (void)pts;

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
         auto kin = best->GetKinematics();
         double Efit = kin.kineticEnergy + mass_pi;
         double pfit = std::sqrt(Efit * Efit - mass_pi * mass_pi);
         double dPhi = kin.phi - phMC;
         while (dPhi > M_PI) dPhi -= 2 * M_PI;
         while (dPhi < -M_PI) dPhi += 2 * M_PI;
         hDP[k]->Fill((pfit - pmc) / pmc);
         hDTh[k]->Fill((kin.theta - thMC) * radToMrad);
         hDPhi[k]->Fill(dPhi * radToMrad);
         hDVy[k]->Fill(best->GetVertex().Y() - yV_mc);
      }
   }

   auto *c = new TCanvas("c", "HYDRA distributions", 1500, 1000);
   c->Divide(2, 2, 0.012, 0.015);

   auto drawSet = [&](int pad, std::vector<TH1F *> &h, const char *title) {
      c->cd(pad);
      gPad->SetGrid();
      double maxY = 0;
      for (auto *hp : h)
         if (hp && hp->Integral() > 0) {
            hp->Scale(1.0 / hp->Integral());
            maxY = std::max(maxY, hp->GetMaximum());
         }
      for (int k = 0; k < (int)h.size(); ++k) {
         if (!h[k]) continue;
         h[k]->SetTitle(title);
         h[k]->SetMaximum(1.25 * maxY);
         h[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
      }
      auto *leg = new TLegend(0.62, 0.55, 0.93, 0.92);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
      for (int k = 0; k < (int)h.size(); ++k)
         if (h[k]) leg->AddEntry(h[k], Form("%d MeV/c", plist[k]), "l");
      leg->Draw();
   };

   drawSet(1, hDP,   "(a) #Deltap / p_{MC};#Deltap / p_{MC};fraction");
   drawSet(2, hDTh,  "(b) #Delta#theta;#Delta#theta (mrad);fraction");
   drawSet(3, hDPhi, "(c) #Delta#varphi;#Delta#varphi (mrad);fraction");
   drawSet(4, hDVy,  "(d) #Deltay_{vtx} (transverse);#Deltay_{vtx} (mm);fraction");

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
