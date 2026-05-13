/// @file display_fits.C
/// @brief Pad-plane projections of 6 example events with the UKF-fitted
/// trajectory overlaid: pad hits (color), MC truth (green line), UKF
/// smoothed positions (red markers + line).
///
/// Reads ./data/output_digi.root, ./data/output_ukf_HYDRA.root and
/// ./data/HYDRAsim.root (the most-recent 1-pion 2 T sandbox run, no suffix).
/// Selects the 6 events with the longest primary-pion trajectory inside
/// drift_volume, matching display_tracks.C.
///
/// Usage: root -b -q analysis/display_fits.C

void display_fits(const char *digiFile = "data/output_digi.root",
                  const char *ukfFile  = "data/output_ukf_HYDRA.root",
                  const char *simFile  = "data/HYDRAsim.root",
                  const char *outPrefix = "data/display_fits_HYDRA")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.08);

   TFile fSim(simFile), fDigi(digiFile), fUKF(ukfFile);
   auto *tSim  = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   auto *tUKF  = (TTree *)fUKF.Get("cbmsim");

   auto *trks  = new TClonesArray("AtMCTrack");
   auto *pts   = new TClonesArray("AtMCPoint");
   auto *evArr = new TClonesArray("AtEvent");
   auto *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tSim->SetBranchAddress("AtTpcPoint", &pts);
   tDigi->SetBranchAddress("AtEventH", &evArr);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   // Transposed-display pad-plane: x (beam) vertical, y (transverse) horizontal.
   const double padXmin = 0., padXmax = 256.; // beam, mm
   const double padYmin = 0., padYmax = 88.;  // transverse, mm

   // Pick the 6 events with the longest primary-pion drift-volume trajectory
   struct Pick { Long64_t entry = -1; int nPiPts = 0; int nh = 0; };
   std::array<Pick, 6> picks;
   Long64_t n = std::min({tSim->GetEntries(), tDigi->GetEntries(), tUKF->GetEntries()});
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      if (evArr->GetEntries() == 0) continue;
      std::set<int> piIds;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (std::abs(mc->GetPdgCode()) == 211 && mc->GetMotherId() == -1) piIds.insert(j);
      }
      int nPiPts = 0;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (piIds.count(p->GetTrackID())) ++nPiPts;
      }
      int nh = ((AtEvent *)evArr->At(0))->GetNumHits();
      for (size_t b = 0; b < picks.size(); ++b) {
         if (nPiPts > picks[b].nPiPts) {
            for (size_t j = picks.size() - 1; j > b; --j) picks[j] = picks[j - 1];
            picks[b] = {i, nPiPts, nh};
            break;
         }
      }
   }

   auto *c = new TCanvas("c", "HYDRA fits", 1700, 1100);
   c->Divide(3, 2, 0.01, 0.02);

   for (size_t b = 0; b < picks.size(); ++b) {
      c->cd(b + 1);
      gPad->SetGrid();
      if (picks[b].entry < 0) continue;

      tSim->GetEntry(picks[b].entry);
      tDigi->GetEntry(picks[b].entry);
      tUKF->GetEntry(picks[b].entry);

      // Build pad-hit histogram (transposed: hist x ← phys y, hist y ← phys x)
      auto *hPad = new TH2F(Form("hPad_%zu", b),
                            Form("evt %lld: %d hits;y (mm);x (mm, beam)",
                                 picks[b].entry, picks[b].nh),
                            44, padYmin, padYmax, 128, padXmin, padXmax);
      hPad->SetDirectory(nullptr);
      auto *ev = (AtEvent *)evArr->At(0);
      for (auto &h : ev->GetHits()) {
         auto pos = h->GetPosition();
         hPad->Fill(pos.Y(), pos.X(), h->GetCharge());
      }
      hPad->Draw("colz");

      // MC primary pion trajectory (green)
      std::set<int> piIds;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (std::abs(mc->GetPdgCode()) == 211 && mc->GetMotherId() == -1) piIds.insert(j);
      }
      auto *gMC = new TGraph();
      gMC->SetLineColor(kGreen + 3); gMC->SetLineWidth(2);
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (!piIds.count(p->GetTrackID())) continue;
         // Transposed: graph x ← phys y, graph y ← phys x  (phys in cm → mm)
         gMC->SetPoint(gMC->GetN(), p->GetY() * 10., p->GetX() * 10.);
      }
      gMC->Draw("L same");

      // UKF smoothed positions (red markers + line)
      if (teArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)teArr->At(0);
         auto &fitted = te->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (best) {
            auto &sp = best->GetSmoothedPositions();
            auto *gUKF = new TGraph();
            for (auto &pp : sp) gUKF->SetPoint(gUKF->GetN(), pp.Y(), pp.X());
            gUKF->SetMarkerStyle(24); gUKF->SetMarkerSize(0.9);
            gUKF->SetMarkerColor(kRed + 1);
            gUKF->SetLineColor(kRed + 1); gUKF->SetLineWidth(2);
            gUKF->Draw("PL same");
         }
      }

      // Pad-plane boundary
      auto *box = new TLine();
      box->SetLineColor(kBlue + 1); box->SetLineWidth(2);
      box->DrawLine(padYmin, padXmin, padYmax, padXmin);
      box->DrawLine(padYmax, padXmin, padYmax, padXmax);
      box->DrawLine(padYmax, padXmax, padYmin, padXmax);
      box->DrawLine(padYmin, padXmax, padYmin, padXmin);
      // MC vertex marker
      const double yMid = 0.5 * (padYmin + padYmax);
      auto *vtxArrow = new TArrow(yMid, padXmin - 15, yMid, padXmin + 5, 0.02, "|>");
      vtxArrow->SetLineColor(kGreen + 3); vtxArrow->SetFillColor(kGreen + 3);
      vtxArrow->SetLineWidth(2);
      vtxArrow->Draw();
   }

   // Single legend on the canvas
   auto *leg = new TLegend(0.40, 0.95, 0.99, 0.99);
   leg->SetBorderSize(0); leg->SetFillStyle(1001); leg->SetFillColor(kWhite);
   leg->SetNColumns(3); leg->SetTextSize(0.025);
   auto *gMCleg = new TGraph(); gMCleg->SetLineColor(kGreen + 3); gMCleg->SetLineWidth(2);
   auto *gUKFleg = new TGraph(); gUKFleg->SetMarkerStyle(24); gUKFleg->SetMarkerColor(kRed + 1);
   gUKFleg->SetLineColor(kRed + 1); gUKFleg->SetLineWidth(2);
   leg->AddEntry((TObject*)0, "pad hits (colz)", "");
   leg->AddEntry(gMCleg, "MC primary", "l");
   leg->AddEntry(gUKFleg, "UKF smoothed", "pl");
   c->cd(0);
   leg->Draw();

   TString png = TString(outPrefix) + ".png";
   TString pdf = TString(outPrefix) + ".pdf";
   c->SaveAs(png); c->SaveAs(pdf);
   std::cout << "Wrote " << outPrefix << ".{png,pdf}\n";
}
