/// @file display_tracks.C
/// @brief Display a grid of example tracks projected onto the HYDRA pad
/// plane (x, y). Reads digi+sim output, picks events spanning a range of
/// θ values, draws hits + the pad-plane boundary.
///
/// Usage:
///   root -b -q 'analysis/display_tracks.C("data/output_digi.root",
///       "data/HYDRAsim.root","data/display_tracks_HYDRA")'

void display_tracks(const char *digiFile = "data/output_digi.root",
                    const char *simFile = "data/HYDRAsim.root",
                    const char *outPrefix = "data/display_tracks_HYDRA")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim(simFile);
   TFile fDigi(digiFile);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   auto *trks = new TClonesArray("AtMCTrack");
   auto *pts = new TClonesArray("AtMCPoint");
   auto *evArr = new TClonesArray("AtEvent");
   auto *patArr = new TClonesArray("AtPatternEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tSim->SetBranchAddress("AtTpcPoint", &pts);
   tDigi->SetBranchAddress("AtEventH", &evArr);
   tDigi->SetBranchAddress("AtPatternEvent", &patArr);

   // HYDRA Prototype MCPP: 256 × 88 mm² active area, 2 mm pads, lower-
   // left vertex anchored at (0, 0). X = beam, Y = transverse.
   const double padXmin = 0.0,   padXmax = 256.0; // 256 mm long (beam)
   const double padYmin = 0.0,   padYmax = 88.0;  // 88 mm transverse
   const double padPitch = 2.0;

   // Pick the 6 events with the LONGEST primary-pion trajectory inside
   // drift_volume — these are clean non-interacting passages. Length is
   // approximated by the count of AtMCPoints whose track id matches a
   // primary π±.
   struct Pick { Long64_t entry = -1; int nPiPts = 0; int nh = 0; };
   std::array<Pick, 6> picks;
   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      if (evArr->GetEntries() == 0) continue;

      std::set<int> piIds;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (std::abs(mc->GetPdgCode()) == 211 && mc->GetMotherId() == -1)
            piIds.insert(j);
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

   auto *c = new TCanvas("c", "HYDRA pad-plane projections", 1700, 1100);
   c->Divide(3, 2, 0.01, 0.02);

   for (size_t b = 0; b < picks.size(); ++b) {
      c->cd(b + 1);
      gPad->SetGrid();

      if (picks[b].entry < 0) {
         auto *txt = new TLatex(0.5, 0.5, "(no event)");
         txt->SetNDC(); txt->SetTextAlign(22);
         txt->Draw();
         continue;
      }

      tSim->GetEntry(picks[b].entry);
      tDigi->GetEntry(picks[b].entry);

      // Count primary pions in this event
      int nPi = 0;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (std::abs(mc->GetPdgCode()) == 211 && mc->GetMotherId() == -1) ++nPi;
      }

      // Display: physical-y (transverse, 88 mm) on horizontal,
      // physical-x (beam, 256 mm) on vertical, so beam goes up.
      auto *hPad = new TH2F(Form("hPad_%zu", b),
                            Form("evt %lld:  %d pions, %d pad hits;y (mm);x (mm, beam)",
                                 picks[b].entry, nPi, picks[b].nh),
                            44, padYmin, padYmax,
                            128, padXmin, padXmax);
      hPad->SetMinimum(0);

      auto *ev = (AtEvent *)evArr->At(0);
      auto &hits = ev->GetHits();
      TGraph *gAll = new TGraph();
      for (size_t i = 0; i < hits.size(); ++i) {
         auto pos = hits[i]->GetPosition(); // mm in digi frame
         // Transposed display: hist x ← phys y, hist y ← phys x.
         hPad->Fill(pos.Y(), pos.X(), hits[i]->GetCharge());
         gAll->SetPoint(i, pos.Y(), pos.X());
      }
      hPad->Draw("colz");
      gAll->SetMarkerStyle(20);
      gAll->SetMarkerColor(kBlack);
      gAll->SetMarkerSize(0.5);
      gAll->Draw("P same");

      // MC TRUTH trajectories — primary pions only (filter out δ-rays
      // and other secondaries). AtMCPoint coords are in cm; pad plane
      // is in mm.
      std::set<int> primaryPiIds;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (std::abs(mc->GetPdgCode()) == 211 && mc->GetMotherId() == -1)
            primaryPiIds.insert(j);
      }
      std::map<int, TGraph *> gMCmap;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         int tid = p->GetTrackID();
         if (primaryPiIds.count(tid) == 0) continue;
         auto it = gMCmap.find(tid);
         if (it == gMCmap.end()) it = gMCmap.emplace(tid, new TGraph()).first;
         // Transposed display: graph x ← phys y, graph y ← phys x.
         it->second->SetPoint(it->second->GetN(), p->GetY() * 10., p->GetX() * 10.);
      }
      const int trkColors[] = {kGreen + 3, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet + 1};
      int colorIdx = 0;
      for (auto &kv : gMCmap) {
         kv.second->SetLineColor(trkColors[colorIdx++ % 5]);
         kv.second->SetLineWidth(2);
         kv.second->Draw("L same");
      }

      // PRA track hits, if any track is found
      if (patArr->GetEntries() > 0) {
         auto *pat = (AtPatternEvent *)patArr->At(0);
         int color = kRed + 1;
         for (auto &tr : pat->GetTrackCand()) {
            TGraph *gTr = new TGraph();
            int j = 0;
            for (auto &hp : tr.GetHitArray()) {
               auto pos = hp->GetPosition();
               gTr->SetPoint(j++, pos.Y(), pos.X());
            }
            gTr->SetMarkerStyle(24);
            gTr->SetMarkerColor(color);
            gTr->SetMarkerSize(0.9);
            gTr->SetLineColor(color);
            gTr->Draw("P same");
            ++color;
         }
      }

      // Pad-plane boundary in transposed display coords (h_x = phys-y,
      // h_y = phys-x).
      auto *box = new TLine();
      box->SetLineColor(kBlue + 1); box->SetLineWidth(2);
      box->DrawLine(padYmin, padXmin, padYmax, padXmin);
      box->DrawLine(padYmax, padXmin, padYmax, padXmax);
      box->DrawLine(padYmax, padXmax, padYmin, padXmax);
      box->DrawLine(padYmin, padXmax, padYmin, padXmin);

      // MC vertex marker (Prototype target is upstream at phys x ≈
      // -400 mm, phys y = 44 mm). Draw a small arrow at the entrance
      // face — bottom of the panel, pointing up.
      const double yMid = 0.5 * (padYmin + padYmax);
      auto *vtxArrow = new TArrow(yMid, padXmin - 15, yMid, padXmin + 5, 0.02, "|>");
      vtxArrow->SetLineColor(kGreen + 3); vtxArrow->SetFillColor(kGreen + 3);
      vtxArrow->SetLineWidth(2);
      vtxArrow->Draw();
   }

   TString png = TString(outPrefix) + ".png";
   TString pdf = TString(outPrefix) + ".pdf";
   c->SaveAs(png);
   c->SaveAs(pdf);
   std::cout << "Wrote " << outPrefix << ".{png,pdf}\n";

   std::cout << "\n=== events chosen (top 6 by primary-π MC point count) ===\n";
   for (size_t b = 0; b < picks.size(); ++b)
      printf("rank %zu: entry=%lld  nPi MC pts=%d  hits=%d\n",
             b + 1, picks[b].entry, picks[b].nPiPts, picks[b].nh);
}
