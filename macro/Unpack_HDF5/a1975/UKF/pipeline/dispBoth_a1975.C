/// @file dispBoth_a1975.C
/// @brief Static (PNG) dual-fitter event display — the WSL-friendly counterpart of
/// run_eve_both_a1975.C (no Eve/OpenGL). For each event it overlays both Kalman
/// filters' fitted trajectories on the reconstructed hits in three projections:
///   XY (pad plane), ZY, ZX.   gray = hits, red = genfit, blue dashed = UKF.
///
/// Frame: hits are digi-frame (z = drift), fitted smoothed positions are lab-frame
/// (z_lab = ZPadPlane - z_digi), mapped back via z_digi = ZPadPlane - z_lab.
///
/// Modes (the `kind` argument):
///   ""         single event `iEvent` (-1 = first event that has both fits)
///   "multi"    first `nWanted` events with >= 2 PRA tracks
///   "backward" first `nWanted` events with a backward PRA track (GeoTheta > 90 deg)
///
///   root -l -b -q 'dispBoth_a1975.C("run_0106","multi",4,"/mnt/f/a1975/reco/","/tmp/")'
///   root -l -b -q 'dispBoth_a1975.C("run_0106","backward",4,"/mnt/f/a1975/reco/","/tmp/")'
///   root -l -b -q 'dispBoth_a1975.C("run_0106","",1,"/mnt/f/a1975/reco/","/tmp/",146)'

namespace {
double labToDigiZ(double zLab, double zPad) { return zPad - zLab; }

// Render one entry to a PNG. trees/branches already set up by the caller.
void renderEvent(Long64_t iEvent, TTree *tr, TTree *tf, TClonesArray *evtArr, TClonesArray *gfArr,
                 TClonesArray *ukArr, TString tag, TString fitDir, TString runName, double zPad)
{
   tr->GetEntry(iEvent);
   tf->GetEntry(iEvent);

   std::vector<double> hx, hy, hz;
   if (auto *ev = (AtEvent *)evtArr->At(0))
      for (const auto &hit : ev->GetHitArray()) {
         auto p = hit.GetPosition();
         hx.push_back(p.X()); hy.push_back(p.Y()); hz.push_back(p.Z());
      }

   auto buildPolylines = [&](TClonesArray *arr, std::vector<std::vector<ROOT::Math::XYZPoint>> &out) {
      auto *te = (AtTrackingEvent *)arr->At(0);
      if (!te) return;
      for (const auto &ftp : te->GetFittedTracks()) {
         std::vector<ROOT::Math::XYZPoint> poly;
         const auto &pV = ftp->GetTrackPropertiesStruct().initialPositionXtr;
         if (pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0)
            poly.emplace_back(pV.X(), pV.Y(), labToDigiZ(pV.Z(), zPad));
         for (auto &sp : ftp->GetSmoothedPositions())
            poly.emplace_back(sp.X(), sp.Y(), labToDigiZ(sp.Z(), zPad));
         if (poly.size() > 1) out.push_back(std::move(poly));
      }
   };
   std::vector<std::vector<ROOT::Math::XYZPoint>> gfPoly, ukPoly;
   buildPolylines(gfArr, gfPoly);
   buildPolylines(ukArr, ukPoly);
   std::cout << "  [ev " << iEvent << "] hits=" << hx.size() << "  genfit=" << gfPoly.size()
             << "  UKF=" << ukPoly.size() << "\n";

   gStyle->SetOptStat(0);
   auto *c = new TCanvas(TString::Format("c_%lld", iEvent), "dual fit", 1500, 520);
   c->Divide(3, 1);
   auto coord = [](const ROOT::Math::XYZPoint &p, int proj, int axis) {
      if (proj == 0) return axis == 0 ? p.X() : p.Y();
      if (proj == 1) return axis == 0 ? p.Z() : p.Y();
      return axis == 0 ? p.Z() : p.X();
   };
   const char *titles[3] = {"XY (pad plane);x [mm];y [mm]", "ZY;z_{drift} [mm];y [mm]", "ZX;z_{drift} [mm];x [mm]"};
   std::vector<TObject *> keep;
   for (int proj = 0; proj < 3; ++proj) {
      c->cd(proj + 1);
      auto *gH = new TGraph();
      for (size_t i = 0; i < hx.size(); ++i) {
         ROOT::Math::XYZPoint p(hx[i], hy[i], hz[i]);
         gH->SetPoint(gH->GetN(), coord(p, proj, 0), coord(p, proj, 1));
      }
      gH->SetTitle(TString::Format("%s  %s", tag.Data(), titles[proj]));
      gH->SetMarkerStyle(1);
      gH->SetMarkerColor(kGray + 1);
      gH->Draw("AP");
      keep.push_back(gH);
      auto drawPolys = [&](std::vector<std::vector<ROOT::Math::XYZPoint>> &polys, Color_t col, Style_t style) {
         for (auto &poly : polys) {
            auto *g = new TGraph();
            for (auto &p : poly) g->SetPoint(g->GetN(), coord(p, proj, 0), coord(p, proj, 1));
            g->SetLineColor(col); g->SetLineWidth(2); g->SetLineStyle(style);
            g->Draw("L SAME");
            keep.push_back(g);
         }
      };
      drawPolys(gfPoly, kRed + 1, 1);
      drawPolys(ukPoly, kAzure + 1, 2);
   }
   TString png = TString::Format("%s%s_both_%s_ev%lld.png", fitDir.Data(), runName.Data(), tag.Data(), iEvent);
   c->SaveAs(png);
   std::cout << "  wrote " << png << "\n";
}
} // namespace

void dispBoth_a1975(TString runName = "run_0106", TString kind = "", Int_t nWanted = 4,
                    TString recoDir = "/mnt/f/a1975/reco/", TString fitDir = "/tmp/", Long64_t iEvent = -1,
                    Double_t zPadPlane = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");

   TFile *fr = TFile::Open(recoDir + runName + "_reco.root");
   TFile *ff = TFile::Open(fitDir + runName + "_both.root");
   if (!fr || !ff) { std::cout << "Missing input file(s)\n"; return; }
   auto *tr = (TTree *)fr->Get("cbmsim");
   auto *tf = (TTree *)ff->Get("cbmsim");

   auto *evtArr = new TClonesArray("AtEvent");
   auto *patArr = new TClonesArray("AtPatternEvent");
   auto *gfArr = new TClonesArray("AtTrackingEvent");
   auto *ukArr = new TClonesArray("AtTrackingEvent");
   tr->SetBranchAddress("AtEventCorrected", &evtArr);
   tr->SetBranchAddress("AtPatternEvent", &patArr);
   tf->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   tf->SetBranchAddress("AtTrackingEventUKF", &ukArr);

   Long64_t n = std::min(tr->GetEntries(), tf->GetEntries());

   // ─── single-event mode ───
   if (kind == "") {
      if (iEvent < 0) {
         for (Long64_t i = 0; i < n; ++i) {
            tf->GetEntry(i);
            auto *ge = (AtTrackingEvent *)gfArr->At(0);
            auto *ue = (AtTrackingEvent *)ukArr->At(0);
            if (ge && ue && !ge->GetFittedTracks().empty() && !ue->GetFittedTracks().empty()) { iEvent = i; break; }
         }
         if (iEvent < 0) { std::cout << "No event with both fits found.\n"; return; }
      }
      renderEvent(iEvent, tr, tf, evtArr, gfArr, ukArr, "ev", fitDir, runName, zPadPlane);
      return;
   }

   // ─── scan modes: collect matching event indices from the PRA tracks ───
   std::vector<Long64_t> hits;
   for (Long64_t i = 0; i < n && (int)hits.size() < nWanted; ++i) {
      tr->GetEntry(i);
      auto *pe = (AtPatternEvent *)patArr->At(0);
      if (!pe) continue;
      auto &tracks = pe->GetTrackCand();
      bool match = false;
      if (kind == "multi") {
         match = tracks.size() >= 2;
      } else if (kind == "backward") {
         for (auto &t : tracks)
            if (t.GetGeoTheta() * TMath::RadToDeg() > 90.0) { match = true; break; }
      }
      if (match) hits.push_back(i);
   }
   std::cout << "kind='" << kind << "' : found " << hits.size() << " matching events: ";
   for (auto i : hits) std::cout << i << " ";
   std::cout << "\n";
   for (auto i : hits)
      renderEvent(i, tr, tf, evtArr, gfArr, ukArr, kind, fitDir, runName, zPadPlane);
}
