/// @file inspect_fit.C
/// @brief 2D-projection inspection viewer for the PUMA UKF fit (mirrors
///        16C_pp/run_ukf_display.C's TCanvas approach — no Eve, no OpenGL,
///        no unit-convention guessing).
///
/// For one event:
///  - Pattern-event hits per track (one color per track)
///  - Per-fit smoothed UKF positions (red open circles)
///  - Per-fit back-extrapolated vertex (red star)
///  - Per-fit kinematic-direction stub from the vertex (50 mm dashed arrow)
///  - MC truth start vertex (black star) and primary-direction stub (gray)
///
/// Three pads: XY, XZ, YZ. Whatever convention is wrong (Z flip, helix
/// curvature side, theta/phi) is visible immediately by comparing the
/// fitted stub to the actual hit cloud and to the truth direction.
///
/// Run: root -b -q 'inspect_fit.C(0)'   // inspect event 0
///      root -l 'inspect_fit.C(7)'      // open canvas interactively

void inspect_fit(int event = 0,
                 TString digiFile = "./data/output_digi.root",
                 TString simFile = "./data/attpcsim.root",
                 TString outPng = "")
{
   auto *fD = TFile::Open(digiFile);
   auto *fS = TFile::Open(simFile);
   if (!fD || fD->IsZombie() || !fS || fS->IsZombie()) {
      std::cerr << "Cannot open inputs\n"; return;
   }
   auto *tD = (TTree *)fD->Get("cbmsim");
   auto *tS = (TTree *)fS->Get("cbmsim");
   if (!tD || !tS) { std::cerr << "Missing cbmsim tree\n"; return; }

   TClonesArray *peArr = nullptr, *teArr = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tD->SetBranchAddress("AtPatternEvent", &peArr);
   if (tD->GetBranch("AtTrackingEvent"))
      tD->SetBranchAddress("AtTrackingEvent", &teArr);
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   if (event >= tD->GetEntries() || event >= tS->GetEntries()) {
      std::cerr << "Event " << event << " out of range\n"; return;
   }
   tD->GetEntry(event);
   tS->GetEntry(event);

   if (!peArr || peArr->GetEntries() == 0) {
      std::cerr << "Event " << event << ": no AtPatternEvent\n"; return;
   }
   auto *pe = (AtPatternEvent *)peArr->At(0);
   auto &tracks = pe->GetTrackCand();
   std::cout << "Event " << event << ": " << tracks.size()
             << " PRA tracks, " << mcPts->GetEntries() << " MC points\n";

   // ---- Build per-pad data buckets ----
   const int nProj = 3;       // XY, XZ, YZ
   const Color_t trackColors[] = {kBlue + 1, kGreen + 2, kOrange + 7,
                                  kViolet + 1, kCyan + 2, kYellow - 2,
                                  kAzure + 7};
   auto col = [&](int i) { return trackColors[i % 7]; };

   // Helper: draw a frame pad and overlay all elements ------------------
   auto pickAxes = [&](int proj, double X, double Y, double Z, double &a, double &b) {
      // proj=0: XY,  proj=1: XZ (X horizontal, Z vertical),  proj=2: YZ
      if (proj == 0) { a = X; b = Y; }
      else if (proj == 1) { a = Z; b = X; }
      else              { a = Z; b = Y; }
   };
   auto axisLabels = [&](int proj) -> std::pair<const char *, const char *> {
      if (proj == 0) return {"X [mm]", "Y [mm]"};
      if (proj == 1) return {"Z [mm]", "X [mm]"};
      return {"Z [mm]", "Y [mm]"};
   };

   auto *c = new TCanvas("c_inspect_fit",
                         TString::Format("PUMA inspect_fit  event %d", event),
                         1700, 600);
   c->Divide(3, 1);

   for (int proj = 0; proj < nProj; ++proj) {
      c->cd(proj + 1);
      gPad->SetGrid();

      // Axis ranges sized to the pad-plane R + drift length, with margins.
      double aLo, aHi, bLo, bHi;
      if (proj == 0) { aLo = bLo = -150; aHi = bHi = 150; }
      else if (proj == 1) { aLo = -50; aHi = 350; bLo = -150; bHi = 150; }
      else                { aLo = -50; aHi = 350; bLo = -150; bHi = 150; }

      auto labels = axisLabels(proj);
      auto *frame = gPad->DrawFrame(aLo, bLo, aHi, bHi,
                                    TString::Format("event %d  %s vs %s ;%s;%s",
                                                    event, labels.second, labels.first,
                                                    labels.first, labels.second));
      (void)frame;

      // --- MC truth points (small grey dots) ---
      auto *gMC = new TGraph();
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         double a, b;
         pickAxes(proj, mp->GetX() * 10.0, mp->GetY() * 10.0, mp->GetZ() * 10.0, a, b);
         gMC->SetPoint(gMC->GetN(), a, b);
      }
      gMC->SetMarkerStyle(7);
      gMC->SetMarkerColor(kGray + 1);
      gMC->Draw("P SAME");

      // --- PRA hits per track ---
      for (size_t i = 0; i < tracks.size(); ++i) {
         auto &tr = tracks[i];
         auto &hits = tr.GetHitArray();
         if (hits.empty()) continue;
         auto *g = new TGraph();
         for (auto &h : hits) {
            const auto &p = h->GetPosition();
            double a, b;
            pickAxes(proj, p.X(), p.Y(), p.Z(), a, b);
            g->SetPoint(g->GetN(), a, b);
         }
         g->SetMarkerStyle(20);
         g->SetMarkerSize(0.5);
         g->SetMarkerColor(col((int)i));
         g->Draw("P SAME");
      }

      // --- Fitted-track elements ---
      if (teArr && teArr->GetEntries() > 0) {
         auto *te = (AtTrackingEvent *)teArr->At(0);
         const auto &fits = te->GetFittedTracks();
         for (size_t i = 0; i < fits.size(); ++i) {
            const auto &ft = *fits[i];
            const auto &props = ft.GetTrackPropertiesStruct();
            const auto &kin = ft.GetKinematics();
            const auto &pinfo = ft.GetParticleInfo();
            const auto &pV = props.initialPositionXtr;
            const auto &p0 = props.initialPosition;

            // Vertex (red star)
            if (pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0) {
               double a, b;
               pickAxes(proj, pV.X(), pV.Y(), pV.Z(), a, b);
               auto *m = new TMarker(a, b, 29);
               m->SetMarkerSize(2.0);
               m->SetMarkerColor(kRed + 1);
               m->Draw();
            }

            // First cluster (small red triangle)
            if (p0.X() != 0 || p0.Y() != 0 || p0.Z() != 0) {
               double a, b;
               pickAxes(proj, p0.X(), p0.Y(), p0.Z(), a, b);
               auto *m = new TMarker(a, b, 22);
               m->SetMarkerSize(1.2);
               m->SetMarkerColor(kRed + 1);
               m->Draw();
            }

            // Smoothed UKF positions (red open circles, joined by line)
            const auto &smoothed = props.fSmoothedPositions;
            if (!smoothed.empty()) {
               auto *gSm = new TGraph();
               for (const auto &sp : smoothed) {
                  double a, b;
                  pickAxes(proj, sp.X(), sp.Y(), sp.Z(), a, b);
                  gSm->SetPoint(gSm->GetN(), a, b);
               }
               gSm->SetMarkerStyle(24);
               gSm->SetMarkerSize(0.7);
               gSm->SetLineColor(kRed + 1);
               gSm->SetLineWidth(2);
               gSm->SetMarkerColor(kRed + 1);
               gSm->Draw("LP SAME");
            }

            // Kinematic direction stub from vertex (50 mm dashed)
            if (kin.kineticEnergy > 0 && pinfo.mass > 0
                && std::isfinite(kin.theta) && std::isfinite(kin.phi)
                && (pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0)) {
               const double L = 50.0; // mm
               const double ux = std::sin(kin.theta) * std::cos(kin.phi);
               const double uy = std::sin(kin.theta) * std::sin(kin.phi);
               const double uz = std::cos(kin.theta);
               double a0, b0, a1, b1;
               pickAxes(proj, pV.X(), pV.Y(), pV.Z(), a0, b0);
               pickAxes(proj, pV.X() + L * ux, pV.Y() + L * uy, pV.Z() + L * uz, a1, b1);
               auto *l = new TLine(a0, b0, a1, b1);
               l->SetLineColor(kRed + 1);
               l->SetLineWidth(2);
               l->SetLineStyle(2);
               l->Draw();
            }
         }
      }

      // --- MC primary truth: vertex (black star) + initial-momentum stub ---
      if (mcTrks->GetEntries() > 0) {
         for (int k = 0; k < std::min(3, mcTrks->GetEntries()); ++k) {
            auto *mt = (AtMCTrack *)mcTrks->At(k);
            int pdg = mt->GetPdgCode();
            // Only primaries with |pdg| in {321, 211, 2212}
            if (std::abs(pdg) != 321 && std::abs(pdg) != 211 && std::abs(pdg) != 2212)
               continue;
            double vx = mt->GetStartX() * 10.0; // cm → mm
            double vy = mt->GetStartY() * 10.0;
            double vz = mt->GetStartZ() * 10.0;
            double a, b;
            pickAxes(proj, vx, vy, vz, a, b);
            auto *m = new TMarker(a, b, 29);
            m->SetMarkerSize(2.0);
            m->SetMarkerColor(kBlack);
            m->Draw();

            // truth direction
            double px = mt->GetPx(), py = mt->GetPy(), pz = mt->GetPz();
            double pmag = std::sqrt(px*px + py*py + pz*pz);
            if (pmag > 0) {
               const double L = 50.0;
               double a1, b1;
               pickAxes(proj, vx + L * px / pmag, vy + L * py / pmag, vz + L * pz / pmag, a1, b1);
               auto *l = new TLine(a, b, a1, b1);
               l->SetLineColor(kBlack);
               l->SetLineWidth(2);
               l->SetLineStyle(3);
               l->Draw();
            }
         }
      }
   }

   if (outPng.Length() > 0) {
      c->SaveAs(outPng);
      std::cout << "Saved " << outPng << "\n";
   }
}
