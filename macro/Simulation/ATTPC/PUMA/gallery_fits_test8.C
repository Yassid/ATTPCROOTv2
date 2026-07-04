/// @file gallery_fits_test8.C
/// @brief Gallery of reconstructed branch-8 events. Each event gets TWO stacked pads:
///        XY (top) and a Z projection (bottom, "XZ" default or "YZ"), with BOTH fitters
///        overlaid (WSL-safe PNG). grey = hits, black open = MC truth,
///        BLUE = UKF fit+vertex, RED = genfit fit+vertex, green star = truth vertex.
///        Pad title shows PRA candidate count and each fitter's fit count.
///
/// Run: root -b -q gallery_fits_test8.C                       // 4x2 events (8), XY+XZ
///      root -b -q 'gallery_fits_test8.C(4, 2, 0, "YZ")'      // YZ instead of XZ
///      root -b -q 'gallery_fits_test8.C(4, 3, 20)'           // 4x3 events from #20
void gallery_fits_test8(int nEvCol = 4, int nEvRow = 2, int firstEvt = 0, TString zView = "XZ",
                        TString digiFile = "./data/output_digi_both8.root",
                        TString simFile = "./data/attpcsim.root", TString outPng = "./data/gallery_fits.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const int vz = (zView == "YZ") ? 2 : 1; // 1 = XZ, 2 = YZ

   TFile fD(digiFile); TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);  TTree *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("missing tree\n"); return; }

   TClonesArray *ukfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *gfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   tD->SetBranchAddress("AtPatternEvent", &patArr);
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   tS->SetBranchAddress("AtTpcPoint", &mcPts);

   // view: 0=XY (pad plane); 1=Z-X, 2=Z-Y (drift time z runs LEFT-TO-RIGHT, so z is the
   // horizontal axis and the transverse coord is vertical).
   auto proj = [](int view, double x, double y, double z, double &a, double &b) {
      if (view == 0) { a = x; b = y; } else if (view == 1) { a = z; b = x; } else { a = z; b = y; }
   };
   // PUMA detector geometry (mm): annular gas pad plane + drift volume along z.
   const double kRin = 62.9, kRout = 121.1, kZpad = 0.0, kZcath = 300.0;
   auto drawDetector = [&](int view) {
      if (view == 0) { // pad plane: inner + outer rings
         for (double r : {kRin, kRout}) {
            auto *el = new TEllipse(0, 0, r, r); el->SetFillStyle(0);
            el->SetLineColor(kGray + 1); el->SetLineStyle(r == kRin ? 2 : 1); el->Draw();
         }
      } else { // drift volume box: pad plane (z=0) .. cathode (z=300), radial walls at +/-Rout
         auto line = [](double x1, double y1, double x2, double y2, int col, int sty) {
            auto *l = new TLine(x1, y1, x2, y2); l->SetLineColor(col); l->SetLineStyle(sty); l->SetLineWidth(2); l->Draw();
         };
         // NOTE: only the OUTER radius is a valid bound in this projection — r<=Rout
         // implies |x|<=Rout. The inner annular hole is RADIAL (r<Rin), which does NOT
         // map to an |x|<Rin band (a hit at small x, large y is still in the gas), so
         // no inner lines are drawn here (they would falsely imply a forbidden strip).
         line(kZpad, -kRout, kZpad, kRout, kAzure + 2, 1);   // pad plane (readout) — emphasized
         line(kZcath, -kRout, kZcath, kRout, kGray + 1, 1);  // cathode
         line(kZpad, kRout, kZcath, kRout, kGray + 1, 1);    // outer wall +R
         line(kZpad, -kRout, kZcath, -kRout, kGray + 1, 1);  // outer wall -R
      }
   };
   auto drawFits = [&](TClonesArray *teArr, int view, int col) {
      if (teArr->GetEntries() == 0) return;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      for (const auto &ft : te->GetFittedTracks()) {
         const auto &sp = ft->GetSmoothedPositions();
         if (sp.size() >= 2) {
            auto *pl = new TPolyLine((int)sp.size());
            for (size_t i = 0; i < sp.size(); ++i) { double a, b; proj(view, sp[i].X(), sp[i].Y(), sp[i].Z(), a, b); pl->SetPoint((int)i, a, b); }
            pl->SetLineColor(col); pl->SetLineWidth(2); pl->Draw();
         }
         const auto &pr = ft->GetTrackPropertiesStruct();
         double vx = pr.initialPositionXtr.X(), vy = pr.initialPositionXtr.Y(), vzz = pr.initialPositionXtr.Z();
         if (std::abs(vx) < 1e-9 && std::abs(vy) < 1e-9 && std::abs(vzz) < 1e-9) { const auto &v = ft->GetVertex(0); vx = v.X(); vy = v.Y(); vzz = v.Z(); }
         double a, b; proj(view, vx, vy, vzz, a, b);
         auto *vm = new TMarker(a, b, 29); vm->SetMarkerColor(col); vm->SetMarkerSize(1.5); vm->Draw();
      }
   };
   // draw one projection pad for the current (already-loaded) event
   auto drawPad = [&](int view, Long64_t e, const char *title) {
      gPad->SetMargin(0.09, 0.02, 0.02, 0.08);
      if (view == 0) // XY: full pad plane
         gPad->DrawFrame(-130, -130, 130, 130, title);
      else           // Z view: drift z (0..cathode) horizontal, transverse vertical
         gPad->DrawFrame(-15, -130, kZcath + 15, 130, title);
      drawDetector(view);
      auto *pat = (AtPatternEvent *)patArr->At(0);
      std::vector<double> ha, hb;
      for (auto &tr : pat->GetTrackCand()) for (auto &h : tr.GetHitArray()) { auto p = h->GetPosition(); double a, b; proj(view, p.X(), p.Y(), p.Z(), a, b); ha.push_back(a); hb.push_back(b); }
      if (!ha.empty()) { auto *g = new TPolyMarker((int)ha.size(), ha.data(), hb.data()); g->SetMarkerColor(kGray + 1); g->SetMarkerStyle(20); g->SetMarkerSize(0.35); g->Draw(); }
      int nMC = mcPts->GetEntries();
      if (nMC) { std::vector<double> ma, mb; for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k); double a, b; proj(view, mp->GetX() * 10, mp->GetY() * 10, mp->GetZ() * 10, a, b); ma.push_back(a); mb.push_back(b); }
         auto *g = new TPolyMarker((int)ma.size(), ma.data(), mb.data()); g->SetMarkerColor(kBlack); g->SetMarkerStyle(24); g->SetMarkerSize(0.4); g->Draw(); }
      drawFits(ukfArr, view, kBlue + 1);
      drawFits(gfArr, view, kRed + 1);
      double a, b; proj(view, 0., 0., 75., a, b);
      auto *tv = new TMarker(a, b, 29); tv->SetMarkerColor(kGreen + 2); tv->SetMarkerSize(1.3); tv->Draw();
   };

   auto *c = new TCanvas("gallery", "gallery of fits", 300 * nEvCol, 300 * nEvRow); // 2 stacked pads/event
   c->Divide(nEvCol, 2 * nEvRow, 0.001, 0.001);
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());

   int placed = 0, e = firstEvt;
   const int nEvents = nEvCol * nEvRow;
   while (placed < nEvents && e < nE) {
      tD->GetEntry(e); tS->GetEntry(e);
      auto *pat = (patArr->GetEntries()) ? (AtPatternEvent *)patArr->At(0) : nullptr;
      if (!pat || pat->GetTrackCand().empty()) { ++e; continue; }

      int nPRA = pat->GetTrackCand().size();
      int nU = ukfArr->GetEntries() ? (int)((AtTrackingEvent *)ukfArr->At(0))->GetFittedTracks().size() : 0;
      int nG = gfArr->GetEntries() ? (int)((AtTrackingEvent *)gfArr->At(0))->GetFittedTracks().size() : 0;

      int evCol = placed % nEvCol, evRow = placed / nEvCol;
      int padXY = (2 * evRow) * nEvCol + evCol + 1;      // top pad of the event's stack
      int padZ = (2 * evRow + 1) * nEvCol + evCol + 1;   // bottom pad
      const char *tCoord = (vz == 2) ? "y" : "x";
      c->cd(padXY); drawPad(0, e, Form("evt %lld  XY  PRA:%d U:%d g:%d;x [mm];y [mm]", e, nPRA, nU, nG));
      c->cd(padZ);  drawPad(vz, e, Form("evt %lld  drift(z) vs %s;z [mm]  (drift / time);%s [mm]", e, tCoord, tCoord));
      ++placed; ++e;
   }
   c->SaveAs(outPng);
   printf("gallery: %d events (XY + %s) -> %s\n", placed, zView.Data(), outPng.Data());
}
