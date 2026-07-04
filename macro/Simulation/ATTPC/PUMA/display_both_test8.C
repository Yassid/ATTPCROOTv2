/// @file display_both_test8.C
/// @brief Static 2D-projection event display overlaying BOTH fitters on the PUMA
///        branch-8 sample (WSL-safe: no Eve/3D — saves PNGs).
///
/// Per event, three pads (XY, XZ, YZ) show:
///   grey dots        — reconstructed hits (AtPatternEvent point cloud)
///   black open circ  — MC truth points (attpcsim.root)
///   BLUE line+dots   — UKF fitted trajectory (smoothed);  blue star = UKF vertex
///   RED  line+dots   — GENFIT fitted trajectory (smoothed); red star = genfit vertex
///   green star       — MC truth vertex (0,0,75) mm
///
/// Run: root -b -q 'display_both_test8.C(8)'            // first 8 events
///      root -b -q 'display_both_test8.C(1, 42)'        // 1 event, starting at 42
void display_both_test8(int nShow = 8, int firstEvt = 0, TString digiFile = "./data/output_digi_both8.root",
                        TString simFile = "./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TFile fD(digiFile);
   TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);
   TTree *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("missing tree\n"); return; }

   TClonesArray *ukfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *gfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   tD->SetBranchAddress("AtPatternEvent", &patArr);
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   tS->SetBranchAddress("AtTpcPoint", &mcPts);

   // projection selector: view 0=XY, 1=XZ, 2=YZ
   auto proj = [](int view, double x, double y, double z, double &a, double &b) {
      if (view == 0) { a = x; b = y; }
      else if (view == 1) { a = x; b = z; }
      else { a = y; b = z; }
   };
   const char *vname[3] = {"XY", "XZ", "YZ"};
   const char *vax[3] = {"x [mm];y [mm]", "x [mm];z [mm]", "y [mm];z [mm]"};

   // draw one fitter's fitted trajectories into pad `view`, colour `col`
   auto drawFits = [&](TClonesArray *teArr, int view, int col) {
      if (teArr->GetEntries() == 0) return;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      for (const auto &ft : te->GetFittedTracks()) {
         const auto &sp = ft->GetSmoothedPositions();
         if (sp.size() >= 2) {
            auto *pl = new TPolyLine((int)sp.size());
            for (size_t i = 0; i < sp.size(); ++i) {
               double a, b; proj(view, sp[i].X(), sp[i].Y(), sp[i].Z(), a, b);
               pl->SetPoint((int)i, a, b);
            }
            pl->SetLineColor(col); pl->SetLineWidth(2); pl->Draw();
            auto *mk = new TPolyMarker((int)sp.size());
            for (size_t i = 0; i < sp.size(); ++i) {
               double a, b; proj(view, sp[i].X(), sp[i].Y(), sp[i].Z(), a, b);
               mk->SetPoint((int)i, a, b);
            }
            mk->SetMarkerColor(col); mk->SetMarkerStyle(20); mk->SetMarkerSize(0.5); mk->Draw();
         }
         // vertex (initialPositionXtr; fall back to GetVertex)
         const auto &pr = ft->GetTrackPropertiesStruct();
         double vx = pr.initialPositionXtr.X(), vy = pr.initialPositionXtr.Y(), vz = pr.initialPositionXtr.Z();
         if (std::abs(vx) < 1e-9 && std::abs(vy) < 1e-9 && std::abs(vz) < 1e-9) {
            const auto &v = ft->GetVertex(0); vx = v.X(); vy = v.Y(); vz = v.Z();
         }
         double a, b; proj(view, vx, vy, vz, a, b);
         auto *vm = new TMarker(a, b, 29); // star
         vm->SetMarkerColor(col); vm->SetMarkerSize(2.2); vm->Draw();
      }
   };

   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   int shown = 0;
   for (Long64_t e = firstEvt; e < nE && shown < nShow; ++e) {
      tD->GetEntry(e);
      tS->GetEntry(e);
      if (patArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patArr->At(0);
      auto &tracks = pat->GetTrackCand();
      if (tracks.empty()) continue;

      // collect reco hits and MC points into projection arrays
      std::vector<double> hx, hy, hz;
      for (auto &tr : tracks)
         for (auto &h : tr.GetHitArray()) { auto p = h->GetPosition(); hx.push_back(p.X()); hy.push_back(p.Y()); hz.push_back(p.Z()); }
      int nMC = mcPts->GetEntries();
      std::vector<double> mx(nMC), my(nMC), mz(nMC);
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k); mx[k] = mp->GetX() * 10; my[k] = mp->GetY() * 10; mz[k] = mp->GetZ() * 10; }

      auto *c = new TCanvas(Form("c%lld", e), Form("event %lld", e), 1500, 520);
      c->Divide(3, 1);
      for (int view = 0; view < 3; ++view) {
         c->cd(view + 1);
         // frame ranges: XY fixed to detector, XZ/YZ zoom near vertex
         double lo1 = -130, hi1 = 130, lo2, hi2;
         if (view == 0) { lo2 = -130; hi2 = 130; } else { lo2 = 40; hi2 = 130; }
         auto *fr = gPad->DrawFrame(lo1, lo2, hi1, hi2, Form("event %lld  %s;%s", e, vname[view], vax[view]));
         (void)fr;
         // reco hits (grey)
         if (!hx.empty()) {
            auto *gh = new TPolyMarker((int)hx.size());
            for (size_t i = 0; i < hx.size(); ++i) { double a, b; proj(view, hx[i], hy[i], hz[i], a, b); gh->SetPoint((int)i, a, b); }
            gh->SetMarkerColor(kGray + 1); gh->SetMarkerStyle(20); gh->SetMarkerSize(0.4); gh->Draw();
         }
         // MC truth points (black open circles)
         if (nMC) {
            auto *gm = new TPolyMarker(nMC);
            for (int k = 0; k < nMC; ++k) { double a, b; proj(view, mx[k], my[k], mz[k], a, b); gm->SetPoint(k, a, b); }
            gm->SetMarkerColor(kBlack); gm->SetMarkerStyle(24); gm->SetMarkerSize(0.5); gm->Draw();
         }
         drawFits(ukfArr, view, kBlue + 1);
         drawFits(gfArr, view, kRed + 1);
         // truth vertex (green star)
         double a, b; proj(view, 0., 0., 75., a, b);
         auto *tv = new TMarker(a, b, 29); tv->SetMarkerColor(kGreen + 2); tv->SetMarkerSize(2.0); tv->Draw();
      }
      // legend on first pad
      c->cd(1);
      auto *lg = new TLegend(0.60, 0.75, 0.98, 0.98);
      lg->SetTextSize(0.032);
      auto *lu = new TMarker(); lu->SetMarkerColor(kBlue + 1); lu->SetMarkerStyle(20);
      auto *lgf = new TMarker(); lgf->SetMarkerColor(kRed + 1); lgf->SetMarkerStyle(20);
      lg->AddEntry(lu, "UKF fit + vertex", "lp");
      lg->AddEntry(lgf, "genfit fit + vertex", "lp");
      lg->AddEntry((TObject *)nullptr, "grey=hits blk=MC grn=truth vtx", "");
      lg->Draw();

      TString png = Form("./data/recon_evt_%02lld.png", e);
      c->SaveAs(png);
      printf("saved %s  (UKF fits=%d, genfit fits=%d)\n", png.Data(),
             ukfArr->GetEntries() ? (int)((AtTrackingEvent *)ukfArr->At(0))->GetFittedTracks().size() : 0,
             gfArr->GetEntries() ? (int)((AtTrackingEvent *)gfArr->At(0))->GetFittedTracks().size() : 0);
      shown++;
   }
   printf("done: %d event PNGs in ./data/\n", shown);
}
