/// @file draw_miss_test8.C
/// @brief Why does one pion leave no hits in PRA=1 events? Draw both primaries'
///        helices (computed from their MC momenta) over the detector, wide XY view.
///        The detected pion's helix should overlay its MC hits (convention check);
///        the other helix shows where the undetected pion actually goes.
/// Run: root -b -q 'draw_miss_test8.C(12)'   // or a 2x2 of PRA=1 events
void draw_miss_test8(int ev0 = 12, int ev1 = 34, int ev2 = 42, int ev3 = 51,
                     TString simFile = "./data/attpcsim.root", TString outPng = "./data/miss_geom.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const double Bz = 4.0, kRin = 62.9, kRout = 121.1;

   TFile fS(simFile); TTree *tS = (TTree *)fS.Get("cbmsim");
   TClonesArray *mcT = new TClonesArray("AtMCTrack"); tS->SetBranchAddress("MCTrack", &mcT);
   TClonesArray *mcP = new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint", &mcP);

   auto *c = new TCanvas("miss", "missing pion geometry", 1200, 1200);
   c->Divide(2, 2);
   int evs[4] = {ev0, ev1, ev2, ev3};

   for (int pad = 0; pad < 4; ++pad) {
      c->cd(pad + 1);
      Long64_t e = evs[pad];
      tS->GetEntry(e);
      gPad->DrawFrame(-700, -700, 700, 700, Form("evt %lld  (wide XY);x [mm];y [mm]", e));
      // detector rings
      for (double r : {kRin, kRout}) { auto *el = new TEllipse(0, 0, r, r); el->SetFillStyle(0); el->SetLineColor(kGray + 2); el->Draw(); }

      // detected pion's MC hits (grey) — convention validation
      for (int k = 0; k < mcP->GetEntries(); ++k) { auto *mp = (AtMCPoint *)mcP->At(k); if (mp->GetTrackID() < 0) continue;
         auto *m = new TMarker(mp->GetX() * 10, mp->GetY() * 10, 20); m->SetMarkerColor(kGray + 1); m->SetMarkerSize(0.4); m->Draw(); }

      // both primaries' helices from MC momenta
      int colr[2] = {kBlue + 1, kRed + 1};
      for (int t = 0; t < std::min(2, mcT->GetEntries()); ++t) {
         auto *mt = (AtMCTrack *)mcT->At(t);
         double px = mt->GetPx(), py = mt->GetPy();           // GeV/c
         double x0 = mt->GetStartX() * 10, y0 = mt->GetStartY() * 10; // mm
         double pt = std::hypot(px, py);
         double R = pt / (0.299792458 * Bz) * 1000.0;          // mm, |q|=1
         int qs = (mt->GetPdgCode() > 0) ? +1 : -1;            // pi+ = +1
         double phip = std::atan2(py, px);
         double phic = phip - qs * M_PI / 2.0;                 // center direction
         double cx = x0 + R * std::cos(phic), cy = y0 + R * std::sin(phic);
         double a0 = std::atan2(y0 - cy, x0 - cx);
         // sweep one full turn; clip drawing when it leaves a generous window
         const int N = 400; auto *pl = new TPolyLine();
         int np = 0;
         for (int i = 0; i <= N; ++i) {
            double s = (2 * M_PI * R) * i / N;                 // arc length
            double a = a0 + qs * (s / R);
            double x = cx + R * std::cos(a), y = cy + R * std::sin(a);
            if (std::hypot(x, y) > 680) break;                 // left the view -> stop (particle long gone)
            pl->SetPoint(np++, x, y);
         }
         pl->SetLineColor(colr[t]); pl->SetLineWidth(2); pl->Draw();
         auto *sm = new TMarker(x0, y0, 29); sm->SetMarkerColor(kGreen + 2); sm->SetMarkerSize(1.6); sm->Draw();
      }
   }
   c->SaveAs(outPng);
   printf("saved %s\n", outPng.Data());
}
