// Static (batch-safe) event display for the a1975 16C+p new-UKF workspace.
// Renders the AtEventH hits of one event to a PNG: XY pad-plane view, XZ and YZ
// side views, and a software-rendered 3D scatter (TGraph2D via TView, no OpenGL).
// Works on WSL where the interactive Eve/OpenGL viewer does not.
//
// Usage:
//   root -b -q 'event_png_UKF.C("run_0116_disp", 0)'       // draw event 0
//   root -b -q 'event_png_UKF.C("run_0116_disp", -1)'      // auto-pick the event with most hits
void event_png_UKF(TString inFile = "run_0116_disp", Long64_t evt = -1, Int_t scanN = 200)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kRainBow);
   gStyle->SetMarkerStyle(20);

   TFile *f = TFile::Open(inFile + ".root");
   if (!f || f->IsZombie()) {
      printf("ERROR: cannot open %s.root\n", inFile.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *evArr = nullptr;
   t->SetBranchAddress("AtEventH", &evArr);

   // Auto-pick the event with the most hits in the first scanN events
   if (evt < 0) {
      Int_t best = -1, bestN = -1;
      Long64_t lim = std::min((Long64_t)scanN, t->GetEntries());
      for (Long64_t i = 0; i < lim; ++i) {
         t->GetEntry(i);
         if (evArr->GetEntries() == 0)
            continue;
         AtEvent *e = (AtEvent *)evArr->At(0);
         if (e && (Int_t)e->GetNumHits() > bestN) {
            bestN = e->GetNumHits();
            best = i;
         }
      }
      evt = best;
      printf("Auto-picked event %lld with %d hits\n", evt, bestN);
   }

   t->GetEntry(evt);
   AtEvent *e = (AtEvent *)evArr->At(0);
   if (!e || e->GetNumHits() == 0) {
      printf("Event %lld has no hits\n", evt);
      return;
   }
   Int_t n = e->GetNumHits();

   TGraph *gXY = new TGraph();
   TGraph *gXZ = new TGraph(); // Z horizontal, X vertical
   TGraph *gYZ = new TGraph();
   TGraph2D *g3 = new TGraph2D();
   g3->SetName("g3");
   for (Int_t i = 0; i < n; ++i) {
      auto &h = e->GetHit(i);
      auto p = h.GetPosition();
      gXY->SetPoint(i, p.X(), p.Y());
      gXZ->SetPoint(i, p.Z(), p.X());
      gYZ->SetPoint(i, p.Z(), p.Y());
      g3->SetPoint(i, p.Z(), p.X(), p.Y());
   }

   auto deco = [](TGraph *g, const char *title, const char *xt, const char *yt) {
      g->SetTitle(Form("%s;%s;%s", title, xt, yt));
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.5);
      g->SetMarkerColor(kAzure + 2);
   };
   deco(gXY, "Pad plane (XY)", "X [mm]", "Y [mm]");
   deco(gXZ, "Top view (XZ)", "Z [mm] (drift)", "X [mm]");
   deco(gYZ, "Side view (YZ)", "Z [mm] (drift)", "Y [mm]");

   TCanvas *c = new TCanvas("c", "event", 1400, 1000);
   c->Divide(2, 2);
   c->cd(1);
   gXY->Draw("AP");
   c->cd(2);
   gXZ->Draw("AP");
   c->cd(3);
   gYZ->Draw("AP");
   c->cd(4);
   g3->SetTitle(Form("3D  (run %s, evt %lld, %d hits);Z [mm];X [mm];Y [mm]", inFile.Data(), evt, n));
   g3->SetMarkerStyle(20);
   g3->SetMarkerSize(0.6);
   g3->SetMarkerColor(kAzure + 2);
   g3->Draw("P0"); // software TView 3D scatter, no OpenGL

   TString png = Form("pipeline/plots/%s_evt%lld.png", inFile.Data(), evt);
   c->SaveAs(png);
   printf("Saved %s  (%d hits)\n", png.Data(), n);
}
