// Batch static event display for the a1975 16C+p new-UKF workspace (WSL-safe, no GL).
// Produces:
//   1) a contact sheet (grid of XY pad-plane views) over [evStart, evEnd) for a quick survey
//   2) individual detailed quad PNGs (XY/XZ/YZ + software 3D) for each event in the range
// Output goes into a "<file>_pngs/" subfolder.
//
// Usage:  root -b -q 'event_png_batch_UKF.C("run_0116_disp", 0, 24)'
void event_png_batch_UKF(TString inFile = "run_0116_disp", Long64_t evStart = 0, Long64_t evEnd = 24)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetMarkerStyle(20);

   TFile *f = TFile::Open(inFile + ".root");
   if (!f || f->IsZombie()) {
      printf("ERROR: cannot open %s.root\n", inFile.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *evArr = nullptr;
   t->SetBranchAddress("AtEventH", &evArr);

   if (evEnd > t->GetEntries())
      evEnd = t->GetEntries();
   Long64_t nEv = evEnd - evStart;
   if (nEv <= 0) {
      printf("empty range\n");
      return;
   }

   TString outDir = inFile + "_pngs";
   gSystem->mkdir(outDir, kTRUE);

   // --- contact sheet grid ---
   Int_t ncol = 6;
   Int_t nrow = (nEv + ncol - 1) / ncol;
   TCanvas *cs = new TCanvas("cs", "contact", 300 * ncol, 250 * nrow);
   cs->Divide(ncol, nrow, 0.001, 0.001);

   Int_t pad = 1;
   for (Long64_t ev = evStart; ev < evEnd; ++ev, ++pad) {
      t->GetEntry(ev);
      AtEvent *e = (evArr->GetEntries() > 0) ? (AtEvent *)evArr->At(0) : nullptr;
      Int_t n = e ? e->GetNumHits() : 0;

      // points for this event
      TGraph *gXY = new TGraph();
      TGraph *gXZ = new TGraph();
      TGraph *gYZ = new TGraph();
      TGraph2D *g3 = new TGraph2D();
      g3->SetName(Form("g3_%lld", ev));
      for (Int_t i = 0; i < n; ++i) {
         auto &h = e->GetHit(i);
         auto p = h.GetPosition();
         gXY->SetPoint(i, p.X(), p.Y());
         gXZ->SetPoint(i, p.Z(), p.X());
         gYZ->SetPoint(i, p.Z(), p.Y());
         g3->SetPoint(i, p.Z(), p.X(), p.Y());
      }

      // contact-sheet pad: XY view
      cs->cd(pad);
      gXY->SetTitle(Form("evt %lld  (%d hits);X [mm];Y [mm]", ev, n));
      gXY->SetMarkerStyle(20);
      gXY->SetMarkerSize(0.4);
      gXY->SetMarkerColor(kAzure + 2);
      if (n > 0)
         gXY->Draw("AP");

      // individual detailed quad
      TCanvas *cd = new TCanvas(Form("cd_%lld", ev), "evt", 1400, 1000);
      cd->Divide(2, 2);
      auto deco = [](TGraph *g, const char *ti, const char *xt, const char *yt) {
         g->SetTitle(Form("%s;%s;%s", ti, xt, yt));
         g->SetMarkerStyle(20);
         g->SetMarkerSize(0.5);
         g->SetMarkerColor(kAzure + 2);
      };
      deco(gXY, "Pad plane (XY)", "X [mm]", "Y [mm]");
      deco(gXZ, "Top view (XZ)", "Z [mm] (drift)", "X [mm]");
      deco(gYZ, "Side view (YZ)", "Z [mm] (drift)", "Y [mm]");
      cd->cd(1);
      if (n > 0)
         gXY->DrawClone("AP");
      cd->cd(2);
      if (n > 0)
         gXZ->Draw("AP");
      cd->cd(3);
      if (n > 0)
         gYZ->Draw("AP");
      cd->cd(4);
      g3->SetTitle(Form("3D evt %lld (%d hits);Z [mm];X [mm];Y [mm]", ev, n));
      g3->SetMarkerStyle(20);
      g3->SetMarkerSize(0.6);
      g3->SetMarkerColor(kAzure + 2);
      if (n > 0)
         g3->Draw("P0");
      cd->SaveAs(Form("%s/evt%04lld.png", outDir.Data(), ev));
      delete cd;
      printf("  evt %lld : %d hits\n", ev, n);
   }

   cs->SaveAs(Form("%s/contact_sheet_%lld_%lld.png", outDir.Data(), evStart, evEnd));
   printf("Done. Output in %s/  (contact sheet + %lld individual PNGs)\n", outDir.Data(), nEv);
}
