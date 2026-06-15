/// @file pid_draw_gate_d2.C
/// @brief Display the D2-target Spyral PID plane and let YOU draw a proton gate by
///        hand, then save it as an AtParticleID json that ex_dp_a1975.C can load.
///
/// D2 difference vs the proton-target pid_draw_gate_a1975.C: the D2 reco has NO
/// separate <run>_pid.root files — the faithful Spyral PID (sqrtdEdx, brho) lives in
/// the AtPIDEvent branch INSIDE the <run>_genfitter_p.root files (GetSpyral()). This
/// macro fills the plane from there, using ALL valid spyral tracks (UNGATED), so you
/// see the full distribution (the auto-derived proton_band_d2.json hugged the wrong
/// ridge — draw the real proton band here).
///
/// Run INTERACTIVELY (needs a display — use `root -l`, NOT `-b`):
///
///   cd .../a1975/D2_UKF
///   root -l 'pid/pid_draw_gate_d2.C("proton_band_d2","pid/proton_band_d2.json")'
///
/// Then on the canvas:
///   • LEFT-CLICK to drop each vertex of the polygon around the proton band
///   • DOUBLE-CLICK to close the polygon
/// The gate is written to <outFile> immediately (xaxis=sqrtdEdx, yaxis=brho, Z=1,A=1
/// — the exact schema ex_dp_a1975.C's AtParticleID::LoadJSON reads), and a PNG
/// snapshot is saved to pid/plots/.
///
/// The filled plane is cached to pid/pid_plane_cache_d2.root, so re-running to redraw
/// is instant. Pass rebuild=true to refill from the <run>_genfitter_p.root files.
///
/// After drawing, re-make the spectrum with the new gate:
///   root -b -q 'ex_dp_a1975.C("run_0016,run_0017,run_0018,run_0019",
///               "/mnt/f/a1975/reco_d2/","","pid/proton_band_d2.json")'

void pid_draw_gate_d2(TString gateName = "proton_band_d2", TString outFile = "pid/proton_band_d2.json",
                      TString runsCSV = "run_0016,run_0017,run_0018,run_0019,run_0020,run_0021,run_0022,run_0023,"
                                        "run_0026,run_0027",
                      TString inDir = "/mnt/f/a1975/reco_d2/", double brMax = 2.0, int Z = 1, int A = 1,
                      bool rebuild = false)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const TString cacheFile = "pid/pid_plane_cache_d2.root";
   TH2F *hS = nullptr;

   if (!rebuild && !gSystem->AccessPathName(cacheFile)) {
      TFile *fc = TFile::Open(cacheFile);
      hS = (TH2F *)fc->Get("hS");
      if (hS)
         hS->SetDirectory(nullptr);
      fc->Close();
      printf("loaded cached PID plane from %s (pass rebuild=true to refill)\n", cacheFile.Data());
   }

   if (!hS) {
      hS = new TH2F("hS", "D2 Spyral PID  (draw your proton gate);#sqrt{dEdx};B#rho [T m]", 300, 0, 50, 400, 0, brMax);
      TObjArray *runs = runsCSV.Tokenize(",");
      long n = 0;
      for (int ri = 0; ri < runs->GetEntries(); ++ri) {
         TString run = ((TObjString *)runs->At(ri))->GetString();
         TString fn = inDir + run + "_genfitter_p.root";
         if (gSystem->AccessPathName(fn)) {
            printf("skip (missing) %s\n", fn.Data());
            continue;
         }
         TFile *f = TFile::Open(fn);
         TTree *t = (TTree *)f->Get("cbmsim");
         TClonesArray *pe = nullptr;
         t->SetBranchAddress("AtPIDEvent", &pe);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (pe->GetEntriesFast() == 0)
               continue;
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (!ev)
               continue;
            for (auto &r : ev->GetSpyral())
               if (r.valid) {
                  hS->Fill(r.sqrtdEdx, r.brho);
                  ++n;
               }
         }
         f->Close();
         printf("processed %s\n", run.Data());
      }
      printf("filled plane with %ld spyral-valid tracks (UNGATED)\n", n);
      TFile *fc = new TFile(cacheFile, "RECREATE");
      hS->Write();
      fc->Close();
      printf("cached plane -> %s\n", cacheFile.Data());
   }

   if (gROOT->IsBatch()) {
      printf("batch mode: cache built, skipping interactive draw. Re-run with `root -l` to draw.\n");
      return;
   }

   TCanvas *c = new TCanvas("cgate", "draw D2 proton PID gate", 950, 780);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   hS->Draw("colz");

   printf("\n\033[1;33m============== DRAW THE PROTON GATE ==============\033[0m\n");
   printf("  • LEFT-CLICK to drop each polygon vertex around the proton band\n");
   printf("  • DOUBLE-CLICK to close the polygon\n");
   printf("  (drawing on the #sqrt{dEdx} vs B#rho plane)\n");
   printf("\033[1;33m=================================================\033[0m\n\n");

   // blocks until you finish drawing a graphical cut named "CUTG"
   TCutG *cut = (TCutG *)gPad->WaitPrimitive("CUTG", "CutG");
   if (!cut || cut->GetN() < 3) {
      printf("\033[1;31mNo gate drawn (need >=3 vertices). Nothing saved.\033[0m\n");
      return;
   }

   int nv = cut->GetN();
   // ROOT closes the TCutG by repeating the first point as the last — drop the dup.
   double x0, y0, xl, yl;
   cut->GetPoint(0, x0, y0);
   cut->GetPoint(nv - 1, xl, yl);
   int nWrite = (nv > 1 && xl == x0 && yl == y0) ? nv - 1 : nv;

   FILE *jf = fopen(outFile.Data(), "w");
   if (!jf) {
      printf("\033[1;31mERROR: could not open %s for writing.\033[0m\n", outFile.Data());
      return;
   }
   fprintf(jf, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdEdx\",\n    \"yaxis\": \"brho\",\n    \"vertices\": [\n",
           gateName.Data());
   for (int i = 0; i < nWrite; ++i) {
      double x, y;
      cut->GetPoint(i, x, y);
      fprintf(jf, "        [%.4f, %.4f]%s\n", x, y, i < nWrite - 1 ? "," : "");
   }
   fprintf(jf, "    ],\n    \"Z\": %d,\n    \"A\": %d\n}\n", Z, A);
   fclose(jf);

   cut->SetLineColor(kRed);
   cut->SetLineWidth(3);
   cut->Draw("same");
   TString png = "pid/plots/drawn_" + gateName + ".png";
   c->SaveAs(png);

   printf("\n\033[1;32mSaved gate '%s' (%d vertices, Z=%d A=%d) -> %s\033[0m\n", gateName.Data(), nWrite, Z, A,
          outFile.Data());
   printf("snapshot -> %s\n", png.Data());
   printf("Re-make the spectrum with it:\n");
   printf("  root -b -q 'ex_dp_a1975.C(\"%s\",\"%s\",\"\",\"%s\")'\n", runsCSV.Data(), inDir.Data(), outFile.Data());
}
