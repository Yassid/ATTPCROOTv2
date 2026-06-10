/// @file pid_draw_gate_a1975.C
/// @brief Display the faithful Spyral PID plane and let YOU draw a gate by hand,
///        then save it as an AtParticleID json the analysis macros can load.
///
/// Run INTERACTIVELY (needs a display — use `root -l`, NOT `-b`):
///
///   root -l 'pid/pid_draw_gate_a1975.C("deuteron_band","pid/deuteron_band.json")'
///
/// Then on the canvas:
///   • LEFT-CLICK to drop each vertex of the polygon around the band
///   • DOUBLE-CLICK to close the polygon
/// The gate is written to <outFile> immediately (xaxis=sqrtdedx, yaxis=brho), and
/// a PNG snapshot is saved to pid/plots/.
///
/// The filled plane is cached to pid/pid_plane_cache.root, so re-running to redraw
/// a gate is instant. Pass rebuild=true to refill from the <run>_pid.root files.
///
/// Args: gateName (json "name"), outFile (json path), runsCSV, inDir, brMax (y-axis
/// max — keep 2.0 so the deuteron band above the proton band is visible), rebuild.

void pid_draw_gate_a1975(TString gateName = "deuteron_band", TString outFile = "pid/deuteron_band.json",
                         TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                           "run_0114,run_0115",
                         TString inDir = "/mnt/f/a1975/reco/", double brMax = 2.0, bool rebuild = false)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const TString cacheFile = "pid/pid_plane_cache.root";
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
      hS = new TH2F("hS", "Faithful Spyral PID  (draw your gate);#sqrt{dEdx};B#rho [T m]", 300, 0, 40, 400, 0, brMax);
      TObjArray *runs = runsCSV.Tokenize(",");
      long n = 0;
      for (int ri = 0; ri < runs->GetEntries(); ++ri) {
         TString run = ((TObjString *)runs->At(ri))->GetString();
         TString fn = inDir + run + "_pid.root";
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
      printf("filled plane with %ld spyral-valid tracks\n", n);
      TFile *fc = new TFile(cacheFile, "RECREATE");
      hS->Write();
      fc->Close();
      printf("cached plane -> %s\n", cacheFile.Data());
   }

   if (gROOT->IsBatch()) {
      printf("batch mode: cache built, skipping interactive draw. Re-run with `root -l` to draw.\n");
      return;
   }

   TCanvas *c = new TCanvas("cgate", "draw PID gate", 950, 780);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   hS->Draw("colz");

   printf("\n\033[1;33m================ DRAW THE GATE ================\033[0m\n");
   printf("  • LEFT-CLICK to drop each polygon vertex around the band\n");
   printf("  • DOUBLE-CLICK to close the polygon\n");
   printf("  (drawing on the #sqrt{dEdx} vs B#rho plane)\n");
   printf("\033[1;33m===============================================\033[0m\n\n");

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
   fprintf(jf, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n    \"vertices\": [\n",
           gateName.Data());
   for (int i = 0; i < nWrite; ++i) {
      double x, y;
      cut->GetPoint(i, x, y);
      fprintf(jf, "        [%.3f, %.3f]%s\n", x, y, i < nWrite - 1 ? "," : "");
   }
   fprintf(jf, "    ]\n}\n");
   fclose(jf);

   cut->SetLineColor(kRed);
   cut->SetLineWidth(3);
   cut->Draw("same");
   TString png = "pid/plots/drawn_" + gateName + ".png";
   c->SaveAs(png);

   printf("\n\033[1;32mSaved gate '%s' (%d vertices) -> %s\033[0m\n", gateName.Data(), nWrite, outFile.Data());
   printf("snapshot -> %s\n", png.Data());
   printf("Apply it: root -b -q 'pid/gate_deuteron_a1975.C(\"...\",\"%s\",\"%s\")'\n", inDir.Data(), outFile.Data());
}
