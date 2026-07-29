/// @file draw_icgate_C15.C
/// @brief INTERACTIVE ion-chamber (beam) gate selection for a2091 15C.
///
/// Shows the combined IC spectrum from the per-run _ic.root summaries and lets you pick the
/// window with two mouse clicks; writes a small JSON that the analysis reads.
///
/// USAGE (interactively, NOT with -b):
///   root -l 'pid/draw_icgate_C15.C("ic_15C")'                 // pick with the mouse
///   root -l 'pid/draw_icgate_C15.C("ic_15C",1000,1270)'       // or set the window directly
///   root -l 'pid/draw_icgate_C15.C("ic_peak2",1280,1500)'     // the second peak
///
/// HOW TO PICK: when the canvas appears, LEFT-CLICK once at the LOW edge of the window and
/// once at the HIGH edge. Only the x positions are used, so the click height is irrelevant.
/// The macro prints the counts inside the window and saves <name>.json in this folder.
///
/// >>> WHY THIS EXISTS: the inherited window 950-1350 is WRONG -- it cuts straight through the
///     1372 ADC peak, mixing two beam species. The combined spectrum has peaks at 1142 (dominant)
///     and 1372 ADC with the valley at ~1270-1290. WHICH ONE IS 15C IS NOT ESTABLISHED: for 12Be
///     the beam of interest was NOT the dominant peak, so gate each peak and check which yields a
///     clean 15C(p,p) elastic peak at Ex = 0.
///
/// Args: name         = gate name -> <name>.json
///       presetLo/Hi  = set the window non-interactively (both > 0 skips the mouse step)
///       singlePulse  = use only npulse==1 events (pile-up rejected), as gate_events_C15.C does
///       lo,hi        = x-axis display range

void draw_icgate_C15(TString name = "ic_15C", double presetLo = -1, double presetHi = -1,
                     TString icDir = "/home/yassid/a2091_C15_ic/", bool singlePulse = true, double lo = 700,
                     double hi = 1800)
{
   gStyle->SetOptStat(0);
   TString dir = gSystem->DirName(__FILE__);
   TString cache = dir + "/plots/ic_combined.root";

   // ---- combined spectrum, cached so repeat runs are instant ----
   TH1F *hAll = nullptr, *hOne = nullptr;
   TFile *fc = TFile::Open(cache, "READ");
   if (fc && !fc->IsZombie() && fc->Get("hAll") && fc->Get("hOne")) {
      hAll = (TH1F *)((TH1F *)fc->Get("hAll"))->Clone("hAllc");
      hOne = (TH1F *)((TH1F *)fc->Get("hOne"))->Clone("hOnec");
      hAll->SetDirectory(nullptr); hOne->SetDirectory(nullptr);
      fc->Close();
      printf("loaded cached combined IC spectrum\n");
   } else {
      if (fc) fc->Close();
      hAll = new TH1F("hAll", "", 520, 0, 2600);
      hOne = new TH1F("hOne", "", 520, 0, 2600);
      void *dp = gSystem->OpenDirectory(icDir);
      if (!dp) { printf("\033[1;31mERROR: cannot open %s\033[0m\n", icDir.Data()); return; }
      std::vector<TString> runs; const char *e;
      while ((e = gSystem->GetDirEntry(dp))) { TString f(e); if (f.EndsWith("_ic.root")) runs.push_back(f); }
      gSystem->FreeDirectory(dp);
      printf("building combined IC spectrum from %zu runs ...\n", runs.size());
      for (auto &rf : runs) {
         TFile *f = TFile::Open(icDir + rf);
         TTree *t = f ? (TTree *)f->Get("ic") : nullptr;
         if (!t) { if (f) f->Close(); continue; }
         Int_t entry, npulse; Float_t icmax;
         t->SetBranchAddress("entry", &entry); t->SetBranchAddress("icmax", &icmax);
         t->SetBranchAddress("npulse", &npulse);
         for (Long64_t k = 0; k < t->GetEntries(); k++) {
            t->GetEntry(k);
            if (icmax < 0) continue;
            hAll->Fill(icmax);
            if (npulse == 1) hOne->Fill(icmax);
         }
         f->Close();
      }
      TFile *fo = new TFile(cache, "RECREATE");
      hAll->Write("hAll"); hOne->Write("hOne"); fo->Close();
      printf("cached -> %s\n", cache.Data());
   }

   TH1F *h = singlePulse ? hOne : hAll;
   h->SetTitle(TString::Format("a2091 15C combined IC%s  --  click the LOW then HIGH edge;IC max [ADC];events",
                               singlePulse ? " (single pulse)" : ""));
   h->GetXaxis()->SetRangeUser(lo, hi);
   h->SetLineWidth(2);

   double gLo = presetLo, gHi = presetHi;

   if (!(presetLo > 0 && presetHi > 0)) {
      TCanvas *c = new TCanvas("cIC", "select IC gate", 1100, 750);
      h->Draw("hist");
      c->Update();
      printf("\n\033[1;33m======================================================================\033[0m\n");
      printf("\033[1;32m Select the %s window:\033[0m\n", name.Data());
      printf("   * LEFT-CLICK at the LOW edge, then LEFT-CLICK at the HIGH edge\n");
      printf("   * only the x position matters; peaks are at 1142 and 1372 ADC,\n");
      printf("     the valley between them is near 1270-1290\n");
      printf("\033[1;33m======================================================================\033[0m\n\n");

      double xs[2] = {-1, -1};
      for (int i = 0; i < 2; i++) {
         auto *m = (TMarker *)gPad->WaitPrimitive("TMarker", "Marker");
         if (!m) break;
         xs[i] = m->GetX();
         printf("   edge %d at %.0f ADC\n", i + 1, xs[i]);
         // CRITICAL: take the marker back out of the pad. WaitPrimitive returns the first
         // primitive in the pad named "TMarker", so if this one is left in place the SECOND
         // call finds it again and returns immediately with the same x -- giving gLo == gHi
         // and "Invalid window" every time. Not deleted, just unlisted (the pad owns it).
         gPad->GetListOfPrimitives()->Remove(m);
         // replace it with a permanent visual marker of the chosen edge
         auto *le = new TLine(xs[i], gPad->GetUymin(), xs[i], gPad->GetUymax());
         le->SetLineColor(kRed);
         le->SetLineWidth(2);
         le->Draw();
         gPad->Modified();
         gPad->Update();
      }
      gROOT->SetEditorMode(); // leave the canvas in normal mode, not stuck placing markers
      if (xs[0] < 0 || xs[1] < 0) {
         // fallback: whatever range the axis is zoomed to
         int b1 = h->GetXaxis()->GetFirst(), b2 = h->GetXaxis()->GetLast();
         gLo = h->GetXaxis()->GetBinLowEdge(b1);
         gHi = h->GetXaxis()->GetBinUpEdge(b2);
         printf("\033[1;33mNo two markers placed -- falling back to the current zoom: %.0f - %.0f\033[0m\n", gLo, gHi);
      } else {
         gLo = std::min(xs[0], xs[1]);
         gHi = std::max(xs[0], xs[1]);
      }
   }
   if (!(gLo > 0 && gHi > gLo)) {
      printf("\033[1;31mInvalid window: got low=%.1f high=%.1f -- nothing saved.\033[0m\n", gLo, gHi);
      if (gHi > 0 && std::fabs(gHi - gLo) < 1e-6)
         printf("  The two edges are identical. Click ONCE at the low edge, then ONCE at the high\n"
                "  edge (two separate single clicks, not a double-click).\n");
      printf("  You can always set the window directly instead:\n"
             "    root -l 'pid/draw_icgate_C15.C(\"%s\",1000,1270)'\n",
             name.Data());
      return;
   }

   // ---- report ----
   double inW = h->Integral(h->FindBin(gLo), h->FindBin(gHi));
   // NB: bin range explicit. TH1::Integral() with no args honours SetRangeUser, so the
   // display zoom above would silently become the denominator (it reported 79.7% instead
   // of the true 69.0% for the 1000-1270 window).
   double tot = h->Integral(1, h->GetNbinsX());
   printf("\n\033[1;32m=== %s : IC in [%.0f, %.0f] ===\033[0m\n", name.Data(), gLo, gHi);
   printf("  events in window : %.0f  (%.1f%% of %s)\n", inW, 100.0 * inW / tot,
          singlePulse ? "single-pulse" : "all-IC");
   // Reference is the tallest bin INSIDE the window, not the global maximum: a window placed
   // on the secondary 1372 peak must be judged against that peak's own height, otherwise the
   // global 1142 peak makes any edge look negligible.
   int bLo = h->FindBin(gLo), bHi = h->FindBin(gHi);
   int bpk = bLo; double pk = -1;
   for (int b = bLo; b <= bHi; b++) if (h->GetBinContent(b) > pk) { pk = h->GetBinContent(b); bpk = b; }
   printf("  peak in window   : %.0f ADC (height %.0f)\n", h->GetBinCenter(bpk), pk);
   double edgeLo = h->GetBinContent(bLo), edgeHi = h->GetBinContent(bHi);
   printf("  edge heights     : %.0f at low, %.0f at high  (%.1f%% / %.1f%% of the in-window peak)\n", edgeLo,
          edgeHi, pk > 0 ? 100 * edgeLo / pk : 0, pk > 0 ? 100 * edgeHi / pk : 0);
   if (edgeLo > 0.2 * pk || edgeHi > 0.2 * pk)
      printf("  \033[1;33mWARNING: an edge sits above 20%% of the peak -- the window is cutting\n"
             "           through a structure rather than sitting in a valley.\033[0m\n");
   else
      printf("  both edges are low -> the window sits in the valleys, good.\n");

   // ---- save ----
   TString out = dir + "/" + name + ".json";
   std::ofstream j(out.Data());
   j << "{\n";
   j << "  \"name\": \"" << name << "\",\n";
   j << "  \"observable\": \"icmax\",\n";
   j << "  \"comment\": \"max of FRIB generic trace[0] over time buckets 1050-1250\",\n";
   j << "  \"singlePulse\": " << (singlePulse ? "true" : "false") << ",\n";
   j << "  \"icLo\": " << gLo << ",\n";
   j << "  \"icHi\": " << gHi << "\n";
   j << "}\n";
   j.close();
   printf("\n\033[1;32mSaved -> %s\033[0m\n", out.Data());
   printf("  use with: gate_events_C15.C(..., icLo=%.0f, icHi=%.0f)\n", gLo, gHi);

   // draw the chosen window on the spectrum for the record
   TCanvas *c2 = new TCanvas("cICsel", "IC gate", 1100, 750);
   h->Draw("hist");
   c2->Update();
   double ymax = gPad->GetUymax();
   auto *l1 = new TLine(gLo, 0, gLo, ymax); l1->SetLineColor(kRed); l1->SetLineWidth(2); l1->Draw();
   auto *l2 = new TLine(gHi, 0, gHi, ymax); l2->SetLineColor(kRed); l2->SetLineWidth(2); l2->Draw();
   c2->SaveAs(dir + "/plots/icgate_" + name + ".png");
   printf("  plot -> %s/plots/icgate_%s.png\n", dir.Data(), name.Data());
}
