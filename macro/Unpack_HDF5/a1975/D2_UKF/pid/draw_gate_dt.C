/// @file draw_gate_dt.C
/// @brief Draw the a1975 (d,t) triton gate BY HAND on the dv 1.10424 Spyral PID plane.
///
/// Reads the points cache written by pid/pid_plane_dt.C -- i.e. the gate is drawn on exactly
/// the numbers the analysis will cut on, not on a redrawn approximation of them, and no fit
/// has to exist first.
///
/// THREE CONVENTIONS, each of which is here because ignoring it cost work before:
///
///   1. AUTOSAVE ON POLYGON CLOSE.  The gate is written the moment the polygon closes, not
///      when the window is shut.  A gate drawn and then dismissed by closing the GUI used to
///      be lost outright.
///   2. Z, A AND A UNIQUE NAME are written into the json.  Without them the loader cannot tell
///      one polygon from another when several are open, and the wrong one gets applied.
///   3. The TCutG is UNLINKED from the pad before anything deletes it.  Deleting a TCutG a pad
///      still owns corrupts the heap (realloc(): invalid next size) and has eaten a finished
///      gate.
///
/// RUN IT INTERACTIVELY -- it needs a display, so `root -l`, never `-b`:
///
///   cd .../a1975/D2_UKF
///   root -l 'pid/draw_gate_dt.C()'
///
/// On the canvas: LEFT-CLICK each vertex around the triton band; DOUBLE-CLICK (or click the
/// first vertex again) to close it.  The json is written at that moment and the polygon is
/// redrawn in green as confirmation.  Re-run to redraw from scratch.

#include <fstream>

namespace
{
TCutG *gCut = nullptr;
TString gOut;
TString gName;
int gZ = 1, gA = 3;
TGraph *gLive = nullptr;
std::vector<double> gX, gY;
// A spurious kButton1Down arrives when the canvas is first mapped -- it put a phantom vertex
// at Brho 1.55 on the first run. Ignore anything in the first second, and anything outside
// the axis ranges, which a real click on the plot cannot be.
Long64_t gT0 = 0;
double gXlo = 0, gXhi = 0, gYlo = 0, gYhi = 0;

void SaveGate()
{
   if (gX.size() < 3)
      return;
   std::ofstream o(gOut.Data());
   o << "{\n \"name\": \"" << gName << "\",\n";
   o << " \"xaxis\": \"sqrtdedx\",\n \"yaxis\": \"brho\",\n";
   o << " \"Z\": " << gZ << ",\n \"A\": " << gA << ",\n";
   o << " \"vertices\": [\n";
   for (size_t i = 0; i < gX.size(); ++i)
      o << "  [" << gX[i] << ", " << gY[i] << "]" << (i + 1 < gX.size() ? "," : "") << "\n";
   o << "  [" << gX[0] << ", " << gY[0] << "]\n ]\n}\n"; // closed
   o.close();
   printf("\n\033[1;32m*** gate saved: %s  (%zu vertices, Z=%d A=%d) ***\033[0m\n", gOut.Data(), gX.size(), gZ, gA);
}
} // namespace

/// click handler: accumulate vertices, close on double-click
void dtGateExec()
{
   int ev = gPad->GetEvent();
   if (ev != kButton1Down && ev != kButton1Double)
      return;
   if (gT0 == 0 || (Long64_t)gSystem->Now() - gT0 < 1000)
      return; // canvas still settling
   double x = gPad->AbsPixeltoX(gPad->GetEventX());
   double y = gPad->AbsPixeltoY(gPad->GetEventY());
   x = gPad->PadtoX(x);
   y = gPad->PadtoY(y);

   if (x < gXlo || x > gXhi || y < gYlo || y > gYhi) {
      printf("  (ignored a click outside the axes: %.3f, %.4f)\n", x, y);
      return;
   }

   if (ev == kButton1Double) {
      SaveGate();
      // draw the closed polygon in green as confirmation
      auto *g = new TGraph((int)gX.size() + 1);
      for (size_t i = 0; i < gX.size(); ++i)
         g->SetPoint((int)i, gX[i], gY[i]);
      g->SetPoint((int)gX.size(), gX[0], gY[0]);
      g->SetLineColor(kGreen + 2);
      g->SetLineWidth(3);
      g->Draw("L same");
      gPad->Update();
      return;
   }
   gX.push_back(x);
   gY.push_back(y);
   printf("  vertex %2zu : sqrtdEdx %8.3f   Brho %7.4f\n", gX.size(), x, y);
   if (!gLive) {
      gLive = new TGraph();
      gLive->SetLineColor(kRed);
      gLive->SetMarkerColor(kRed);
      gLive->SetMarkerStyle(20);
      gLive->SetMarkerSize(0.8);
      gLive->SetLineWidth(2);
      gLive->Draw("LP same");
   }
   gLive->SetPoint(gLive->GetN(), x, y);
   gPad->Modified();
   gPad->Update();
}

void draw_gate_dt(TString cache = "pid/pid_plane_dt_dv1104.root", TString outJson = "pid/triton_d2_dv1104.json",
                  TString gateName = "triton_d2_dv1104", int Z = 1, int A = 3, TString overlay = "",
                  double sqrtMax = 60, double brhoMax = 2.0)
{
   gOut = outJson;
   gName = gateName;
   gZ = Z;
   gA = A;
   gX.clear();
   gY.clear();
   gLive = nullptr;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("cannot open %s -- run pid/pid_plane_dt.C first\n", cache.Data());
      return;
   }
   auto *t = (TTree *)f->Get("pts");
   if (!t) {
      printf("no tree 'pts' in %s\n", cache.Data());
      return;
   }
   printf("=== %lld PID points from %s ===\n", t->GetEntries(), cache.Data());

   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   auto *c = new TCanvas("cdt", "a1975 (d,t) PID -- draw the triton gate", 1100, 850);
   c->SetLogz();
   c->SetRightMargin(0.13);
   t->Draw(Form("brho:sqrtdedx>>hdt(300,0,%g,300,0,%g)", sqrtMax, brhoMax), "", "colz");
   auto *h = (TH2F *)gDirectory->Get("hdt");
   h->SetTitle("a1975 ^{16}C(d,t) Spyral PID, dv 1.10424;#sqrt{dE/dx};B#rho [T m]");

   if (overlay.Length()) { // optional: show an old gate for reference, NOT to copy
      std::ifstream in(overlay.Data());
      if (in) {
         std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
         auto p = s.find('[', s.find("vertices"));
         if (p != std::string::npos) {
            std::vector<double> v;
            const char *q = s.c_str() + p, *e = s.c_str() + s.size();
            int depth = 0;
            while (q < e) {
               if (*q == '[')
                  depth++;
               if (*q == ']') {
                  depth--;
                  if (depth <= 0)
                     break;
               }
               char *np = nullptr;
               double d = strtod(q, &np);
               if (np != q) {
                  v.push_back(d);
                  q = np;
               } else
                  ++q;
            }
            auto *go = new TGraph();
            for (size_t i = 0; i + 1 < v.size(); i += 2)
               go->SetPoint(go->GetN(), v[i], v[i + 1]);
            go->SetLineColor(kGray + 2);
            go->SetLineStyle(2);
            go->SetLineWidth(2);
            go->Draw("L same");
            printf("overlaid (dashed grey, REFERENCE ONLY -- it belongs to the old dv): %s\n", overlay.Data());
         }
      }
   }

   gXlo = 0;
   gXhi = sqrtMax;
   gYlo = 0;
   gYhi = brhoMax;
   c->Update();
   gT0 = (Long64_t)gSystem->Now();
   c->AddExec("dtgate", "dtGateExec()");
   printf("\n  LEFT-CLICK each vertex around the triton band.\n");
   printf("  DOUBLE-CLICK to close the polygon -- it is saved to %s at that moment.\n\n", outJson.Data());
}
