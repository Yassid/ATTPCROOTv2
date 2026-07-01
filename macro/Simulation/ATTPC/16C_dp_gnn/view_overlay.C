/// @file view_overlay.C
/// @brief Event-by-event scanner for the 16C(d,p) GNN overlay training set.
///        Draws each labeled point cloud (blue=proton, red=17C, gray=beam/noise) in
///        XY (pad plane) + ZY (drift side) pads, with on-canvas navigation buttons.
///
/// Build the input once:  ~/Spyral/venv/bin/python convert_overlay_root.py
/// Run (interactive):     source build/config.sh ; root -l 'view_overlay.C(0)'
///   buttons: |< First / < Prev / Next > / Last >| / "Next spiral >"
///   or from the prompt: vnext() vprev() vfirst() vlast() vnext_spiral() view_overlay(N)

#include <set>
#include <vector>

namespace {
   TFile *gF = nullptr;
   TTree *gT = nullptr;
   std::set<int> gEvents;
   int gEv = -1;
   TCanvas *gC = nullptr;
   TPad *gMain = nullptr;
   const char *OV = "data/overlay.root";
   // label -> (color, name)
   const int kCol[3] = {kAzure + 1, kRed + 1, kGray + 1};
   const char *kName[3] = {"proton", "17C", "beam/noise"};
}

void view_overlay_init() {
   if (gT) return;
   gF = TFile::Open(OV);
   if (!gF || gF->IsZombie()) { printf("\033[1;31mcannot open %s — run convert_overlay_root.py\033[0m\n", OV); return; }
   gT = (TTree *)gF->Get("ov");
   int ev; gT->SetBranchAddress("event", &ev);
   for (Long64_t i = 0; i < gT->GetEntries(); ++i) { gT->GetEntry(i); gEvents.insert(ev); }
   gT->ResetBranchAddresses();
   printf("view_overlay: %zu events loaded\n", gEvents.size());
}

// draw one (vy:vx style) projection for a given event, colored per label, on current pad
static void drawProj(int ev, const char *vexp, const char *title, const char *ax, const char *ay) {
   auto *mg = new TMultiGraph();
   mg->SetTitle(Form("%s;%s;%s", title, ax, ay));
   for (int L = 2; L >= 0; --L) {                       // gray first, signal on top
      Long64_t n = gT->Draw(vexp, Form("event==%d && label==%d", ev, L), "goff");
      if (n <= 0) continue;
      auto *g = new TGraph((int)n, gT->GetV2(), gT->GetV1());
      g->SetMarkerColor(kCol[L]); g->SetMarkerStyle(L == 2 ? 1 : 20);
      g->SetMarkerSize(L == 2 ? 0.4 : 0.5); g->SetTitle(kName[L]);
      mg->Add(g, "P");
   }
   mg->Draw("AP");
}

static void buildUI() {
   if (gC) return;
   gC = new TCanvas("cOV", "GNN overlay scanner", 1100, 920);
   auto *bar = new TPad("bar", "", 0, 0.94, 1, 1.0); bar->Draw();
   gMain = new TPad("main", "", 0, 0, 1, 0.94); gMain->Draw(); gMain->Divide(1, 2);
   bar->cd();
   struct B { const char *lbl, *cmd; double x0, x1; };
   B btns[] = {{"|< First", "vfirst();", 0.02, 0.13}, {"< Prev", "vprev();", 0.14, 0.25},
               {"Next >", "vnext();", 0.26, 0.37}, {"Last >|", "vlast();", 0.38, 0.49},
               {"Next spiral >", "vnext_spiral();", 0.66, 0.98}};
   for (auto &b : btns) {
      auto *tb = new TButton(b.lbl, b.cmd, b.x0, 0.1, b.x1, 0.9);
      tb->SetFillColor(kAzure - 9); tb->SetTextSize(0.40); tb->Draw();
   }
   gC->cd();
}

void view_overlay(int ev = 0) {
   view_overlay_init();
   if (!gT) return;
   buildUI();
   gEv = ev;
   gStyle->SetOptStat(0);
   int nP = gT->Draw("1", Form("event==%d && label==0", ev), "goff");
   int n17 = gT->Draw("1", Form("event==%d && label==1", ev), "goff");
   int nBg = gT->Draw("1", Form("event==%d && label==2", ev), "goff");
   gMain->cd(1); gPad->Clear();
   drawProj(ev, "y:x", Form("overlay ev %d  (proton %d / 17C %d / beam-noise %d)", ev, nP, n17, nBg),
            "x [mm]", "y [mm]");
   gMain->cd(2); gPad->Clear();
   drawProj(ev, "y:z", "side view", "z [mm] drift", "y [mm]");
   gC->Update();
   printf("event %d : proton %d  17C %d  beam/noise %d\n", ev, nP, n17, nBg);
}

void vnext()  { auto it = gEvents.upper_bound(gEv); if (it != gEvents.end()) view_overlay(*it); else printf("no next\n"); }
void vprev()  { auto it = gEvents.lower_bound(gEv); if (it != gEvents.begin()) { --it; view_overlay(*it); } else printf("no prev\n"); }
void vfirst() { view_overlay_init(); if (!gEvents.empty()) view_overlay(*gEvents.begin()); }
void vlast()  { view_overlay_init(); if (!gEvents.empty()) view_overlay(*gEvents.rbegin()); }

// jump to the next event whose proton track has > 150 hits (a spiral / long track)
void vnext_spiral() {
   view_overlay_init();
   for (auto it = gEvents.upper_bound(gEv); it != gEvents.end(); ++it) {
      if (gT->Draw("1", Form("event==%d && label==0", *it), "goff") > 150) { view_overlay(*it); return; }
   }
   printf("no further spiral events\n");
}
