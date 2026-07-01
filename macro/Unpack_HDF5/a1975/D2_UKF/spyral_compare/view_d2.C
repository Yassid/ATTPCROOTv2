/// @file view_d2.C
/// @brief Side-by-side point-cloud viewer: Spyral (left) vs ATTPCROOT-PRA (right),
///        for 16C(d,p) run_0016. Like the previous a1975 viewers — runs in ROOT.
///
/// Spyral clouds come from spyral_clouds.root (TTree "spyral", built by
/// scripts/spyral_to_root.py); ATTPCROOT comes straight from run_0016_reco.root
/// (AtEventCorrected hits + AtPatternEvent tracks). Gray = full point cloud each
/// framework saw; colored = clusters (Spyral HDBSCAN label) / tracks (PRA track id).
///
/// Layout per event: 2 columns (Spyral | ATTPCROOT) x 2 rows (XY pad-plane / ZY side).
///
/// Usage (interactive):
///    source build/config.sh
///    root -l 'view_d2.C(15)'      // show event 15
///    root> vnext()                 // next event present in Spyral tree
///    root> vprev()
///    root> view_d2(11090)          // jump to any event
///
/// NOTE: z (drift) is REFLECTED between the two frameworks (handedness) — compare
/// shapes, not absolute z. The ZY panels are each in native coordinates.

#include <set>
#include <vector>

namespace {
   TFile *gSpy = nullptr, *gReco = nullptr;
   TTree *gSpyT = nullptr, *gRecoT = nullptr;
   std::set<int> gSpyEvents;
   TClonesArray *gCorr = nullptr, *gPat = nullptr;
   int gEv = -1;
   TCanvas *gC = nullptr;
   TPad *gMain = nullptr;
   const int kCol[8] = {kRed+1,kBlue+1,kGreen+2,kMagenta+1,kOrange+7,kCyan+2,kViolet,kSpring-6};

   const char* SPY  = "/home/yassid/spyral_d2/spyral_clouds.root";
   const char* RECO = "/mnt/f/a1975/reco_d2/run_0016_reco.root";
}

void view_d2_init() {
   if (gSpyT) return;
   gSpy = TFile::Open(SPY);
   gSpyT = (TTree*)gSpy->Get("spyral");
   // index which events are present (flat tree: one row per point)
   int ev; gSpyT->SetBranchAddress("event", &ev);
   for (Long64_t i = 0; i < gSpyT->GetEntries(); ++i) { gSpyT->GetEntry(i); gSpyEvents.insert(ev); }
   gSpyT->ResetBranchAddresses();
   gReco = TFile::Open(RECO);
   gRecoT = (TTree*)gReco->Get("cbmsim");
   gRecoT->SetBranchAddress("AtEventCorrected", &gCorr);
   gRecoT->SetBranchAddress("AtPatternEvent",   &gPat);
   printf("view_d2: %zu Spyral events, %lld ATTPCROOT entries loaded\n",
          gSpyEvents.size(), gRecoT->GetEntries());
}

// fill a TGraph from a TTree::Draw selection ("vy:vx"); returns nullptr if empty
static TGraph* graphFromDraw(TTree* t, const char* vexp, const char* cut, int color, double msize, int mstyle) {
   Long64_t n = t->Draw(vexp, cut, "goff");
   if (n <= 0) return nullptr;
   auto* g = new TGraph((int)n, t->GetV2(), t->GetV1());   // vexp = "y:x" -> V1=y, V2=x
   g->SetMarkerColor(color); g->SetMarkerSize(msize); g->SetMarkerStyle(mstyle);
   return g;
}

// ---- Spyral pad: background cloud (gray) + one colored graph per HDBSCAN cluster ----
static int drawSpyralPad(int ev, const char* title, const char* vexp, const char* ax, const char* ay) {
   auto* mg = new TMultiGraph(); mg->SetTitle(Form("%s;%s;%s", title, ax, ay));
   if (auto* g = graphFromDraw(gSpyT, vexp, Form("event==%d && kind==1", ev), kGray+1, 0.3, 1)) mg->Add(g, "P");
   // labels present
   Long64_t nl = gSpyT->Draw("lab", Form("event==%d && kind==0", ev), "goff");
   std::set<int> labs; for (Long64_t i = 0; i < nl; ++i) labs.insert((int)gSpyT->GetV1()[i]);
   int k = 0;
   for (int L : labs)
      if (auto* g = graphFromDraw(gSpyT, vexp, Form("event==%d && kind==0 && lab==%d", ev, L), kCol[k++ % 8], 0.5, 20))
         mg->Add(g, "P");
   mg->Draw("AP");
   return labs.size();
}

// ---- ATTPCROOT pad: corrected hits (gray) + one colored graph per PRA track ----
static int drawAtPad(const char* title, bool side, const char* ax, const char* ay, std::vector<double>& thetaOut) {
   auto* mg = new TMultiGraph(); mg->SetTitle(Form("%s;%s;%s", title, ax, ay));
   if (gCorr && gCorr->GetEntriesFast() > 0) {
      AtEvent* e = (AtEvent*)gCorr->At(0);
      int nh = e->GetNumHits();
      std::vector<double> hx, hy;
      for (int i = 0; i < nh; ++i) { auto p = e->GetHit(i).GetPosition();
         hx.push_back(side ? p.Z() : p.X()); hy.push_back(p.Y()); }
      if (nh) { auto* g = new TGraph(nh, hx.data(), hy.data());
         g->SetMarkerColor(kGray+1); g->SetMarkerStyle(1); mg->Add(g, "P"); }
   }
   int ntrk = 0;
   if (gPat && gPat->GetEntriesFast() > 0) {
      AtPatternEvent* pe = (AtPatternEvent*)gPat->At(0);
      int k = 0;
      for (auto& tr : pe->GetTrackCand()) {
         if (!side) thetaOut.push_back(tr.GetGeoTheta()*TMath::RadToDeg());
         std::vector<double> tx, ty;
         for (auto& hp : tr.GetHitArray()) { auto p = hp->GetPosition();
            tx.push_back(side ? p.Z() : p.X()); ty.push_back(p.Y()); }
         if (!tx.empty()) { auto* g = new TGraph(tx.size(), tx.data(), ty.data());
            g->SetMarkerColor(kCol[k % 8]); g->SetMarkerSize(0.5); g->SetMarkerStyle(20); mg->Add(g, "P"); }
         ++k; ++ntrk;
      }
   }
   mg->Draw("AP");
   return ntrk;
}

// build the canvas once: a thin button strip on top + a 2x2 main pad below
static void buildUI() {
   if (gC) return;
   gC = new TCanvas("cD2", "Spyral vs ATTPCROOT", 1300, 940);
   auto* bar = new TPad("bar", "", 0.0, 0.94, 1.0, 1.0); bar->Draw();
   gMain = new TPad("main", "", 0.0, 0.0, 1.0, 0.94); gMain->Draw(); gMain->Divide(2, 2);
   bar->cd();
   struct B { const char* lbl; const char* cmd; double x0; double x1; };
   B btns[] = {{"|< First","vfirst();",0.02,0.12}, {"< Prev","vprev();",0.13,0.23},
               {"Next >","vnext();",0.24,0.34}, {"Last >|","vlast();",0.35,0.45},
               {"Next only-PRA-back >","vnext_onlypra();",0.60,0.98}};
   for (auto& b : btns) {
      auto* tb = new TButton(b.lbl, b.cmd, b.x0, 0.1, b.x1, 0.9);
      tb->SetFillColor(kAzure - 9); tb->SetTextSize(0.40); tb->Draw();
   }
   gC->cd();
}

void view_d2(int ev = 15) {
   view_d2_init();
   buildUI();
   gEv = ev;
   gStyle->SetOptStat(0);
   if (ev >= 0 && ev < gRecoT->GetEntries()) gRecoT->GetEntry(ev);
   std::vector<double> theta;

   gMain->cd(1); gPad->Clear(); int nSp = drawSpyralPad(ev, Form("Spyral ev %d", ev), "y:x", "x [mm]", "y [mm]");
   gMain->cd(2); gPad->Clear(); int nAt = drawAtPad(Form("ATTPCROOT ev %d", ev), false, "x [mm]", "y [mm]", theta);
   gMain->cd(3); gPad->Clear(); drawSpyralPad(ev, "Spyral  side view", "y:z", "z [mm] drift", "y [mm]");
   gMain->cd(4); gPad->Clear(); drawAtPad("ATTPCROOT  side view", true, "z [mm] drift", "y [mm]", theta);
   gMain->cd(1)->SetTitle(Form("Spyral ev %d  (%d clusters)", ev, nSp));
   gC->Update();

   TString th; for (double t : theta) th += Form("%.0f ", t);
   printf("event %d : Spyral %d clusters | ATTPCROOT %d tracks  theta=[ %s]\n", ev, nSp, nAt, th.Data());
}

// ---- navigation (callable from buttons or the prompt) ----
void vnext()  { auto it = gSpyEvents.upper_bound(gEv); if (it != gSpyEvents.end()) view_d2(*it); else printf("no next\n"); }
void vprev()  { auto it = gSpyEvents.lower_bound(gEv); if (it != gSpyEvents.begin()) { --it; view_d2(*it); } else printf("no prev\n"); }
void vfirst() { view_d2_init(); if (!gSpyEvents.empty()) view_d2(*gSpyEvents.begin()); }
void vlast()  { view_d2_init(); if (!gSpyEvents.empty()) view_d2(*gSpyEvents.rbegin()); }

// jump to the next "only-PRA-backward" event (ATTPCROOT has a backward 90<theta<165
// track but Spyral kept nothing backward) — the key disagreement class to inspect.
void vnext_onlypra() {
   view_d2_init();
   for (auto it = gSpyEvents.upper_bound(gEv); it != gSpyEvents.end(); ++it) {
      int ev = *it;
      if (ev < 0 || ev >= gRecoT->GetEntries()) continue;
      gRecoT->GetEntry(ev);
      bool praBack = false;
      if (gPat && gPat->GetEntriesFast() > 0)
         for (auto& tr : ((AtPatternEvent*)gPat->At(0))->GetTrackCand()) {
            double th = tr.GetGeoTheta()*TMath::RadToDeg();
            if (th > 90 && th < 165 && tr.GetHitArray().size() >= 20) { praBack = true; break; }
         }
      if (!praBack) continue;
      // require Spyral has NO backward cluster: approximate via no cluster point with theta>90
      Long64_t nb = gSpyT->Draw("1",
         Form("event==%d && kind==0 && (sqrt(x*x+y*y) > 0) && (atan2(sqrt(x*x+y*y),z)*57.2958 > 90)", ev), "goff");
      if (nb == 0) { view_d2(ev); return; }
   }
   printf("no further only-PRA-backward events\n");
}
