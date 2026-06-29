// Interactive Spyral-vs-GNN event viewer (ROOT macro, 2D TGX11 GUI -- works under WSLg;
// ROOT Eve/3D does not). Next/Prev buttons + event jump + z-y/x-y toggle. Two pads:
// Spyral clustering (left) vs GNN clustering (right) on the SAME points, colored by cluster.
//
// Build the CSVs first (once):  viewer/view.sh   (or run viewer/spyral_extract.py + gnn_on_points.py)
// Then:   root -l 'viewer/view_events.C'                       // default run_0300, embed_self2
//   or:   root -l 'viewer/view_events.C("data/spyral_clusters.csv","data/gnn_clusters.csv")'
#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TString.h>
#include <TApplication.h>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

struct Pt { double x, y, z; int cl; };
typedef std::map<int, std::vector<Pt>> EvMap;

static void loadCSV(const char *fn, EvMap &m) {
   std::ifstream in(fn);
   if (!in) { printf("ERROR: cannot open %s\n", fn); return; }
   std::string line; std::getline(in, line);                 // header
   while (std::getline(in, line)) {
      std::stringstream ss(line); std::string t;
      std::getline(ss, t, ','); int ev = std::stoi(t);
      std::getline(ss, t, ','); double x = std::stod(t);
      std::getline(ss, t, ','); double y = std::stod(t);
      std::getline(ss, t, ','); double z = std::stod(t);
      std::getline(ss, t, ',');                              // q (unused)
      std::getline(ss, t, ','); int cl = std::stoi(t);
      m[ev].push_back({x, y, z, cl});
   }
}

class EvViewer : public TGMainFrame {
   TRootEmbeddedCanvas *fEC;
   TGNumberEntry *fNum;
   TGLabel *fInfo;
   EvMap fSp, fGn;
   std::vector<int> fIds;
   int fIdx, fProj;                                          // fProj: 0=z-y, 1=x-y
   std::vector<TObject *> fGarb;
public:
   EvViewer(const char *spcsv, const char *gncsv);
   void Next() { if (!fIds.empty()) { fIdx = (fIdx + 1) % fIds.size(); DrawEv(); } }
   void Prev() { if (!fIds.empty()) { fIdx = (fIdx - 1 + fIds.size()) % fIds.size(); DrawEv(); } }
   void Goto() { int id = fNum->GetIntNumber();
                 for (size_t i = 0; i < fIds.size(); ++i) if (fIds[i] == id) { fIdx = i; DrawEv(); return; } }
   void SetZY() { fProj = 0; DrawEv(); }
   void SetXY() { fProj = 1; DrawEv(); }
   void DrawEv();
   void DrawPad(EvMap &m, const char *tag, int evid);
   void CloseWindow() override { gApplication->Terminate(0); }
   ~EvViewer() override { Cleanup(); }
   ClassDef(EvViewer, 0)
};

static int kCol(int c) {                                     // cluster -> ROOT color; noise=gray
   static int p[] = {kRed, kGreen + 2, kBlue, kOrange + 1, kViolet, kCyan + 1,
                     kMagenta, kSpring + 4, kAzure + 1, kPink, kYellow + 1, kTeal};
   return c < 0 ? kGray : p[c % 12];
}

void EvViewer::DrawPad(EvMap &m, const char *tag, int evid) {
   auto it = m.find(evid);
   double a0 = 1e9, a1 = -1e9, y0 = 1e9, y1 = -1e9;
   std::map<int, std::vector<Pt>> byc;
   if (it != m.end())
      for (auto &p : it->second) {
         double a = fProj == 0 ? p.z : p.x;
         a0 = std::min(a0, a); a1 = std::max(a1, a);
         y0 = std::min(y0, p.y); y1 = std::max(y1, p.y);
         byc[p.cl].push_back(p);
      }
   if (a0 > a1) { a0 = 0; a1 = 1; y0 = 0; y1 = 1; }
   double ma = 0.05 * (a1 - a0) + 1, my = 0.05 * (y1 - y0) + 1;
   int ncl = 0; for (auto &kv : byc) if (kv.first >= 0) ncl++;
   TH2F *fr = new TH2F(Form("fr_%s_%d", tag, evid), Form("%s  event %d  (%d clusters)", tag, evid, ncl),
                       10, a0 - ma, a1 + ma, 10, y0 - my, y1 + my);
   fr->SetStats(0); fr->GetXaxis()->SetTitle(fProj == 0 ? "z (mm)" : "x (mm)");
   fr->GetYaxis()->SetTitle("y (mm)"); fr->Draw(); fGarb.push_back(fr);
   for (auto &kv : byc) {
      TGraph *g = new TGraph();
      for (auto &p : kv.second) g->SetPoint(g->GetN(), fProj == 0 ? p.z : p.x, p.y);
      g->SetMarkerStyle(20); g->SetMarkerSize(0.6); g->SetMarkerColor(kCol(kv.first));
      g->Draw("P same"); fGarb.push_back(g);
   }
}

void EvViewer::DrawEv() {
   for (auto o : fGarb) delete o; fGarb.clear();
   int ev = fIds[fIdx];
   TCanvas *c = fEC->GetCanvas();
   c->Clear(); c->Divide(2, 1);
   c->cd(1); DrawPad(fSp, "ATTPCROOT (triplclust)", ev);
   c->cd(2); DrawPad(fGn, "Spyral (HDBSCAN)", ev);
   c->Update();
   fNum->SetIntNumber(ev);
   int nsp = fSp.count(ev) ? fSp[ev].size() : 0;
   fInfo->SetText(Form(" #%d / %d   event %d   (%d hits) ", fIdx + 1, (int)fIds.size(), ev, nsp));
}

EvViewer::EvViewer(const char *spcsv, const char *gncsv)
   : TGMainFrame(gClient->GetRoot(), 1200, 640), fIdx(0), fProj(0) {
   loadCSV(spcsv, fSp); loadCSV(gncsv, fGn);
   for (auto &kv : fSp) fIds.push_back(kv.first);
   printf("loaded %d events (Spyral) / %d (GNN)\n", (int)fSp.size(), (int)fGn.size());

   TGHorizontalFrame *bar = new TGHorizontalFrame(this, 1200, 40);
   TGTextButton *bp = new TGTextButton(bar, " < Prev ");
   TGTextButton *bn = new TGTextButton(bar, " Next > ");
   bp->Connect("Clicked()", "EvViewer", this, "Prev()");
   bn->Connect("Clicked()", "EvViewer", this, "Next()");
   bar->AddFrame(bp, new TGLayoutHints(kLHintsLeft, 4, 2, 4, 4));
   bar->AddFrame(bn, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   bar->AddFrame(new TGLabel(bar, "event:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
   fNum = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger);
   fNum->Connect("ValueSet(Long_t)", "EvViewer", this, "Goto()");
   fNum->GetNumberEntry()->Connect("ReturnPressed()", "EvViewer", this, "Goto()");
   bar->AddFrame(fNum, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   TGTextButton *bzy = new TGTextButton(bar, " z-y ");
   TGTextButton *bxy = new TGTextButton(bar, " x-y ");
   bzy->Connect("Clicked()", "EvViewer", this, "SetZY()");
   bxy->Connect("Clicked()", "EvViewer", this, "SetXY()");
   bar->AddFrame(bzy, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   bar->AddFrame(bxy, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   fInfo = new TGLabel(bar, "  Spyral (left)  vs  GNN (right) ");
   bar->AddFrame(fInfo, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 8, 8, 4, 4));
   AddFrame(bar, new TGLayoutHints(kLHintsExpandX));

   fEC = new TRootEmbeddedCanvas("ec", this, 1200, 590);
   AddFrame(fEC, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   SetWindowName("Spyral vs GNN  -  16C(d,p) run_0016");
   MapSubwindows(); Resize(GetDefaultSize()); MapWindow();
   if (!fIds.empty()) DrawEv();
}

void view_events(const char *spcsv = "data/spyral_clusters.csv",
                 const char *gncsv = "data/gnn_clusters.csv") {
   new EvViewer(spcsv, gncsv);
}
