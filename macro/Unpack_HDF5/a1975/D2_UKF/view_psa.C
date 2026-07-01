// Navigate events comparing ATTPCROOT vs Spyral point clouds (ROOT GUI, works under WSLg).
// Left  = ATTPCROOT AtPSAMultiFit aligned to Spyral (thr40, prom20, sep50, Spyral z-calibration).
// Right = Spyral.  Both now in the SAME z frame (no flip).  Prev/Next + event jump + z-y/x-y toggle.
//
// Build the CSVs first (already done):  mfs_cloud.csv (ATTPCROOT), spyral_cloud.csv (Spyral).
// Run:   root -l 'view_psa.C+'      (the + compiles it once with ACLiC)
#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TApplication.h>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

struct Pt { double x, y, z; };
typedef std::map<int, std::vector<Pt>> EvMap;

static void loadCSV(const char *fn, EvMap &m, double zflip)
{
   std::ifstream in(fn);
   if (!in) { printf("ERROR: cannot open %s\n", fn); return; }
   std::string line; std::getline(in, line);                 // header: eid,x,y,z,q
   while (std::getline(in, line)) {
      std::stringstream ss(line); std::string t;
      std::getline(ss, t, ','); int ev = std::stoi(t);
      std::getline(ss, t, ','); double x = std::stod(t);
      std::getline(ss, t, ','); double y = std::stod(t);
      std::getline(ss, t, ','); double z = std::stod(t);
      if (zflip > 0) z = zflip - z;
      m[ev].push_back({x, y, z});
   }
}

class PsaViewer : public TGMainFrame {
   TRootEmbeddedCanvas *fEC;
   TGNumberEntry *fNum;
   TGLabel *fInfo;
   EvMap fAt, fSp;
   std::vector<int> fIds;
   int fIdx, fProj;                                          // fProj: 0=z-y, 1=x-y
   std::vector<TObject *> fGarb;
public:
   PsaViewer(const char *atCsv, const char *spCsv);
   void Next() { if (!fIds.empty()) { fIdx = (fIdx + 1) % fIds.size(); DrawEv(); } }
   void Prev() { if (!fIds.empty()) { fIdx = (fIdx - 1 + fIds.size()) % fIds.size(); DrawEv(); } }
   void Goto() { int id = fNum->GetIntNumber();
                 for (size_t i = 0; i < fIds.size(); ++i) if (fIds[i] == id) { fIdx = i; DrawEv(); return; } }
   void SetZY() { fProj = 0; DrawEv(); }
   void SetXY() { fProj = 1; DrawEv(); }
   void DrawEv();
   void DrawPad(EvMap &m, const char *tag, int evid, int color);
   void CloseWindow() override { gApplication->Terminate(0); }
   ~PsaViewer() override { Cleanup(); }
   ClassDef(PsaViewer, 0)
};

void PsaViewer::DrawPad(EvMap &m, const char *tag, int evid, int color)
{
   auto it = m.find(evid);
   int nh = (it != m.end()) ? it->second.size() : 0;
   double xlo = (fProj == 0) ? 50 : -280, xhi = (fProj == 0) ? 1100 : 280;
   TH2F *fr = new TH2F(Form("fr_%s_%d", tag, evid), Form("%s  event %d  (%d hits)", tag, evid, nh),
                       10, xlo, xhi, 10, -280, 280);
   fr->SetStats(0); fr->GetXaxis()->SetTitle(fProj == 0 ? "z (mm)" : "x (mm)");
   fr->GetYaxis()->SetTitle("y (mm)"); fr->Draw(); fGarb.push_back(fr);
   if (nh > 0) {
      TGraph *g = new TGraph();
      for (auto &p : it->second) g->SetPoint(g->GetN(), fProj == 0 ? p.z : p.x, p.y);
      g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(color);
      g->Draw("P same"); fGarb.push_back(g);
   }
}

void PsaViewer::DrawEv()
{
   for (auto o : fGarb) delete o; fGarb.clear();
   int ev = fIds[fIdx];
   TCanvas *c = fEC->GetCanvas();
   c->Clear(); c->Divide(2, 1);
   c->cd(1); DrawPad(fAt, "ATTPCROOT", ev, kAzure + 1);
   c->cd(2); DrawPad(fSp, "SPYRAL", ev, kOrange + 7);
   c->Update();
   fNum->SetIntNumber(ev);
   int na = fAt.count(ev) ? fAt[ev].size() : 0, ns = fSp.count(ev) ? fSp[ev].size() : 0;
   fInfo->SetText(Form(" #%d / %d   event %d   ATTPCROOT %d  vs  Spyral %d hits ", fIdx + 1, (int)fIds.size(), ev, na, ns));
}

PsaViewer::PsaViewer(const char *atCsv, const char *spCsv)
   : TGMainFrame(gClient->GetRoot(), 1300, 660), fIdx(0), fProj(0)
{
   loadCSV(atCsv, fAt, -1);   // ATTPCROOT already uses Spyral z-calibration now (no flip)
   loadCSV(spCsv, fSp, -1);
   for (auto &kv : fAt) if (fSp.count(kv.first)) fIds.push_back(kv.first); // events in both
   printf("loaded %d ATTPCROOT / %d Spyral events; %d matched\n", (int)fAt.size(), (int)fSp.size(), (int)fIds.size());

   TGHorizontalFrame *bar = new TGHorizontalFrame(this, 1300, 40);
   TGTextButton *bp = new TGTextButton(bar, " < Prev ");
   TGTextButton *bn = new TGTextButton(bar, " Next > ");
   bp->Connect("Clicked()", "PsaViewer", this, "Prev()");
   bn->Connect("Clicked()", "PsaViewer", this, "Next()");
   bar->AddFrame(bp, new TGLayoutHints(kLHintsLeft, 4, 2, 4, 4));
   bar->AddFrame(bn, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   bar->AddFrame(new TGLabel(bar, "event:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
   fNum = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger);
   fNum->Connect("ValueSet(Long_t)", "PsaViewer", this, "Goto()");
   fNum->GetNumberEntry()->Connect("ReturnPressed()", "PsaViewer", this, "Goto()");
   bar->AddFrame(fNum, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   TGTextButton *bzy = new TGTextButton(bar, " z-y ");
   TGTextButton *bxy = new TGTextButton(bar, " x-y ");
   bzy->Connect("Clicked()", "PsaViewer", this, "SetZY()");
   bxy->Connect("Clicked()", "PsaViewer", this, "SetXY()");
   bar->AddFrame(bzy, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   bar->AddFrame(bxy, new TGLayoutHints(kLHintsLeft, 2, 8, 4, 4));
   fInfo = new TGLabel(bar, " left: ATTPCROOT   right: Spyral ");
   bar->AddFrame(fInfo, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 8, 8, 4, 4));
   AddFrame(bar, new TGLayoutHints(kLHintsExpandX));

   fEC = new TRootEmbeddedCanvas("ec", this, 1300, 610);
   AddFrame(fEC, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   SetWindowName("ATTPCROOT vs Spyral PSA  -  run_0300");
   MapSubwindows(); Resize(GetDefaultSize()); MapWindow();
   if (!fIds.empty()) DrawEv();
}

void view_psa(const char *atCsv = "mfs_cloud.csv", const char *spCsv = "spyral_cloud.csv")
{
   new PsaViewer(atCsv, spCsv);
}
