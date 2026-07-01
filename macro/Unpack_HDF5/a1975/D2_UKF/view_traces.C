// Interactive pad-TRACE viewer (ROOT GUI, works under WSLg). Step through events and pads and
// see each pad's raw ADC trace, the threshold line, and where the PSA put hits.
// Run:   root -l 'view_traces.C+'                       // default run_0300_psa_max.root
//   or:  root -l 'view_traces.C+("/path/to/reco.root", 20)'   // file, threshold line
#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TApplication.h>
#include <map>
#include <vector>
#include <algorithm>

#include "AtRawEvent.h"
#include "AtEvent.h"
#include "AtPad.h"
#include "AtHit.h"

class TraceViewer : public TGMainFrame {
   TRootEmbeddedCanvas *fEC;
   TGNumberEntry *fEvNum, *fPadNum;
   TGLabel *fInfo;
   TFile *fFile;
   TTree *fTree;
   TClonesArray *fRaw, *fEv;
   Long64_t fNEv, fEvIdx;
   int fPadIdx;
   double fThr;
   std::vector<int> fPads;                 // pad numbers in current event (sorted)
   std::map<int, std::vector<double>> fHitTB; // PSA hit time buckets per pad
   std::vector<TObject *> fGarb;
public:
   TraceViewer(const char *file, double thr);
   void LoadEvent();
   void NextEv() { fEvIdx = (fEvIdx + 1) % fNEv; LoadEvent(); }
   void PrevEv() { fEvIdx = (fEvIdx - 1 + fNEv) % fNEv; LoadEvent(); }
   void GotoEv() { fEvIdx = std::min((Long64_t)fEvNum->GetIntNumber(), fNEv - 1); LoadEvent(); }
   void NextPad() { if (!fPads.empty()) { fPadIdx = (fPadIdx + 1) % fPads.size(); Draw(); } }
   void PrevPad() { if (!fPads.empty()) { fPadIdx = (fPadIdx - 1 + fPads.size()) % fPads.size(); Draw(); } }
   void GotoPad();
   void Draw();
   void CloseWindow() override { gApplication->Terminate(0); }
   ~TraceViewer() override { Cleanup(); }
   ClassDef(TraceViewer, 0)
};

void TraceViewer::LoadEvent()
{
   fTree->GetEntry(fEvIdx);
   fPads.clear(); fHitTB.clear();
   auto raw = (AtRawEvent *)fRaw->At(0);
   if (raw)
      for (const auto &p : raw->GetPads())
         fPads.push_back(p->GetPadNum());
   std::sort(fPads.begin(), fPads.end());
   auto ev = (AtEvent *)fEv->At(0);
   if (ev)
      for (const auto &h : ev->GetHits())
         fHitTB[h->GetPadNum()].push_back(h->GetTimeStamp());
   fPadIdx = 0;
   fEvNum->SetIntNumber(fEvIdx);
   Draw();
}

void TraceViewer::GotoPad()
{
   int want = fPadNum->GetIntNumber();
   for (size_t i = 0; i < fPads.size(); ++i)
      if (fPads[i] >= want) { fPadIdx = i; Draw(); return; }
}

void TraceViewer::Draw()
{
   for (auto o : fGarb) delete o; fGarb.clear();
   TCanvas *c = fEC->GetCanvas();
   c->cd();
   if (fPads.empty()) { c->Clear(); c->Update(); return; }
   int pad = fPads[fPadIdx];
   auto raw = (AtRawEvent *)fRaw->At(0);
   const AtPad *thePad = nullptr;
   for (const auto &p : raw->GetPads())
      if (p->GetPadNum() == pad) { thePad = p.get(); break; }
   bool hasHit = fHitTB.count(pad);
   TH2F *fr = new TH2F(Form("fr_%lld_%d", fEvIdx, pad), Form("event %lld   pad %d   (%lu/%lu in event)%s",
                       fEvIdx, pad, (unsigned long)(fPadIdx + 1), (unsigned long)fPads.size(),
                       hasHit ? "   <-- HAS PSA HIT" : ""),
                       10, 0, 512, 10, -50, 400);
   fr->SetStats(0); fr->GetXaxis()->SetTitle("time bucket"); fr->GetYaxis()->SetTitle("ADC");
   fr->Draw(); fGarb.push_back(fr);
   if (thePad) {
      auto g = new TGraph();
      const auto &adc = thePad->GetADC();
      double mx = 60;
      for (int i = 0; i < 512; ++i) { g->SetPoint(i, i, adc[i]); if (adc[i] > mx) mx = adc[i]; }
      fr->GetYaxis()->SetRangeUser(-50, mx * 1.1);
      g->SetLineColor(kBlack); g->SetLineWidth(1); g->Draw("L same"); fGarb.push_back(g);
   }
   auto thr = new TLine(0, fThr, 512, fThr); thr->SetLineColor(kRed); thr->SetLineStyle(2); thr->Draw(); fGarb.push_back(thr);
   if (hasHit)
      for (double tb : fHitTB[pad]) {
         auto l = new TLine(tb, -50, tb, 400); l->SetLineColor(kGreen + 2); l->SetLineStyle(3); l->Draw(); fGarb.push_back(l);
      }
   c->Update();
   fInfo->SetText(Form(" event %lld/%lld   pad %d   %s ", fEvIdx, fNEv - 1, pad, hasHit ? "(green = PSA hit)" : "(no hit)"));
}

TraceViewer::TraceViewer(const char *file, double thr)
   : TGMainFrame(gClient->GetRoot(), 1200, 640), fEvIdx(0), fPadIdx(0), fThr(thr)
{
   fFile = TFile::Open(file);
   fTree = (TTree *)fFile->Get("cbmsim");
   fRaw = nullptr; fEv = nullptr;
   fTree->SetBranchAddress("AtRawEvent", &fRaw);
   fTree->SetBranchAddress("AtEventH", &fEv);
   fNEv = fTree->GetEntries();
   printf("opened %s : %lld events\n", file, fNEv);

   TGHorizontalFrame *bar = new TGHorizontalFrame(this, 1200, 40);
   auto add = [&](const char *lbl, const char *slot) {
      auto b = new TGTextButton(bar, lbl);
      b->Connect("Clicked()", "TraceViewer", this, slot);
      bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
   };
   add(" < Ev ", "PrevEv()"); add(" Ev > ", "NextEv()");
   bar->AddFrame(new TGLabel(bar, "ev:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 2, 4, 4));
   fEvNum = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger);
   fEvNum->Connect("ValueSet(Long_t)", "TraceViewer", this, "GotoEv()");
   fEvNum->GetNumberEntry()->Connect("ReturnPressed()", "TraceViewer", this, "GotoEv()");
   bar->AddFrame(fEvNum, new TGLayoutHints(kLHintsLeft, 2, 12, 4, 4));
   add(" < Pad ", "PrevPad()"); add(" Pad > ", "NextPad()");
   bar->AddFrame(new TGLabel(bar, "pad:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 4, 2, 4, 4));
   fPadNum = new TGNumberEntry(bar, 0, 6, -1, TGNumberFormat::kNESInteger);
   fPadNum->Connect("ValueSet(Long_t)", "TraceViewer", this, "GotoPad()");
   fPadNum->GetNumberEntry()->Connect("ReturnPressed()", "TraceViewer", this, "GotoPad()");
   bar->AddFrame(fPadNum, new TGLayoutHints(kLHintsLeft, 2, 12, 4, 4));
   fInfo = new TGLabel(bar, " step events & pads; red = threshold, green = PSA hit ");
   bar->AddFrame(fInfo, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 8, 8, 4, 4));
   AddFrame(bar, new TGLayoutHints(kLHintsExpandX));

   fEC = new TRootEmbeddedCanvas("ec", this, 1200, 590);
   AddFrame(fEC, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   SetWindowName("AT-TPC pad trace viewer");
   MapSubwindows(); Resize(GetDefaultSize()); MapWindow();
   LoadEvent();
}

void view_traces(const char *file = "/mnt/f/a1975/reco_d2/run_0300_psa_max.root", double thr = 20)
{
   new TraceViewer(file, thr);
}
