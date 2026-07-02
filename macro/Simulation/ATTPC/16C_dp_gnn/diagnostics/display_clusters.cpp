// display_clusters.cpp  (standalone, links libAtData)
// Render HDBSCAN clusters (AtPatternEvent -> AtTrack) for EXP vs SIM, each cluster a
// distinct colour. Columns: [EXP x-y, EXP z-x, SIM x-y, SIM z-x], rows = events.
//   ./display_clusters <exp_hdb.root> <sim_hdb.root> <nRows> <minTracks> <out.png>
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtHit.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cstdio>

struct Trk { std::vector<double> x, y, z; };
struct Evt { std::vector<Trk> trks; };

std::vector<Evt> grab(const char *fn, int nRows, int minTracks)
{
   std::vector<Evt> out;
   TFile *f = TFile::Open(fn);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", fn); return out; }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *arr = nullptr;
   t->SetBranchAddress("AtPatternEvent", &arr);
   long n = t->GetEntries();
   for (long i = 0; i < n && (int)out.size() < nRows; ++i) {
      t->GetEntry(i);
      if (!arr || arr->GetEntriesFast() == 0) continue;
      auto *pe = (AtPatternEvent *)arr->At(0);
      if (!pe) continue;
      auto &tracks = pe->GetTrackCand();
      if ((int)tracks.size() < minTracks) continue;
      Evt e;
      for (auto &tr : tracks) {
         Trk k;
         auto &ha = tr.GetHitArray();
         for (auto &h : ha) {
            if (!h) continue;
            auto p = h->GetPosition();
            k.x.push_back(p.X()); k.y.push_back(p.Y()); k.z.push_back(p.Z());
         }
         if (!k.x.empty()) e.trks.push_back(std::move(k));
      }
      if (out.empty())
         printf("  [%s] evt %ld: %zu tracks, first-track hits=%zu\n", fn, i, tracks.size(),
                e.trks.empty() ? 0 : e.trks[0].x.size());
      out.push_back(std::move(e));
   }
   f->Close();
   return out;
}

int COL[10] = {kRed+1, kAzure+2, kGreen+2, kOrange+7, kMagenta+1, kCyan+2, kSpring+4, kViolet, kPink+7, kYellow+2};

void drawEvt(const Evt &e, const char *title, bool zx)
{
   double R = 260, Zmax = 1000;
   TH2F *fr = zx ? new TH2F(Form("f%p", (void *)&e), title, 10, 0, Zmax, 10, -R, R)
                 : new TH2F(Form("f%p", (void *)&e), title, 10, -R, R, 10, -R, R);
   fr->SetStats(0); fr->Draw();
   int ci = 0;
   for (auto &k : e.trks) {
      const std::vector<double> &a = zx ? k.z : k.x;
      const std::vector<double> &b = zx ? k.x : k.y;
      TGraph *g = new TGraph(a.size(), a.data(), b.data());
      g->SetMarkerStyle(20); g->SetMarkerSize(0.35); g->SetMarkerColor(COL[ci % 10]);
      g->Draw("P SAME");
      ci++;
   }
}

int main(int argc, char **argv)
{
   if (argc < 6) { printf("usage: %s exp_hdb.root sim_hdb.root nRows minTracks out.png\n", argv[0]); return 1; }
   int nRows = atoi(argv[3]), minTr = atoi(argv[4]);
   gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   auto E = grab(argv[1], nRows, minTr);
   auto S = grab(argv[2], nRows, minTr);
   int rows = std::min((int)E.size(), (int)S.size());
   if (rows == 0) { printf("no events with >=%d tracks\n", minTr); return 1; }
   printf("drawing %d rows (exp %zu, sim %zu)\n", rows, E.size(), S.size());
   TCanvas *c = new TCanvas("c", "clusters", 1600, 300 * rows);
   c->Divide(4, rows, 0.001, 0.001);
   for (int r = 0; r < rows; ++r) {
      c->cd(r * 4 + 1); drawEvt(E[r], Form("EXP #%d x-y (%zu clusters)", r, E[r].trks.size()), false);
      c->cd(r * 4 + 2); drawEvt(E[r], Form("EXP #%d z-x", r), true);
      c->cd(r * 4 + 3); drawEvt(S[r], Form("SIM #%d x-y (%zu clusters)", r, S[r].trks.size()), false);
      c->cd(r * 4 + 4); drawEvt(S[r], Form("SIM #%d z-x", r), true);
   }
   c->SaveAs(argv[5]);
   printf("wrote %s\n", argv[5]);
   return 0;
}
