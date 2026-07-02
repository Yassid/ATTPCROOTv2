// display_compare.cpp  (standalone, links libAtData)
// Side-by-side event displays: EXP vs SIM point clouds (AtEventH).
// Rows = sample events; columns = [EXP x-y, EXP z-x, SIM x-y, SIM z-x].
// Build with build_width.sh-style flags. Run:
//   ./display_compare <exp.root> <sim.root> <nRows> <minHits> <out.png>
#include "AtEvent.h"
#include "AtHit.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLatex.h>

#include <vector>
#include <string>
#include <cstdio>

// collect up to nRows events (with >= minHits) from a file; returns graphs per event
struct Ev { std::vector<double> x, y, z; };

std::vector<Ev> grab(const char *fn, int nRows, int minHits)
{
   std::vector<Ev> out;
   TFile *f = TFile::Open(fn);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", fn); return out; }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *arr = nullptr;
   t->SetBranchAddress("AtEventH", &arr);
   long n = t->GetEntries();
   for (long i = 0; i < n && (int)out.size() < nRows; ++i) {
      t->GetEntry(i);
      if (!arr) continue;
      for (int k = 0; k < arr->GetEntriesFast() && (int)out.size() < nRows; ++k) {
         auto *ev = (AtEvent *)arr->At(k);
         if (!ev) continue;
         const auto &hits = ev->GetHits();
         if ((int)hits.size() < minHits) continue;
         Ev e;
         for (const auto &h : hits) {
            auto p = h->GetPosition();
            e.x.push_back(p.X()); e.y.push_back(p.Y()); e.z.push_back(p.Z());
         }
         out.push_back(std::move(e));
      }
   }
   // keep file open? we copied data out, safe to close
   f->Close();
   return out;
}

void drawScatter(const std::vector<double> &a, const std::vector<double> &b, const char *title,
                 double axmin, double axmax, double aymin, double aymax)
{
   TH2F *fr = new TH2F(Form("fr_%p", (void *)&a), title, 10, axmin, axmax, 10, aymin, aymax);
   fr->SetStats(0);
   fr->Draw();
   if (!a.empty()) {
      TGraph *g = new TGraph(a.size(), a.data(), b.data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.35);
      g->SetMarkerColor(kAzure + 2);
      g->Draw("P SAME");
   }
}

int main(int argc, char **argv)
{
   if (argc < 6) { printf("usage: %s exp.root sim.root nRows minHits out.png\n", argv[0]); return 1; }
   const char *expF = argv[1], *simF = argv[2], *outPng = argv[5];
   int nRows = atoi(argv[3]), minHits = atoi(argv[4]);

   gStyle->SetOptStat(0);
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

   auto expEv = grab(expF, nRows, minHits);
   auto simEv = grab(simF, nRows, minHits);
   int rows = std::min((int)expEv.size(), (int)simEv.size());
   if (rows == 0) { printf("no events with >=%d hits\n", minHits); return 1; }
   printf("drawing %d rows (exp %zu, sim %zu found)\n", rows, expEv.size(), simEv.size());

   TCanvas *c = new TCanvas("c", "exp vs sim", 1600, 300 * rows);
   c->Divide(4, rows, 0.001, 0.001);

   double R = 260, Zmax = 1000;
   for (int r = 0; r < rows; ++r) {
      c->cd(r * 4 + 1);
      drawScatter(expEv[r].x, expEv[r].y, Form("EXP #%d  x-y (%zu hits);x [mm];y [mm]", r, expEv[r].x.size()),
                  -R, R, -R, R);
      c->cd(r * 4 + 2);
      drawScatter(expEv[r].z, expEv[r].x, Form("EXP #%d  z-x;z [mm];x [mm]", r), 0, Zmax, -R, R);
      c->cd(r * 4 + 3);
      drawScatter(simEv[r].x, simEv[r].y, Form("SIM #%d  x-y (%zu hits);x [mm];y [mm]", r, simEv[r].x.size()),
                  -R, R, -R, R);
      c->cd(r * 4 + 4);
      drawScatter(simEv[r].z, simEv[r].x, Form("SIM #%d  z-x;z [mm];x [mm]", r), 0, Zmax, -R, R);
   }
   c->SaveAs(outPng);
   printf("wrote %s\n", outPng);
   return 0;
}
