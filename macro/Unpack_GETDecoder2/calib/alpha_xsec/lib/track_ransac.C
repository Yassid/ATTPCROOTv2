// Tracking with ATTPCROOT's own pattern recognition.
//
// Runs SampleConsensus::AtSampleConsensus -- the method of Ayyad et al., NIM A 880 (2018)
// 166 -- over the unpacked events and writes one CSV of tracked hits per run.  This is the
// framework's own track finder, as opposed to the standalone Python clustering in
// lib/clustering.py; running both over the same runs is the cross-check that the tracking
// step otherwise lacks.
//
//   source <ATTPCROOT>/setup_fr19port.sh      # puts the headers on ROOT_INCLUDE_PATH
//   root -l -b -q 'track_ransac.C("in.root","out.csv")'
//   root -l -b -q 'track_ransac.C("in.root","out.csv",500)'   # first 500 events
//
// Output: event,track,pad,tb,x,y,q
//   one row per hit, with `track` saying which track the finder assigned it to.
//   Downstream code groups by (event, track).
#include "AtEvent.h"
#include "AtHit.h"
#include "AtSampleConsensus.h"
#include "AtTrack.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <vector>

void track_ransac(const char *in, const char *out, Long64_t nEvents = -1,
                  int iterations = 500, int minHits = 20, double distThresh = 10.0)
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("ERROR no cbmsim tree\n"); return; }

   TClonesArray *evArray = nullptr;
   t->SetBranchAddress("AtEventH", &evArray);
   Long64_t nAll = t->GetEntries();
   Long64_t nRun = (nEvents < 0 || nEvents > nAll) ? nAll : nEvents;

   SampleConsensus::AtSampleConsensus sc;
   sc.SetNumIterations(iterations);
   sc.SetMinHitsPattern(minHits);
   sc.SetDistanceThreshold(distThresh);

   std::ofstream o(out);
   o << "event,track,pad,tb,x,y,q\n";
   Long64_t nTracked = 0, nTracks = 0;
   for (Long64_t i = 0; i < nRun; ++i) {
      t->GetEntry(i);
      if (!evArray || evArray->GetEntriesFast() == 0) continue;
      auto *event = (AtEvent *)evArray->At(0);
      if (!event) continue;

      Int_t nh = event->GetNumHits();
      if (nh < minHits) continue;

      // AtSampleConsensus takes observer pointers; the hits stay owned by the event.
      std::vector<const AtHit *> hits;
      hits.reserve(nh);
      for (Int_t j = 0; j < nh; ++j)
         if (auto *h = event->GetHit(j)) hits.push_back(h);
      if (hits.size() < (size_t)minHits) continue;

      auto tracks = sc.SolveTracks(hits);
      if (tracks.empty()) continue;
      ++nTracked;
      for (size_t k = 0; k < tracks.size(); ++k) {
         auto *ha = tracks[k].GetHitArray();
         if (!ha) continue;
         ++nTracks;
         for (auto &h : *ha) {
            TVector3 p = h.GetPosition();
            o << i << "," << k << "," << h.GetHitPadNum() << "," << h.GetTimeStamp() << ","
              << Form("%.3f,%.3f,%.4g", p.X(), p.Y(), h.GetCharge()) << "\n";
         }
      }
      if (i % 500 == 0) { printf("\r  %lld/%lld  tracked %lld   ", i, nRun, nTracked); fflush(stdout); }
   }
   o.close();
   printf("\r  %s : %lld events, %lld with tracks, %lld tracks -> %s\n",
          gSystem->BaseName(in), nRun, nTracked, nTracks, gSystem->BaseName(out));
}

void track_ransac() { printf("usage: root -l -b -q 'track_ransac.C(\"in.root\",\"out.csv\")'\n"); }
