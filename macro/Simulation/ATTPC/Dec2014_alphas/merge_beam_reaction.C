// Merge the alternating beam / reaction entries into events that look like recorded ones.
//
// AtTPC2Body via AtVertexPropagator alternates: EVEN entries hold the pre-generated beam
// traverse used to choose the vertex, ODD entries hold the reaction products.  So a
// reaction event carries the two outgoing alphas and no beam segment at all, while a real
// recorded event has the beam running from the window down to the vertex and then the two
// arms.  Any analysis that wants beam + two arms -- as the alpha cross-section chain does --
// must put the halves back together.
//
// For each pair (2i, 2i+1) this writes one event holding
//    * the beam hits from entry 2i, TRUNCATED at the vertex, since in reality the beam
//      stops there, and
//    * all the product hits from entry 2i+1.
// The vertex is taken as the highest-z product hit: the products start at the vertex and
// travel towards the pad plane, so that is where the beam ended.
//
//   root -l -b -q 'merge_beam_reaction.C("output_digi.root","merged.root")'
#include "AtEvent.h"
#include "AtHit.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>

void merge_beam_reaction(const char *in, const char *out, double zPad = 5.0)
{
   TFile *fi = TFile::Open(in);
   if (!fi || fi->IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *ti = (TTree *)fi->Get("cbmsim");
   if (!ti) { printf("ERROR no cbmsim tree\n"); return; }
   TClonesArray *inArr = nullptr;
   ti->SetBranchAddress("AtEventH", &inArr);

   TFile fo(out, "RECREATE");
   TTree *to = new TTree("cbmsim", "beam and reaction merged");
   auto *outArr = new TClonesArray("AtEvent");
   to->Branch("AtEventH", &outArr);

   Long64_t n = ti->GetEntries(), made = 0;
   long beamKept = 0, beamCut = 0, prod = 0;
   for (Long64_t i = 0; i + 1 < n; i += 2) {
      // --- the reaction half first: it defines the vertex ---
      ti->GetEntry(i + 1);
      if (!inArr || inArr->GetEntriesFast() == 0) continue;
      auto *rEv = (AtEvent *)inArr->At(0);
      if (!rEv || rEv->GetNumHits() == 0) continue;

      std::vector<AtHit> products;
      double zVertex = -1e9;
      for (Int_t j = 0; j < rEv->GetNumHits(); ++j) {
         auto *h = rEv->GetHit(j);
         if (!h) continue;
         products.push_back(*h);
         zVertex = std::max(zVertex, h->GetPosition().Z());
      }
      if (products.empty()) continue;

      // --- the beam half, cut at the vertex ---
      std::vector<AtHit> beam;
      ti->GetEntry(i);
      if (inArr && inArr->GetEntriesFast() > 0) {
         if (auto *bEv = (AtEvent *)inArr->At(0)) {
            for (Int_t j = 0; j < bEv->GetNumHits(); ++j) {
               auto *h = bEv->GetHit(j);
               if (!h) continue;
               // keep only what lies upstream of the vertex: the beam stopped there
               if (h->GetPosition().Z() > zVertex - zPad) { beam.push_back(*h); ++beamKept; }
               else ++beamCut;
            }
         }
      }

      outArr->Clear("C");
      auto *ev = new ((*outArr)[0]) AtEvent();
      ev->SetEventID(made);
      for (auto &h : beam) ev->AddHit(&h);   // fr19port AtEvent::AddHit takes a pointer
      for (auto &h : products) { ev->AddHit(&h); ++prod; }
      to->Fill();
      ++made;
   }
   to->Write();
   fo.Close();
   printf("  %lld pairs -> %lld merged events\n", n / 2, made);
   printf("  beam hits kept %ld, dropped past the vertex %ld; product hits %ld\n",
          beamKept, beamCut, prod);
   printf("  mean %.0f hits/event  (real data: 154)\n", made ? double(beamKept + prod) / made : 0);
}

void merge_beam_reaction() { printf("usage: merge_beam_reaction.C(\"in.root\",\"out.root\")\n"); }

/* RESULT, 2026-10-23, 200 simulated events -> 100 merged:
 *
 *     mean hits/event        131   (real data 154)
 *     beam hits kept         7.4 per event, none dropped past the vertex
 *     events with >=3 tracks   9/100 at minHits=10, against 32% in real data
 *
 * The merge works, but THE SIMULATED BEAM IS TOO THIN: it leaves ~7 hits where the real
 * beam leaves ~21-27, so it falls below the track finder's minimum and is not found as a
 * segment at all.  Lowering minHits confirms it is nothing else:
 *
 *     minHits 10 ->  9/100     minHits 6 -> 32/100     minHits 4 -> 61/100
 *
 * and 32/100 at minHits=6 is exactly what the real data gives at minHits=10.
 *
 * The likely cause is the beam being generated as a pencil: SetSpotRadius(0, -100, 0) in
 * He4He4_sim_el_11p4.C gives it no spot size and no divergence, so it crosses far fewer
 * pads than a real beam with finite emittance.  Give it the measured spot and divergence
 * before taking any efficiency from this -- a beam that the tracker cannot see makes the
 * three-segment efficiency look far worse than it is.
 */
