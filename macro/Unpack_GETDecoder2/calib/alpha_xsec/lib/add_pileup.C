// Overlay an UNREACTED beam track on the merged events, to emulate pile-up.
//
// Why this exists.  A merged simulated event contains exactly one beam particle: the one
// that reacted, truncated at its own vertex.  A real recorded event does not.  The trigger
// opens a window on a continuous beam, so other beam particles cross the chamber in the
// same frame and run the FULL length, straight through the region where the reacting
// beam's short stub sits.
//
// That difference is not cosmetic -- it is the one failure mode the efficiency measurement
// was missing.  In the data, 5.6% of selected events come out ABOVE the kinematic ceiling
// E_beam/2, reaching 9.86 MeV; the simulation puts 0.0% there.  Those events have a median
// beam path of 74.9 mm against 349 mm and only 38% carry a beam track against 90%: near the
// window the beam stub is too short to be found, the clustering absorbs it into one of the
// arms, and that arm comes out too long.  Extra beam hits through that region are exactly
// what would cause it, and the simulation has none.
//
// So eps_reco measured without this is an UPPER BOUND.  Running the chain with and without
// it brackets the systematic.
//
//   root -l -b -q 'add_pileup.C("merged.root","digi.root","pileup.root",1.0)'
//
//   frac = 0.0  reproduces the input;  1.0 puts one extra beam in every event.
#include "AtEvent.h"
#include "AtHit.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

void add_pileup(const char *merged, const char *digi, const char *out,
                double frac = 1.0, int seed = 12345)
{
   TFile *fm = TFile::Open(merged);
   if (!fm || fm->IsZombie()) { printf("ERROR cannot open %s\n", merged); return; }
   TTree *tm = (TTree *)fm->Get("cbmsim");
   TFile *fd = TFile::Open(digi);
   if (!fd || fd->IsZombie()) { printf("ERROR cannot open %s\n", digi); return; }
   TTree *td = (TTree *)fd->Get("cbmsim");
   if (!tm || !td) { printf("ERROR missing cbmsim tree\n"); return; }

   TClonesArray *mArr = nullptr, *dArr = nullptr;
   tm->SetBranchAddress("AtEventH", &mArr);
   td->SetBranchAddress("AtEventH", &dArr);

   TFile fo(out, "RECREATE");
   TTree *to = new TTree("cbmsim", "merged events with beam pile-up");
   auto *oArr = new TClonesArray("AtEvent");
   to->Branch("AtEventH", &oArr);

   TRandom3 rnd(seed);
   Long64_t nm = tm->GetEntries(), nd = td->GetEntries();
   Long64_t nBeamEv = nd / 2;                 // the EVEN entries are the beam-only halves
   long added = 0, extraHits = 0, baseHits = 0;

   for (Long64_t i = 0; i < nm; ++i) {
      tm->GetEntry(i);
      if (!mArr || mArr->GetEntriesFast() == 0) continue;
      auto *mev = (AtEvent *)mArr->At(0);
      if (!mev) continue;

      oArr->Clear("C");
      auto *ev = new ((*oArr)[0]) AtEvent();
      ev->SetEventID(i);
      for (Int_t j = 0; j < mev->GetNumHits(); ++j)
         if (auto *h = mev->GetHit(j)) { ev->AddHit(h); ++baseHits; }

      if (frac > 0 && rnd.Uniform() < frac && nBeamEv > 1) {
         // a DIFFERENT beam traverse, kept at full length -- it did not react
         Long64_t k = (i + 1 + (Long64_t)(rnd.Uniform() * (nBeamEv - 1))) % nBeamEv;
         td->GetEntry(2 * k);
         if (dArr && dArr->GetEntriesFast() > 0) {
            if (auto *bev = (AtEvent *)dArr->At(0)) {
               for (Int_t j = 0; j < bev->GetNumHits(); ++j)
                  if (auto *h = bev->GetHit(j)) { ev->AddHit(h); ++extraHits; }
               ++added;
            }
         }
      }
      to->Fill();
   }
   to->Write();
   fo.Close();
   printf("  %lld events, pile-up beam added to %ld (%.0f%%)\n", nm, added,
          nm ? 100.0 * added / nm : 0.0);
   printf("  hits/event: %.0f base + %.0f pile-up = %.0f   (real data: 140)\n",
          nm ? double(baseHits) / nm : 0.0, nm ? double(extraHits) / nm : 0.0,
          nm ? double(baseHits + extraHits) / nm : 0.0);
}

void add_pileup() { printf("usage: add_pileup.C(\"merged.root\",\"digi.root\",\"out.root\",frac)\n"); }
