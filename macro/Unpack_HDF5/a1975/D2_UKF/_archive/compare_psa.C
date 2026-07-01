// Compare AtPSAMultiFit vs AtPSAMax point clouds on real run_0016.
// Tallies hits + per-pad multiplicity, dumps a few events' clouds to CSV. Run: root -l 'compare_psa.C'
#include <map>
#include <fstream>

static void process(TString tag, TString file, std::ofstream &csv)
{
   auto f = TFile::Open(file);
   auto tree = (TTree *)f->Get("cbmsim");
   TClonesArray *evArr = nullptr;
   tree->SetBranchAddress("AtEventH", &evArr);

   long nEv = tree->GetEntries();
   long totHits = 0, totMultiPads = 0, totPads = 0, maxMult = 0;
   for (long i = 0; i < nEv; ++i) {
      tree->GetEntry(i);
      auto ev = (AtEvent *)evArr->At(0);
      if (!ev) continue;
      std::map<int, int> padMult;
      const auto &hits = ev->GetHits();
      for (const auto &h : hits) {
         padMult[h->GetPadNum()]++;
         // dump first 6 events with >=50 hits for visual continuity check
         if (i < 6)
            csv << tag << "," << i << "," << h->GetPosition().X() << "," << h->GetPosition().Y() << ","
                << h->GetPosition().Z() << "," << h->GetPadNum() << "\n";
      }
      totHits += hits.size();
      totPads += padMult.size();
      for (auto &kv : padMult) {
         if (kv.second > 1) totMultiPads++;
         if (kv.second > maxMult) maxMult = kv.second;
      }
   }
   printf("%-10s | events %3ld | total hits %6ld | pads-with-hits %6ld | MULTI-hit pads %5ld (%.1f%%) | "
          "hits/pad %.3f | max mult %ld\n",
          tag.Data(), nEv, totHits, totPads, totMultiPads, 100.0 * totMultiPads / std::max(1L, totPads),
          (double)totHits / std::max(1L, totPads), maxMult);
   f->Close();
}

void compare_psa()
{
   gSystem->Load("libAtReconstruction.so");
   std::ofstream csv("psa_clouds.csv");
   csv << "psa,event,x,y,z,pad\n";
   printf("\n=== AtPSAMax vs AtPSAMultiFit on real a1975 run_0016 ===\n");
   process("max", "/mnt/f/a1975/reco_d2/run_0016_psa_max.root", csv);
   process("multifit", "/mnt/f/a1975/reco_d2/run_0016_psa_multifit.root", csv);
   csv.close();
   printf("\nwrote psa_clouds.csv (first 6 events, both PSAs) for visual comparison\n");
}
