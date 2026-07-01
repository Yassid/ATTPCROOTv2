// Find events where AtPSAMultiFit recovered a 2nd hit on shared pads (overlapping tracks).
// Dumps: best event cloud (with primary/secondary flag) + one overlapping pad's raw trace+hits.
// Run: root -l 'hunt_overlaps.C'
#include <map>
#include <vector>
#include <fstream>
#include <algorithm>

void hunt_overlaps()
{
   gSystem->Load("libAtReconstruction.so");
   auto f = TFile::Open("/mnt/f/a1975/reco_d2/run_0016_psa_multifit.root");
   auto t = (TTree *)f->Get("cbmsim");
   TClonesArray *evArr = nullptr, *rawArr = nullptr;
   t->SetBranchAddress("AtEventH", &evArr);
   t->SetBranchAddress("AtRawEvent", &rawArr);
   long n = t->GetEntries();

   // find event with most multi-hit pads
   int bestEv = -1, bestMulti = -1;
   for (long i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto ev = (AtEvent *)evArr->At(0);
      if (!ev) continue;
      std::map<int, int> pm;
      for (const auto &h : ev->GetHits()) pm[h->GetPadNum()]++;
      int m = 0;
      for (auto &kv : pm) if (kv.second > 1) m++;
      if (m > bestMulti) { bestMulti = m; bestEv = i; }
   }
   printf("\nBest event = %d  (%d multi-hit pads)\n", bestEv, bestMulti);

   t->GetEntry(bestEv);
   auto ev = (AtEvent *)evArr->At(0);
   auto raw = (AtRawEvent *)rawArr->At(0);

   // group hits by pad
   std::map<int, std::vector<const AtHit *>> byPad;
   for (const auto &h : ev->GetHits()) byPad[h->GetPadNum()].push_back(h.get());

   // dump cloud with rank (0 = primary/highest amp on the pad, >=1 = recovered extra)
   std::ofstream cl("overlap_cloud.csv");
   cl << "x,y,z,pad,mult,rank\n";
   for (auto &kv : byPad) {
      auto hits = kv.second;
      std::sort(hits.begin(), hits.end(), [](auto a, auto b) { return a->GetCharge() > b->GetCharge(); });
      for (size_t r = 0; r < hits.size(); ++r) {
         auto h = hits[r];
         cl << h->GetPosition().X() << "," << h->GetPosition().Y() << "," << h->GetPosition().Z() << ","
            << h->GetPadNum() << "," << hits.size() << "," << r << "\n";
      }
   }
   cl.close();

   // pick the 2-hit pad with the LARGEST peak separation (clearest overlap) for a trace plot
   int exPad = -1;
   double bestSep = -1;
   for (auto &kv : byPad)
      if (kv.second.size() == 2) {
         double sep = std::abs(kv.second[0]->GetTimeStamp() - kv.second[1]->GetTimeStamp());
         if (sep > bestSep) { bestSep = sep; exPad = kv.first; }
      }
   printf("Example overlapping pad = %d  (2 hits, peak sep %.0f TB)\n", exPad, bestSep);

   if (exPad >= 0) {
      std::ofstream tr("overlap_trace.csv");
      tr << "tb,adc\n";
      for (const auto &p : raw->GetPads())
         if (p->GetPadNum() == exPad) {
            const auto &a = p->GetADC();
            for (int k = 0; k < 512; ++k) tr << k << "," << a[k] << "\n";
         }
      tr.close();
      std::ofstream hi("overlap_hits.csv");
      hi << "peakTB,amp,z\n";
      for (auto h : byPad[exPad]) hi << h->GetTimeStamp() << "," << h->GetCharge() << "," << h->GetPosition().Z() << "\n";
      hi.close();
      printf("dumped overlap_cloud.csv, overlap_trace.csv, overlap_hits.csv\n");
   }
}
