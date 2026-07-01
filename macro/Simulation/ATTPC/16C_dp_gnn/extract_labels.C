// Build the labeled point cloud for GNN training. For each PSA hit, resolve its
// truth track ID via:  hit.padNum -> AtRawEvent SimMCPointMap -> mcPointID ->
// AtTpcPoint(AtMCPoint).GetTrackID().  Also records the particle (A,Z) of that track.
// Output CSV: event,x,y,z,charge,padNum,trackID,A,Z,shared
//   shared=1 if the pad had hits from >1 distinct track (ambiguous label).
//
//   root -l -b -q 'extract_labels.C("data/output_digi.root","data/labeled.csv")'
#include <set>
#include <map>
void extract_labels(const char *in = "data/output_digi.root", const char *out = "data/labeled.csv",
                    const char *mcBranch = "AtTpcPoint_1")
{
   TFile *f = TFile::Open(in);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *rawArr = nullptr, *evtArr = nullptr, *mcArr = nullptr;
   t->SetBranchAddress("AtRawEvent", &rawArr);
   t->SetBranchAddress("AtEventH", &evtArr);
   t->SetBranchAddress(mcBranch, &mcArr);

   FILE *fp = fopen(out, "w");
   // px,py,pz = MC truth momentum [GeV] of the contributing point (for lab-angle cuts)
   fprintf(fp, "event,x,y,z,charge,padNum,trackID,A,Z,shared,px,py,pz\n");
   Long64_t N = t->GetEntries();
   long nhits = 0, nlabeled = 0;
   std::set<int> trackset;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!rawArr || rawArr->GetEntriesFast() == 0 || !evtArr || evtArr->GetEntriesFast() == 0) continue;
      auto *raw = (AtRawEvent *)rawArr->At(0);
      auto *evt = (AtEvent *)evtArr->At(0);
      auto &mcmap = raw->GetSimMCPointMap();              // pad -> mcPointID (multimap)
      for (int h = 0; h < evt->GetNumHits(); ++h) {
         auto &hit = evt->GetHit(h);
         int pad = hit.GetPadNum();
         auto pos = hit.GetPosition();
         // gather distinct track IDs contributing to this pad
         std::map<int, int> trkCount; int A = -1, Z = -1, trackID = -1;
         double px = 0, py = 0, pz = 0;
         auto range = mcmap.equal_range(pad);
         for (auto it = range.first; it != range.second; ++it) {
            int mcid = it->second;
            if (mcid < 0 || mcid >= mcArr->GetEntriesFast()) continue;
            auto *mc = (AtMCPoint *)mcArr->At(mcid);
            int tid = mc->GetTrackID();
            trkCount[tid]++;
            if (trackID < 0) {
               trackID = tid; A = mc->GetMassNum(); Z = mc->GetAtomicNum();
               px = mc->GetPx(); py = mc->GetPy(); pz = mc->GetPz();
            }
         }
         int shared = trkCount.size() > 1 ? 1 : 0;
         ++nhits; if (trackID >= 0) { ++nlabeled; trackset.insert(trackID); }
         fprintf(fp, "%lld,%.3f,%.3f,%.3f,%.2f,%d,%d,%d,%d,%d,%.5f,%.5f,%.5f\n",
                 i, pos.X(), pos.Y(), pos.Z(), hit.GetCharge(), pad, trackID, A, Z, shared, px, py, pz);
      }
   }
   fclose(fp);
   printf("\033[1;32mDONE\033[0m %s : %lld events, %ld hits (%ld labeled = %.1f%%), %zu distinct trackIDs\n",
          out, N, nhits, nlabeled, 100.0 * nlabeled / std::max(1L, nhits), trackset.size());
}
