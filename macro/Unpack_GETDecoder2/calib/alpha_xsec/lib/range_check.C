// Is the range TABLE right, or is the gas density wrong?
//
// In the data these two are degenerate: both scale every range by one factor, so no
// measurement there can separate them.  The simulation can, because it runs at a KNOWN
// density (299.5 torr, set by the geometry file) with Geant4's own stopping power, and it
// records both the initial energy of each track and the path length it actually travelled
// (AtTpcPoint::fLength).
//
// So: take each simulated product, read its true initial energy and its true range, and
// compare with what the analysis range table predicts at 299.5 torr.
//   * agree  -> the table is right, and the data's longer arms mean a genuinely thinner gas
//   * differ -> the table is wrong, and the density is not the problem
//
//   root -l -b -q 'range_check.C("sim.root","range_check.csv")'
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <fstream>
#include <map>
#include <vector>

void range_check(const char *in, const char *out)
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("ERROR no cbmsim\n"); return; }
   t->SetEstimate(-1);
   auto col = [&](const char *e) {
      Long64_t n = t->Draw(e, "", "goff");
      if (n > t->GetEstimate()) { t->SetEstimate(n + 1); n = t->Draw(e, "", "goff"); }
      return std::vector<double>(t->GetV1(), t->GetV1() + n);
   };
   auto pev = col("Entry$+0*AtTpcPoint.fLength");
   auto pid = col("AtTpcPoint.fTrackID");
   auto plen = col("AtTpcPoint.fLength");
   auto mev = col("Entry$+0*MCTrack.fPx");
   auto px = col("MCTrack.fPx"), py = col("MCTrack.fPy"), pz = col("MCTrack.fPz");

   // longest path reached by each (entry, trackID) = that track's RANGE
   std::map<std::pair<long, long>, double> rng;
   for (size_t i = 0; i < pev.size(); ++i) {
      auto k = std::make_pair((long)pev[i], (long)pid[i]);
      double L = plen[i];
      auto it = rng.find(k);
      if (it == rng.end() || L > it->second) rng[k] = L;
   }
   // initial energy of each MCTrack, indexed in order within its entry
   const double M = 3727.379;  // MeV
   std::map<std::pair<long, long>, double> ene;
   std::map<long, long> idx;
   for (size_t i = 0; i < mev.size(); ++i) {
      long e = (long)mev[i];
      long j = idx[e]++;
      double p = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]) * 1000.0;  // GeV -> MeV
      ene[std::make_pair(e, j)] = sqrt(p*p + M*M) - M;
   }
   std::ofstream o(out);
   o << "entry,trackid,E_MeV,range_raw\n";
   long n = 0;
   for (auto &kv : rng) {
      auto it = ene.find(kv.first);
      if (it == ene.end()) continue;
      if (it->second <= 0.05) continue;
      o << kv.first.first << "," << kv.first.second << ","
        << Form("%.5f,%.5f", it->second, kv.second) << "\n";
      ++n;
   }
   o.close();
   printf("  %s : %ld tracks -> %s\n", gSystem->BaseName(in), n, gSystem->BaseName(out));
}
