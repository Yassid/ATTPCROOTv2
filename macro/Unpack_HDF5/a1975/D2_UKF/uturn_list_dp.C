/// @file uturn_list_dp.C
/// @brief Write a GUI list file of the WORST z U-TURNS, worst first, so the failure can be paged
///        through in gui_fit_dp.C instead of argued about one event at a time.
///
/// Same metric as zuturn_dp.C: total variation of z along the STORED cluster array over the z
/// span.  1.0 is a monotonic traversal (what the simulation gives, median 1.00), 2.0 is one full
/// U-turn.  The fitted chi2/ndf is carried alongside so it is visible whether the ordering damage
/// actually reached the fit -- a track can U-turn and still fit, and those are as interesting as
/// the ones that fail.
///
/// SORTED BY tv ASCENDING WITHIN A BAND, not worst-first.  Ranking by tv put tv = 36, 10 and 9 at
/// the top -- tracks that are pathological rather than typical, and looking at those tells you
/// about the tail, not about the failure.  The median of the U-turning population is ~1.9, so the
/// band defaults to 1.7-2.5 and the list opens on a REPRESENTATIVE one.  `c2Min` additionally
/// demands that a fit exists and is bad, since a track with no fit object shows nothing in the
/// GUI's fitted-trajectory panels.
///
///   root -l -b -q 'uturn_list_dp.C()'                                  typical failures
///   root -l -b -q 'uturn_list_dp.C(...,800,0,1e9,-1)'                  everything, in tv order
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
void uturn_list_dp(TString listFile = "/mnt/f/a1975/caches/.explorer/longspiral_list.txt",
                   TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                   TString outFile = "/mnt/f/a1975/caches/.explorer/uturn_list.txt",
                   int maxTracks = 800, double tvLo = 1.7, double tvHi = 2.5,
                   double c2Min = 5.0, double vtxMin = -1)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while (in >> r >> e >> ti) jobs.push_back({TString(r.c_str()), e, ti});

   struct Row { TString run; Long64_t ev; int tid; double tv; int ncl; double c2; double vtx; int npt; };
   std::vector<Row> v;
   TString open; TFile *f = nullptr; TTree *t = nullptr; TClonesArray *a = nullptr;
   int done = 0;
   for (auto &j : jobs) {
      if (done >= maxTracks) break;
      TString rn = std::get<0>(j); Long64_t en = std::get<1>(j); int id = std::get<2>(j);
      TString rt = "run_" + rn + "_multifit";
      if (rt != open) {
         if (f) f->Close();
         f = TFile::Open(gfDir + rt + "_genfitter_p.root");
         if (!f || f->IsZombie()) { f = nullptr; continue; }
         t = (TTree *)f->Get("cbmsim"); a = nullptr;
         t->SetBranchAddress("AtTrackingEvent", &a); open = rt;
      }
      t->GetEntry(en);
      auto *ev = (AtTrackingEvent *)a->At(0); if (!ev) continue;
      // the fitted chi2/ndf for this track, if the fit exists at all
      double c2 = -1, vtx = -1; int npt = -1;
      for (auto &ft : ev->GetFittedTracks())   // unique_ptr, so dereference it
         if (ft && ft->GetTrackID() == id) {
            // chi2/ndf lives in the metadata, not on AtFittedTrack itself -- same access
            // ex_dp_dv1104.C uses, so this column means the same thing as the cache's
            double ndf = ft->GetTrackMetadata()->GetNdf();
            if (ndf > 0) c2 = ft->GetTrackMetadata()->GetChi2() / ndf;
            // THE CRITERION CHI2 CANNOT FAKE.  chi2/ndf falls mechanically as clusters are
            // removed, so after truncation it stops being a failure detector; the back-extrapolated
            // vertex has to land on the BEAM AXIS whatever the fit was done on, and a truncated fit
            // that is merely self-consistent misses it just as badly as the full-length one did.
            auto v = ft->GetVertex();
            vtx = std::hypot(v.X(), v.Y());
            npt = (int)ft->GetSmoothedPositions().size();
         }
      for (auto &tr : ev->GetTrackArray()) {
         if (tr.GetTrackID() != id) continue;
         auto *hc = tr.GetHitClusterArray(); const int n = hc->size();
         if (n < 20) continue;
         double s = 0, lo = 1e9, hi = -1e9;
         for (int i = 0; i < n; ++i) { double z = (*hc)[i].GetPosition().Z();
            lo = std::min(lo, z); hi = std::max(hi, z);
            if (i) s += std::fabs(z - (*hc)[i-1].GetPosition().Z()); }
         if (hi - lo < 20) continue;
         v.push_back({rn, en, id, s/(hi-lo), n, c2, vtx, npt});
         ++done;
      }
   }
   if (f) f->Close();
   std::vector<Row> sel;
   for (auto &x : v)
      if (x.tv >= tvLo && x.tv <= tvHi &&
          (x.c2 >= c2Min || (vtxMin > 0 && x.vtx >= vtxMin))) sel.push_back(x);
   std::sort(sel.begin(), sel.end(), [](const Row &A, const Row &B){ return A.tv < B.tv; });
   std::ofstream o(outFile.Data());
   for (auto &x : sel) o << x.run << " " << x.ev << " " << x.tid << "\n";
   o.close();
   printf("\n  %zu tracks measured; %zu in the band tv %.2f-%.2f with chi2/ndf >= %.1f -> %s\n",
          v.size(), sel.size(), tvLo, tvHi, c2Min, outFile.Data());
   std::vector<Row> &w = sel.empty() ? v : sel;
   printf("\n  the first 15 (this is the order the GUI will page through):\n");
   printf("    %-6s %8s %4s %8s %7s %10s %8s %8s\n", "run", "entry", "tid", "z tv", "ncl",
          "chi2/ndf", "fitpts", "vtx|xy|");
   for (size_t i = 0; i < w.size() && i < 15; ++i)
      printf("    %-6s %8lld %4d %8.2f %7d %10.2f %8d %8.1f\n", w[i].run.Data(), w[i].ev,
             w[i].tid, w[i].tv, w[i].ncl, w[i].c2, w[i].npt, w[i].vtx);
   // where the agreed test case sits in this ranking
   for (size_t i = 0; i < v.size(); ++i)
      if (v[i].run == "0016" && v[i].ev == 31316)
         printf("\n  run_0016 e31316 (the agreed case) is rank %zu of %zu, tv %.2f, chi2/ndf %.2f\n",
                i+1, v.size(), v[i].tv, v[i].c2);
}
