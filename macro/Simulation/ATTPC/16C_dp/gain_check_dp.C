/// @file gain_check_dp.C
/// @brief Validate the SIMULATION GAIN for a1975 16C(d,p)17C against the data, by comparing what
/// a track is actually made of: hits per track and clusters per track.
///
/// WHY GAIN AND WHY THIS OBSERVABLE. Gain decides how many electrons reach a pad, hence whether a
/// pad crosses the PSA threshold, hence how many hits a track keeps -- which is exactly what an
/// ACCEPTANCE measures. Too low and short or low-dE/dx tracks are silently deleted and the
/// acceptance comes out too small with nothing in the output looking wrong. The (p,d) simulation
/// ran at 150000 and gave 117 hits/track against 46 in data, a factor 2.54 -- the long-standing
/// "sim and data diverge structurally" item turned out to be gain, not clustering.
///
/// THE (d,t) RESULT THIS STARTS FROM (run_reco_batch_dt.sh header, measured on gs_s3001 vs
/// run_0016, non-beam tracks):
///       gain      median hits/track   sim/data
///       10000            39             0.51     <- the production par's value
///       35000            75             0.99     <- ADOPTED for (d,t)
///      150000           118             1.55     <- the (p,d) par's value
/// PROTONS ARE NOT TRITONS -- lower Z, different dE/dx, different pad occupancy -- so 35000 is a
/// starting point for (d,p), not an inherited answer. That is what this macro is for.
///
/// SAME FOOTING, and it matters:
///   * NON-BEAM tracks only. The beam is a different object and dominates the hit count.
///   * tracks above maxHits (default 500) removed -- the (d,t) scan's convention, to drop beam
///     remnants and merged objects that no gain setting controls.
///   * MEDIANS, never means. A fit-failure or merge tail makes the mean say whatever it likes
///     (see feedback_use_medians).
///
///   root -l 'gain_check_dp.C("/mnt/f/a1975_C16_dp_sim/gs_s4001_reco.root")'

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <algorithm>
#include <vector>
#include <cstdio>

namespace {

/// median of a vector, by value
double med(std::vector<double> v)
{
   if (v.empty()) return -1;
   std::sort(v.begin(), v.end());
   const size_t n = v.size();
   return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
double quant(std::vector<double> v, double q)
{
   if (v.empty()) return -1;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
}

/// hits-per-track and clusters-per-track over the pattern tracks of one file
struct Prof {
   std::vector<double> hits, clus;
   Long64_t events = 0, tracks = 0, dropped = 0;
};

Prof profile(const char *fname, int maxHits, Long64_t maxEvents)
{
   Prof P;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) { printf("  ERROR: cannot open %s\n", fname); return P; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("  ERROR: no cbmsim in %s\n", fname); return P; }
   TClonesArray *pat = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pat);
   const Long64_t N = (maxEvents > 0) ? TMath::Min(maxEvents, t->GetEntries()) : t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!pat || pat->GetEntriesFast() == 0) continue;
      auto *pe = (AtPatternEvent *)pat->At(0);
      if (!pe) continue;
      ++P.events;
      auto &trks = pe->GetTrackCand();
      for (auto &tr : trks) {
         const int nh = (int)tr.GetHitArray().size();
         const int nc = (int)tr.GetHitClusterArray()->size();
         ++P.tracks;
         if (maxHits > 0 && nh > maxHits) { ++P.dropped; continue; }
         P.hits.push_back(nh);
         P.clus.push_back(nc);
      }
   }
   f->Close();
   return P;
}

} // namespace

void gain_check_dp(TString simFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_reco.root",
                   TString dataFile = "/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root",
                   Int_t maxHits = 500, Long64_t maxEventsData = 5000, Long64_t maxEventsSim = 0)
{
   printf("\n=== gain_check_dp : is the simulation gain right for PROTONS? ===\n");
   printf("  sim   %s\n  data  %s\n", simFile.Data(), dataFile.Data());
   printf("  non-beam pattern tracks, tracks above %d hits dropped, MEDIANS\n\n", maxHits);

   Prof S = profile(simFile, maxHits, maxEventsSim);
   Prof D = profile(dataFile, maxHits, maxEventsData);
   if (S.hits.empty() || D.hits.empty()) { printf("  nothing to compare\n"); return; }

   const double sh = med(S.hits), dh = med(D.hits);
   const double sc = med(S.clus), dc = med(D.clus);

   printf("                        events    tracks   dropped    median    p25    p75\n");
   printf("  hits/track   SIM     %7lld   %7lld   %7lld   %7.1f  %5.0f  %5.0f\n",
          S.events, S.tracks, S.dropped, sh, quant(S.hits, 0.25), quant(S.hits, 0.75));
   printf("  hits/track   DATA    %7lld   %7lld   %7lld   %7.1f  %5.0f  %5.0f\n",
          D.events, D.tracks, D.dropped, dh, quant(D.hits, 0.25), quant(D.hits, 0.75));
   printf("                                                 sim/data = %.2f\n\n", dh > 0 ? sh / dh : -1);
   printf("  clusters/trk SIM                               %7.1f  %5.0f  %5.0f\n",
          sc, quant(S.clus, 0.25), quant(S.clus, 0.75));
   printf("  clusters/trk DATA                              %7.1f  %5.0f  %5.0f\n",
          dc, quant(D.clus, 0.25), quant(D.clus, 0.75));
   printf("                                                 sim/data = %.2f\n\n", dc > 0 ? sc / dc : -1);

   const double r = dh > 0 ? sh / dh : -1;
   if (r < 0) printf("  VERDICT: could not compare.\n");
   else if (r > 0.8 && r < 1.25) printf("  VERDICT: hits/track within 25%% -- gain is in the right place.\n");
   else if (r < 0.8) printf("  VERDICT: SIM IS THIN (%.2f). Gain too LOW: pads are failing the PSA\n"
                            "           threshold, so short and low-dE/dx tracks are being deleted and\n"
                            "           any acceptance from this will be too small. RAISE the gain.\n", r);
   else printf("  VERDICT: SIM IS FAT (%.2f). Gain too HIGH -- the (p,d) failure mode exactly.\n"
               "           LOWER the gain.\n", r);
   printf("\n  CAVEAT: the data sample here is EVERY non-beam track in the run, not protons -- there\n"
          "  is no species selection without the PID gate. It is the same footing the (d,t) gain was\n"
          "  chosen on, and it is a comparison of chain output, not of a physics sample.\n\n");
}
