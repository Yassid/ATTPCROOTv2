/// @file perf_vertexend.C
/// @brief Diagnose the vertex-end mis-identification behind the distThres drops.
/// For each gated deuteron, compare the distance-to-beam-axis of the CHOSEN
/// vertex-end cluster (AtGenfit: back()/front() by GeoTheta) vs the track's CLOSEST
/// cluster to the axis. If the chosen end is far while a closer cluster exists, the
/// vertex (and seed) is mis-identified -> the real fix is the seed, not the cut.
///   root -b -q 'pd/perf/perf_vertexend.C("run_0106")'
#include <map>
void perf_vertexend(TString run = "run_0106")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   auto pid = AtTools::AtParticleID::LoadJSON("pid/deuteron_band.json");
   TFile *fr = TFile::Open("/mnt/f/a1975/reco/" + run + "_reco.root");
   TTree *tr = (TTree *)fr->Get("cbmsim"); TClonesArray *pe = nullptr; tr->SetBranchAddress("AtPatternEvent", &pe);
   TFile *fp = TFile::Open("/mnt/f/a1975/reco_gf/" + run + "_pid.root");
   TTree *tp = (TTree *)fp->Get("cbmsim"); TClonesArray *pide = nullptr; tp->SetBranchAddress("AtPIDEvent", &pide);

   long gated = 0, chosenFar = 0, chosenFar_closeExists = 0;
   double sChosen = 0, sMin = 0;
   Long64_t N = std::min(tr->GetEntries(), tp->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      tr->GetEntry(i); tp->GetEntry(i);
      if (pe->GetEntries() == 0 || pide->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)pe->At(0); auto *pv = (AtPIDEvent *)pide->At(0); if (!pat || !pv) continue;
      std::map<int, const AtTools::AtSpyralResult *> pm; for (auto &sr : pv->GetSpyral()) pm[sr.trackID] = &sr;
      for (auto &tk : pat->GetTrackCand()) {
         auto ip = pm.find(tk.GetTrackID()); if (ip == pm.end()) continue;
         auto *r = ip->second; if (!r->valid || !pid.IsInside(r->sqrtdEdx, r->brho)) continue;
         auto *hc = tk.GetHitClusterArray(); if (hc->size() < 3) continue;
         double thc = tk.GetGeoTheta() * TMath::RadToDeg();
         auto chosen = (thc < 90) ? hc->back().GetPosition() : hc->front().GetPosition();
         double dChosen = std::sqrt(chosen.X() * chosen.X() + chosen.Y() * chosen.Y());
         double dMin = 1e9; for (auto &cl : *hc) { auto p = cl.GetPosition(); double d = std::sqrt(p.X() * p.X() + p.Y() * p.Y()); if (d < dMin) dMin = d; }
         ++gated; sChosen += dChosen; sMin += dMin;
         if (dChosen > 150) { ++chosenFar; if (dMin < 150) ++chosenFar_closeExists; }
      }
   }
   printf("RESULT gated=%ld  <chosen-end dist>=%.0fmm  <min-cluster dist>=%.0fmm\n", gated, sChosen / gated, sMin / gated);
   printf("  chosen-end >150mm: %ld (%.0f%%)  -- of those, a cluster <150mm EXISTS: %ld (%.0f%%) => mis-identified vertex\n",
          chosenFar, 100.0 * chosenFar / gated, chosenFar_closeExists, 100.0 * chosenFar_closeExists / std::max(1L, chosenFar));
}
