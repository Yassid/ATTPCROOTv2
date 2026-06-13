/// @file perf_distcut.C
/// @brief Confirm AtGenfit's vertex-to-Z distance cut (distThres=150mm) is what
/// drops gated deuterons. Replicates AtGenfit's dist = |xy of the vertex-end
/// cluster| and correlates dist<150 vs >150 with whether genfit fitted the track.
///   root -b -q 'pd/perf/perf_distcut.C("run_0106")'
#include <set>
#include <map>
void perf_distcut(TString run = "run_0106")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   auto pid = AtTools::AtParticleID::LoadJSON("pid/deuteron_band.json");

   std::set<std::pair<int, int>> gffit;
   {
      TFile *f = TFile::Open("/mnt/f/a1975/reco_gf/" + run + "_genfit.root");
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *te = nullptr; t->SetBranchAddress("AtTrackingEvent", &te);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i); if (te->GetEntries() == 0) continue; auto *ev = (AtTrackingEvent *)te->At(0); if (!ev) continue;
         for (auto &ft : ev->GetFittedTracks()) if (ft) { double ke = ft->GetKinematics().kineticEnergy;
            if (ke > 0 && ke < 1000) gffit.insert({(int)i, ft->GetTrackID()}); }
      }
      f->Close();
   }
   TFile *fr = TFile::Open("/mnt/f/a1975/reco/" + run + "_reco.root");
   TTree *tr = (TTree *)fr->Get("cbmsim"); TClonesArray *pe = nullptr; tr->SetBranchAddress("AtPatternEvent", &pe);
   TFile *fp = TFile::Open("/mnt/f/a1975/reco_gf/" + run + "_pid.root");
   TTree *tp = (TTree *)fp->Get("cbmsim"); TClonesArray *pide = nullptr; tp->SetBranchAddress("AtPIDEvent", &pide);

   long gated = 0, dlt = 0, dgt = 0, dlt_fit = 0, dgt_fit = 0;
   Long64_t N = std::min(tr->GetEntries(), tp->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      tr->GetEntry(i); tp->GetEntry(i);
      if (pe->GetEntries() == 0 || pide->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)pe->At(0); auto *pv = (AtPIDEvent *)pide->At(0); if (!pat || !pv) continue;
      std::map<int, const AtTools::AtSpyralResult *> pm; for (auto &sr : pv->GetSpyral()) pm[sr.trackID] = &sr;
      for (auto &tk : pat->GetTrackCand()) {
         int tid = tk.GetTrackID(); auto ip = pm.find(tid); if (ip == pm.end()) continue;
         auto *r = ip->second; if (!r->valid || !pid.IsInside(r->sqrtdEdx, r->brho)) continue;
         auto *hc = tk.GetHitClusterArray(); if (hc->size() < 3) continue;
         double thc = tk.GetGeoTheta() * TMath::RadToDeg();
         ROOT::Math::XYZPoint ini = (thc < 90) ? hc->back().GetPosition() : hc->front().GetPosition();
         double dist = std::sqrt(ini.X() * ini.X() + ini.Y() * ini.Y());
         ++gated; bool fit = gffit.count({(int)i, tid});
         if (dist < 150) { ++dlt; if (fit) ++dlt_fit; }
         else { ++dgt; if (fit) ++dgt_fit; }
      }
   }
   printf("RESULT gated=%ld | dist<150mm: %ld (genfit-fit %ld=%.0f%%) | dist>150mm: %ld (genfit-fit %ld=%.0f%%)\n",
          gated, dlt, dlt_fit, 100.0 * dlt_fit / std::max(1L, dlt), dgt, dgt_fit, 100.0 * dgt_fit / std::max(1L, dgt));
}
