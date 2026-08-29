/// @file dp_purity_C14.C
/// @brief What are the extra hits in the "fat" backward tracks?
///
///   root -b -q 'dp_purity_C14.C()'
///
/// Two backward protons of the same energy and angle came back with 173 and 1617 hits, and the one
/// with 1617 had its pattern radius 11 % low while the one with 173 was exact. Same physics cannot
/// give a ten-fold difference in hit count, so the extra hits are either not the proton's or not
/// real. Every digitised hit carries AtHit::MCSimPoint (A, Z, trackID), so this can be answered
/// rather than argued.

#include <algorithm>
#include <map>
#include <vector>

static double pu_q(std::vector<double> v, double p)
{
   if (v.empty()) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_purity_C14(TString simFile = "/mnt/f/a1954_C14dp_hf/sims_b285/gs_s8001_sim.root",
                   TString recoFile = "/mnt/f/a1954_C14dp_matfx/ab_reco.root", Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   const double u = 931.49401, mp = 1.007825 * u;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   TFile *fr = TFile::Open(recoFile);
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   // split backward tracks by hit count and report what the hits are made of
   const int NB = 3;
   const double hLo[NB] = {0, 400, 1000}, hHi[NB] = {400, 1000, 1e9};
   std::vector<double> pur[NB], rat[NB], frP[NB], frC[NB], frOther[NB], nh[NB];

   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pm = std::sqrt(px * px + py * py + pz * pz);
         if (pm <= 0) continue;
         keT = std::sqrt(pm * pm + mp * mp) - mp;
         thT = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (keT <= 0 || thT < 92 || thT > 140) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      const AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = &t; }
      if (!best || nb < 30) continue;

      long nProton = 0, nCarbon = 0, nOther = 0, nNoTruth = 0;
      std::map<int, long> byTrack;
      for (const auto &h : best->GetHitArray()) {
         const auto &m = h->GetMCSimPointArray();
         if (m.empty()) { ++nNoTruth; continue; }
         ++byTrack[m[0].trackID];
         if (m[0].Z == 1 && m[0].A == 1) ++nProton;
         else if (m[0].Z == 6) ++nCarbon;
         else ++nOther;
      }
      long tot = nProton + nCarbon + nOther + nNoTruth;
      if (tot < 30) continue;
      long biggest = 0;
      for (auto &kv : byTrack) biggest = std::max(biggest, kv.second);

      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField);
      int b = -1;
      for (int j = 0; j < NB; ++j) if ((double)tot >= hLo[j] && (double)tot < hHi[j]) b = j;
      if (b < 0) continue;
      pur[b].push_back((double)biggest / tot);
      rat[b].push_back(rTrue > 0 ? best->GetGeoRadius() / rTrue : NAN);
      frP[b].push_back(100.0 * nProton / tot);
      frC[b].push_back(100.0 * nCarbon / tot);
      frOther[b].push_back(100.0 * (nOther + nNoTruth) / tot);
      nh[b].push_back((double)tot);
   }

   printf("\n  BACKWARD tracks (true theta 92-140 deg), split by how many hits the pattern gave them\n");
   printf("  %-14s %6s %8s | %9s %9s %9s | %10s %10s\n", "hits in track", "n", "<hits>", "%% proton", "%% carbon",
          "%% other", "purity", "R_geo/R_true");
   for (int b = 0; b < NB; ++b) {
      if (pur[b].size() < 5) { printf("  %5.0f-%-8.0f %6zu   (too few)\n", hLo[b], hHi[b], pur[b].size()); continue; }
      printf("  %5.0f-%-8.0f %6zu %8.0f | %9.1f %9.1f %9.1f | %10.3f %10.3f\n", hLo[b], hHi[b], pur[b].size(),
             pu_q(nh[b], .5), pu_q(frP[b], .5), pu_q(frC[b], .5), pu_q(frOther[b], .5), pu_q(pur[b], .5),
             pu_q(rat[b], .5));
   }
   printf("\n  'purity' is the largest single MC track's share of the cluster. If the fat tracks are pure\n"
          "  proton then the extra hits are the proton's own and the problem is the trajectory model;\n"
          "  if they are not, the clustering is merging something in.\n\n");
   fs->Close(); fr->Close();
}
