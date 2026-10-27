/// @file find_rescue_dp.C
/// @brief Tracks where FIND (gf_dp_find805) genuinely RESCUES a failed full-length fit (gf_dp_cateloss).
///
/// chi2/ndf alone cannot certify a rescue -- it falls mechanically as points are removed. So a
/// candidate must also pass checks chi2 cannot fake:
///   * the full fit FAILS clearly (chi2/ndf >= c2Bad) and find805 PASSES (< c2Good)
///   * find805 keeps a LONG prefix (fit points / clusters >= minKeep; FIND stores no findPct, so this
///     proxy is a lower bound)
///   * the find805 vertex lies within rMax of the beam axis AND no farther than the full fit's
///     (the physical criterion specified for FIND v3, see project_a1975_dp_truncation)
/// and they are RANKED by an independent referee: the KE implied by the Spyral brho of the same
/// track. A real rescue agrees with Spyral where the full-length fit does not.
///
/// Matched on (entry, trackID); duplicate tree entries skipped via the valid-Spyral requirement.
///
///   root -l -b -q 'find_rescue_dp.C("0016,0020,0026,0034", "rescue.csv")'
#include <fstream>
#include <map>
#include <vector>
void find_rescue_dp(TString runsCSV = "0016,0020,0026,0034", TString outCsv = "find_rescue.csv",
                    int minCl = 150, double c2Bad = 1.0, double c2Good = 0.1, double minKeep = 0.5,
                    double rMax = 10.0, TString dirFind = "/mnt/f/a1975/gf_dp_find805/",
                    TString dirFull = "/mnt/f/a1975/gf_dp_cateloss/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   const double mp = 938.272088;
   auto keFromBrho = [&](double brho) { double p = 299.792458 * brho; return std::sqrt(p * p + mp * mp) - mp; };
   struct Row { TString run; Long64_t e; int tid, ncl, np; double c2F, c2B, keF, keB, keS, thF, thB, rF, rB; };
   std::vector<Row> rows;
   long nMatched = 0, nFail = 0, nPass = 0, nLong = 0, nVtx = 0;
   TObjArray *toks = runsCSV.Tokenize(",");
   for (int r = 0; r < toks->GetEntries(); ++r) {
      TString run = ((TObjString*)toks->At(r))->GetString();
      auto *fF = TFile::Open(dirFind + "run_" + run + "_multifit_genfitter_p.root");
      auto *fB = TFile::Open(dirFull + "run_" + run + "_multifit_genfitter_p.root");
      if (!fF || fF->IsZombie() || !fB || fB->IsZombie()) { printf("  skip run %s\n", run.Data()); continue; }
      auto *tF = (TTree*)fF->Get("cbmsim"); auto *tB = (TTree*)fB->Get("cbmsim");
      TClonesArray *aF = nullptr, *pF = nullptr, *aB = nullptr;
      tF->SetBranchAddress("AtTrackingEvent", &aF); tF->SetBranchAddress("AtPIDEvent", &pF);
      tB->SetBranchAddress("AtTrackingEvent", &aB);
      const Long64_t N = std::min(tF->GetEntries(), tB->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         tF->GetEntry(i);
         if (!pF || pF->GetEntries() == 0 || !aF || aF->GetEntries() == 0) continue;
         auto *pid = (AtPIDEvent*)pF->At(0); auto *ev = (AtTrackingEvent*)aF->At(0);
         if (!pid || !ev) continue;
         bool any = false;
         for (auto &sr : pid->GetSpyral()) any |= sr.valid;
         if (!any) continue;
         tB->GetEntry(i);
         auto *evB = (AtTrackingEvent*)aB->At(0);
         if (!evB) continue;
         std::map<int, int> ncl;
         for (auto &tr : ev->GetTrackArray()) ncl[tr.GetTrackID()] = tr.GetHitClusterArray()->size();
         std::map<int, AtFittedTrack*> mF, mB;
         for (auto &f : ev->GetFittedTracks()) if (f) mF[f->GetTrackID()] = f.get();
         for (auto &f : evB->GetFittedTracks()) if (f) mB[f->GetTrackID()] = f.get();
         for (auto &sr : pid->GetSpyral()) {
            if (!sr.valid || !mF.count(sr.trackID) || !mB.count(sr.trackID) || !ncl.count(sr.trackID)) continue;
            auto *f = mF[sr.trackID], *b = mB[sr.trackID];
            int n = ncl[sr.trackID], np = f->GetSmoothedPositions().size();
            if (n < minCl) continue;
            double nF = f->GetTrackMetadata()->GetNdf(), nB = b->GetTrackMetadata()->GetNdf();
            if (nF <= 0 || nB <= 0) continue;
            ++nMatched;
            double c2F = f->GetTrackMetadata()->GetChi2() / nF, c2B = b->GetTrackMetadata()->GetChi2() / nB;
            if (c2B < c2Bad) continue;
            ++nFail;
            if (c2F >= c2Good) continue;
            ++nPass;
            if ((double)np / n < minKeep) continue;
            ++nLong;
            double rF = std::hypot(f->GetVertex().X(), f->GetVertex().Y());
            double rB = std::hypot(b->GetVertex().X(), b->GetVertex().Y());
            if (rF > rMax || rF > rB) continue;
            ++nVtx;
            auto &kF = f->GetKinematicsXtr(); auto &kB = b->GetKinematicsXtr();
            rows.push_back({run, i, sr.trackID, n, np, c2F, c2B, kF.kineticEnergy, kB.kineticEnergy,
                            keFromBrho(sr.brho), kF.theta * TMath::RadToDeg(), kB.theta * TMath::RadToDeg(), rF, rB});
         }
      }
      fF->Close(); fB->Close();
   }
   // referee: relative KE distance to Spyral. Best = find agrees, full does not.
   auto score = [](const Row &r) {
      double dF = std::fabs(r.keF - r.keS) / r.keS, dB = std::fabs(r.keB - r.keS) / r.keS;
      return dB - dF;
   };
   std::sort(rows.begin(), rows.end(), [&](const Row &a, const Row &b) { return score(a) > score(b); });
   std::ofstream o(outCsv.Data());
   o << "run,entry,tid,ncl,fitpts,kept_pct,chi2_full,chi2_find,ke_full,ke_find,ke_spyral,theta_full,theta_find,"
        "rvtx_full,rvtx_find,score\n";
   for (auto &r : rows)
      o << TString::Format("%s,%lld,%d,%d,%d,%.0f,%.3f,%.4f,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%.1f,%.3f\n", r.run.Data(), r.e,
                           r.tid, r.ncl, r.np, 100.0 * r.np / r.ncl, r.c2B, r.c2F, r.keB, r.keF, r.keS, r.thB,
                           r.thF, r.rB, r.rF, score(r)).Data();
   printf("  tracks >= %d clusters matched in both: %ld\n", minCl, nMatched);
   printf("   full fit fails (chi2/ndf >= %.1f)      : %ld\n", c2Bad, nFail);
   printf("   ... find805 passes (< %.2f)            : %ld\n", c2Good, nPass);
   printf("   ... keeping >= %.0f %% of the track      : %ld\n", 100 * minKeep, nLong);
   printf("   ... vertex <= %.0f mm and not worse     : %ld   -> %s\n", rMax, nVtx, outCsv.Data());
}
