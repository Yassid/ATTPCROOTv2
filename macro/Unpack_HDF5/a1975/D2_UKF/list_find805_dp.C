/// @file list_find805_dp.C
/// @brief Stratified track list from gf_dp_find805 for stepping through gui_fit_dp.C (prev/next).
///
/// Four classes, nPer tracks each, spread evenly over the runs given (not the first nPer found):
///   A  full length kept, chi2/ndf < c2Max          control: "find" should leave these alone
///   B  truncated to 25-75 %, passed                 the typical FIND case
///   C  truncated below 15 %, passed                 suspicious: few points can fake a low chi2
///   D  chi2/ndf >= c2Max                            nothing passed, minimum-chi2 rung kept
/// Kept fraction = stored fit points / clusters (FIND writes no findPct), a LOWER BOUND.
/// Duplicate tree entries are skipped by requiring a valid Spyral record (dump_tracks_dp.C).
///
/// Writes <out> as "run entry tid" (gui_fit_dp.C's list format, grouped A,B,C,D) and <out>.csv
/// with the class and numbers, so a frame in the GUI can be looked up.
///
///   root -l -b -q 'list_find805_dp.C("0016,0020,0026,0034", 10, "find805_list.txt")'
#include <fstream>
#include <map>
#include <vector>
void list_find805_dp(TString runsCSV = "0016,0020,0026,0034", int nPer = 10,
                     TString out = "find805_list.txt", int minCl = 60, double c2Max = 0.1,
                     TString gfDir = "/mnt/f/a1975/gf_dp_find805/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   struct Row { TString run; Long64_t e; int tid, ncl, np; double c2, ke, th; };
   std::map<char, std::vector<Row>> cls;
   TObjArray *toks = runsCSV.Tokenize(",");
   for (int r = 0; r < toks->GetEntries(); ++r) {
      TString run = ((TObjString*)toks->At(r))->GetString();
      auto *f = TFile::Open(gfDir + "run_" + run + "_multifit_genfitter_p.root");
      if (!f || f->IsZombie()) { printf("  skip run %s\n", run.Data()); continue; }
      auto *t = (TTree*)f->Get("cbmsim");
      TClonesArray *a = nullptr, *p = nullptr;
      t->SetBranchAddress("AtTrackingEvent", &a); t->SetBranchAddress("AtPIDEvent", &p);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!p || p->GetEntries() == 0 || !a || a->GetEntries() == 0) continue;
         auto *pid = (AtPIDEvent*)p->At(0); auto *ev = (AtTrackingEvent*)a->At(0);
         if (!pid || !ev) continue;
         std::map<int, int> ncl;
         for (auto &tr : ev->GetTrackArray()) ncl[tr.GetTrackID()] = tr.GetHitClusterArray()->size();
         std::map<int, AtFittedTrack*> fm;
         for (auto &ft : ev->GetFittedTracks()) if (ft) fm[ft->GetTrackID()] = ft.get();
         for (auto &sr : pid->GetSpyral()) {
            if (!sr.valid) continue;
            auto it = fm.find(sr.trackID); if (it == fm.end() || !ncl.count(sr.trackID)) continue;
            auto *ft = it->second;
            int n = ncl[sr.trackID], np = ft->GetSmoothedPositions().size();
            double ndf = ft->GetTrackMetadata()->GetNdf();
            if (n < minCl || np < 3 || ndf <= 0) continue;
            double c2 = ft->GetTrackMetadata()->GetChi2() / ndf, fr = (double)np / n;
            auto &k = ft->GetKinematicsXtr();
            Row row{run, i, sr.trackID, n, np, c2, k.kineticEnergy, k.theta * TMath::RadToDeg()};
            char c = 0;
            if (c2 >= c2Max) c = 'D';
            else if (fr > 0.93) c = 'A';
            else if (fr >= 0.25 && fr <= 0.75) c = 'B';
            else if (fr < 0.15) c = 'C';
            if (c) cls[c].push_back(row);
         }
      }
      f->Close();
   }
   std::ofstream fl(out.Data()), fc((out + ".csv").Data());
   fc << "idx,class,run,entry,tid,ncl,fitpts,kept_pct,chi2ndf,ke,theta\n";
   int idx = 0;
   for (char c : {'A', 'B', 'C', 'D'}) {
      auto &v = cls[c];
      const size_t stride = v.size() > (size_t)nPer ? v.size() / nPer : 1;
      int took = 0;
      for (size_t k = 0; k < v.size() && took < nPer; k += stride, ++took) {
         auto &r = v[k];
         fl << r.run << " " << r.e << " " << r.tid << "\n";
         fc << TString::Format("%d,%c,%s,%lld,%d,%d,%d,%.0f,%.4f,%.3f,%.1f\n", ++idx, c, r.run.Data(), r.e,
                               r.tid, r.ncl, r.np, 100.0 * r.np / r.ncl, r.c2, r.ke, r.th).Data();
      }
      printf("  class %c: %6zu tracks, took %d\n", c, v.size(), took);
   }
   printf("  wrote %s and %s.csv\n", out.Data(), out.Data());
}
