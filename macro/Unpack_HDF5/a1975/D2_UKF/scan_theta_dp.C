/// @file scan_theta_dp.C
/// @brief Replay AtPRA's theta determination over many tracks: how much of each track actually
///        reaches the arc-vs-z fit, how steep the line is, and how close the forward/backward
///        sign is to a coin flip.
///
/// Writes TWO products:
///   csvOut  -- one row per track (see the header below), for the distributions
///   listOut -- a GUI list file (<run> <entry> <tid>) so gui_fit_dp.C can step through the same
///              tracks with Next/Prev
///
///   root -b -q 'scan_theta_dp.C("0020,0016",300)'
#include "pra_theta.C"
#include <fstream>

void scan_theta_dp(TString runsCSV = "0020", int minHits = 300, int maxPerRun = 400,
                   TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                   TString csvOut = "/tmp/claude-1000/-home-yassid/d509a965-524d-4548-9401-35c7b60d6646/scratchpad/theta_scan.csv",
                   TString listOut = "/mnt/f/a1975/caches/.explorer/longspiral_list.txt")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ofstream o(csvOut.Data()), l(listOut.Data());
   o << "run,entry,tid,nhits,ncl,stored_theta,replay_theta,sign,narc,ninlier,fracHitsUsed,"
        "fracInlier,slope,radius,arcSpan\n";
   TObjArray *toks = runsCSV.Tokenize(",");
   long total = 0;
   for (int r = 0; r < toks->GetEntries(); ++r) {
      TString num = ((TObjString *)toks->At(r))->GetString();
      TString fn = gfDir + "run_" + num + "_multifit_genfitter_p.root";
      auto *f = TFile::Open(fn);
      if (!f || f->IsZombie()) { printf("  skip %s\n", fn.Data()); continue; }
      auto *t = (TTree *)f->Get("cbmsim");
      TClonesArray *a = nullptr; t->SetBranchAddress("AtTrackingEvent", &a);
      int n = 0;
      for (Long64_t i = 0; i < t->GetEntries() && n < maxPerRun; ++i) {
         t->GetEntry(i);
         auto *ev = (AtTrackingEvent *)a->At(0); if (!ev) continue;
         for (auto &tr : ev->GetTrackArray()) {
            const int nh = tr.GetHitArray().size();
            if (nh < minHits) continue;
            auto p = praThetaHits(tr);
            if (!p.ok) continue;
            const int nin = std::count(p.inlier.begin(), p.inlier.end(), 1);
            const double slope = (std::fabs(p.dirX) > 1e-9) ? p.dirY / p.dirX : 1e9;
            const double alo = *std::min_element(p.arc.begin(), p.arc.end());
            const double ahi = *std::max_element(p.arc.begin(), p.arc.end());
            o << num << "," << i << "," << tr.GetTrackID() << "," << nh << ","
              << tr.GetHitClusterArray()->size() << ","
              << tr.GetGeoTheta() * TMath::RadToDeg() << "," << p.angleDeg << "," << p.sign << ","
              << p.arc.size() << "," << nin << "," << (double)p.arc.size() / nh << ","
              << (double)nin / p.arc.size() << "," << slope << "," << p.radius << ","
              << (ahi - alo) << "\n";
            l << num << " " << i << " " << tr.GetTrackID() << "\n";
            ++n; ++total;
            if (n >= maxPerRun) break;
         }
      }
      f->Close();
      printf("  run_%s: %d tracks\n", num.Data(), n);
   }
   o.close(); l.close();
   printf("\n  %ld tracks -> %s\n  GUI list -> %s\n", total, csvOut.Data(), listOut.Data());
}
