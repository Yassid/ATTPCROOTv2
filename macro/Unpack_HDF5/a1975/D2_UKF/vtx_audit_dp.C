/// @file vtx_audit_dp.C
/// @brief How much of the ACCEPTED (d,p) sample has a fitted vertex that is not on the beam axis?
///
/// The production cut is chi2/ndf, and chi2/ndf falls mechanically as clusters are removed, so
/// after truncation it stops being a failure detector: a fit to 20 of 391 clusters can sit at 0.38
/// and pass everything.  The back-extrapolated vertex cannot be faked that way -- the reaction
/// happened on the beam axis whatever subset of the track was fitted -- so this counts, among the
/// tracks the analysis KEEPS, how many put their vertex somewhere the beam never was.
///
/// Run it on two productions to see whether truncation created the problem or merely failed to
/// fix it.  Whole runs, not a curated list, because a list built from a failure metric cannot
/// answer "what fraction".
#include <cmath>
#include <vector>
void vtx_audit_dp(TString gfDir = "/mnt/f/a1975/gf_dp_find805/",
                  TString runs = "run_0016_multifit,run_0020_multifit",
                  double c2Cut = 5.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TObjArray *ra = runs.Tokenize(",");
   std::vector<double> vAcc, vAll; std::vector<int> nptAcc;
   long nFit = 0, nAcc = 0;
   std::vector<double> prevSig;
   for (int i = 0; i < ra->GetEntries(); ++i) {
      TString rt = ((TObjString *)ra->At(i))->GetString();
      TFile *f = TFile::Open(gfDir + rt + "_genfitter_p.root");
      if (!f || f->IsZombie()) { printf("  skip %s\n", rt.Data()); continue; }
      auto *t = (TTree *)f->Get("cbmsim");
      TClonesArray *a = nullptr; t->SetBranchAddress("AtTrackingEvent", &a);
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         if (a->GetEntries() == 0) continue;
         auto *ev = (AtTrackingEvent *)a->At(0); if (!ev) continue;
         // DUPLICATE-ENTRY GUARD: ~22 % of entries repeat the previous one byte for byte, and
         // without this every duplicated track is counted twice.
         std::vector<double> sig;
         for (auto &tr : ev->GetTrackArray()) sig.push_back(tr.GetHitClusterArray()->size());
         if (!sig.empty() && sig == prevSig) continue;
         prevSig = sig;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft) continue;
            double ndf = ft->GetTrackMetadata()->GetNdf();
            if (!(ndf > 0)) continue;
            double c2 = ft->GetTrackMetadata()->GetChi2() / ndf;
            auto v = ft->GetVertex();
            double r = std::hypot(v.X(), v.Y());
            ++nFit; vAll.push_back(r);
            if (c2 < c2Cut) { ++nAcc; vAcc.push_back(r); nptAcc.push_back((int)ft->GetSmoothedPositions().size()); }
         }
      }
      f->Close();
   }
   auto q = [](std::vector<double> v, double f) { if (v.empty()) return -1.0; std::sort(v.begin(), v.end());
        return v[std::min(v.size()-1, (size_t)(f*v.size()))]; };
   auto frac = [](const std::vector<double> &v, double thr) { if (v.empty()) return -1.0; long n = 0;
        for (double x : v) if (x > thr) ++n; return 100.0*n/v.size(); };
   printf("\n  %s\n  %ld fitted tracks, %ld pass chi2/ndf < %.1f (%.1f %%)\n",
          gfDir.Data(), nFit, nAcc, c2Cut, nFit ? 100.0*nAcc/nFit : 0.0);
   printf("\n  fitted vertex |xy| [mm]      median      p90      p99   >10mm   >20mm   >50mm\n");
   printf("  all fitted            %11.1f %8.1f %8.1f %6.1f%% %6.1f%% %6.1f%%\n",
          q(vAll,0.5), q(vAll,0.90), q(vAll,0.99), frac(vAll,10), frac(vAll,20), frac(vAll,50));
   printf("  ACCEPTED (chi2 < %.1f) %11.1f %8.1f %8.1f %6.1f%% %6.1f%% %6.1f%%\n", c2Cut,
          q(vAcc,0.5), q(vAcc,0.90), q(vAcc,0.99), frac(vAcc,10), frac(vAcc,20), frac(vAcc,50));
   // is the off-axis subset the SHORT-fit subset?  that is the signature of chi2 being bought
   // with clusters rather than earned
   std::vector<double> nptOk, nptBad;
   for (size_t i = 0; i < vAcc.size(); ++i)
      (vAcc[i] > 20 ? nptBad : nptOk).push_back(nptAcc[i]);
   printf("\n  stored fit points, ACCEPTED tracks    n      p25   median      p75\n");
   printf("    vertex <= 20 mm            %8zu %8.0f %8.0f %8.0f\n", nptOk.size(),
          q(nptOk,0.25), q(nptOk,0.5), q(nptOk,0.75));
   printf("    vertex >  20 mm            %8zu %8.0f %8.0f %8.0f\n", nptBad.size(),
          q(nptBad,0.25), q(nptBad,0.5), q(nptBad,0.75));
}
