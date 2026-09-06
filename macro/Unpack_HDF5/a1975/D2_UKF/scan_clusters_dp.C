/// @file scan_clusters_dp.C
/// @brief Fit ONE track repeatedly using different fractions of its clusters.
///
/// Two scans, because "fewer clusters" can mean two different things and they have different
/// causes:
///   PREFIX      the first X% along the stored order -- less of the track, and for a stopping
///               spiral that means less of the DECELERATED end. Sensitive to energy loss.
///   DOWNSAMPLE  every k-th cluster, so the full length is kept and only the DENSITY drops.
///               Sensitive to the number of measurements alone.
/// If the fit recovers under downsampling, the count is the problem. If it recovers only under
/// the prefix, the problem is the stopping tail, not the count.
///
///   root -b -q 'scan_clusters_dp.C("run_0020_multifit",9789,0,0)'   // 0 = z-sort (production)
void scan_clusters_dp(TString runTag = "run_0020_multifit", Long64_t entry = 9789, int tid = 0,
                      int orderMode = 0,
                      TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                      TString geoName = "ATTPC_D300torr_v2_geomanager.root",
                      TString eLossTable = "proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf = TFile::Open(dir + "/geometry/" + geoName); gf->Get("FAIRGeom"); }
   TString elossPath = dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable;

   TFile *f = TFile::Open(gfDir + runTag + "_genfitter_p.root");
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr; t->SetBranchAddress("AtTrackingEvent", &te);
   t->GetEntry(entry);
   auto *ev = (AtTrackingEvent *)te->At(0);
   auto trks = ev->GetTrackArray();
   AtTrack *src = nullptr;
   for (auto &tr : trks) if (tr.GetTrackID() == tid) src = &tr;
   if (!src) { printf("track not found\n"); return; }
   std::vector<AtHitCluster> all = *src->GetHitClusterArray();
   const int N = all.size();
   printf("\n  %s entry %lld track %d : %d clusters, ordering %s\n\n", runTag.Data(), entry, tid, N,
          orderMode == 0 ? "z-sort (production)" : "cluster order");

   auto runFit = [&](const std::vector<AtHitCluster> &use, double &ke, double &th, double &c2n,
                     double &ndfOut, double &rms, int &npts) {
      auto fitter = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 1.00782503207, 1,
                                                            elossPath.Data(), kFALSE, 2, 5);
      fitter->SetCatimaMaterial(kTRUE, kTRUE);
      fitter->SetCatimaELoss(kTRUE, kFALSE);
      fitter->SetELossHybrid(kTRUE, 6.61e-5);
      fitter->SetZPadPlane(1000.0);
      fitter->SetMeasSigma(4.0);
      fitter->SetSeedFromSpyral(kFALSE);
      fitter->SetMatEffectsFallback(kFALSE);
      fitter->SetThetaWindow(10.0, 170.0);
      fitter->SetBackwardSeedFix(kTRUE);
      fitter->SetBackExtrapToAxis(kTRUE);
      fitter->SetUseClusterOrder(orderMode != 0);
      fitter->Init();
      AtTrack track = *src;
      *track.GetHitClusterArray() = use;
      AtPatternEvent pe; pe.AddTrack(track);
      AtTrackingEvent tev;
      fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
      if (tev.GetFittedTracks().empty()) { ke = th = c2n = ndfOut = rms = -1; npts = 0; return false; }
      auto *ft = tev.GetFittedTracks().front().get();
      auto &k = ft->GetKinematicsXtr();
      double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
      auto &sp = ft->GetSmoothedPositions();
      ke = k.kineticEnergy; th = k.theta * TMath::RadToDeg();
      c2n = ndf > 0 ? c2 / ndf : 1e9; ndfOut = ndf; npts = sp.size();
      double s2 = 0;
      for (auto &c : use) { auto q = c.GetPosition(); double bd = 1e30;
         for (auto &m : sp) { double dz = q.Z() - (1000.0 - m.Z());
            double d = (q.X()-m.X())*(q.X()-m.X()) + (q.Y()-m.Y())*(q.Y()-m.Y()) + dz*dz;
            if (d < bd) bd = d; }
         s2 += bd; }
      rms = std::sqrt(s2 / use.size());
      return true;
   };

   // THIRD SCAN: a charge threshold. The tail that breaks this track is LOW-charge junk
   // (mean 421 against 3114 over the first 80%, jumping 30 mm in z per step against 1.7 mm),
   // so a charge cut should remove it without discarding any real track -- unlike a prefix cut,
   // which also throws away the genuine high-charge Bragg region.
   printf("  \033[1mCHARGE CUT  (drop clusters below q, full track otherwise)\033[0m\n");
   printf("  %7s %7s | %9s %8s %10s %7s %8s\n", "q_min", "nclus", "KE[MeV]", "theta", "chi2/ndf", "ndf", "rms[mm]");
   for (double qmin : {0.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0}) {
      std::vector<AtHitCluster> use;
      for (auto &c : all) if (c.GetCharge() >= qmin) use.push_back(c);
      if ((int)use.size() < 4) continue;
      double ke, th, c2n, ndf, rms; int npts;
      if (!runFit(use, ke, th, c2n, ndf, rms, npts))
         printf("  %7.0f %7zu | %9s %8s %10s %7s %8s\n", qmin, use.size(), "-", "-", "NO FIT", "-", "-");
      else
         printf("  %7.0f %7zu | %9.3f %8.2f %10.2f %7.0f %8.1f\n", qmin, use.size(), ke, th, c2n, ndf, rms);
   }
   printf("\n");

   const int fracs[] = {5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};
   for (int scan = 0; scan < 2; ++scan) {
      printf("  \033[1m%s\033[0m\n", scan == 0 ? "PREFIX  (first X% of the stored order)"
                                               : "DOWNSAMPLE  (every k-th cluster, full length)");
      printf("  %5s %7s | %9s %8s %10s %7s %8s\n", "X%", "nclus", "KE[MeV]", "theta", "chi2/ndf", "ndf", "rms[mm]");
      for (int fr : fracs) {
         std::vector<AtHitCluster> use;
         if (scan == 0) {
            int m = std::max(4, (int)std::lround(N * fr / 100.0));
            for (int i = 0; i < m && i < N; ++i) use.push_back(all[i]);
         } else {
            int stride = std::max(1, (int)std::lround(100.0 / fr));
            for (int i = 0; i < N; i += stride) use.push_back(all[i]);
         }
         if ((int)use.size() < 4) continue;
         double ke, th, c2n, ndf, rms; int npts;
         bool ok = runFit(use, ke, th, c2n, ndf, rms, npts);
         if (!ok) { printf("  %5d %7zu | %9s %8s %10s %7s %8s\n", fr, use.size(), "-", "-", "NO FIT", "-", "-"); continue; }
         printf("  %5d %7zu | %9.3f %8.2f %10.2f %7.0f %8.1f\n", fr, use.size(), ke, th, c2n, ndf, rms);
      }
      printf("\n");
   }
}
