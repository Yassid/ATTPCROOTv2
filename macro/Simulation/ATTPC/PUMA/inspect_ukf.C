/// @file inspect_ukf.C
/// @brief Summarize a multi-hypothesis UKF output: per-event hypothesis
///        breakdown, reduced chi^2, kinetic energy, vertex distance.
///
/// Run: root -b -q 'inspect_ukf.C("data/output_ukf_multi_riemann_1k.root")'

void inspect_ukf(TString file = "data/output_ukf_multi_riemann_1k.root")
{
   auto *f = TFile::Open(file);
   if (!f || f->IsZombie()) {
      std::cerr << "Cannot open " << file << "\n";
      return;
   }

   auto *t = (TTree *)f->Get("cbmsim");
   if (!t) {
      std::cerr << "No cbmsim tree in " << file << "\n";
      return;
   }

   TClonesArray *tracking = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &tracking);

   const Long64_t N = t->GetEntries();
   std::cout << "Events: " << N << "\n";

   std::map<std::string, int> hypoCount;
   int nFitted = 0, nWithMeta = 0, nConverged = 0;
   std::vector<double> chi2red, KE, vtxR, theta;

   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!tracking || tracking->GetEntries() == 0)
         continue;
      auto *te = (AtTrackingEvent *)tracking->At(0);
      const auto &fitted = te->GetFittedTracks();

      for (const auto &ft : fitted) {
         ++nFitted;
         const auto &pi = ft->GetParticleInfo();
         std::string name = pi.idPDG.Data();
         if (name.empty()) name = "(unset)";
         hypoCount[name]++;

         const auto &meta = ft->GetTrackMetadata();
         if (meta) {
            ++nWithMeta;
            if (meta->GetFitConverged()) ++nConverged;
            int ndf = meta->GetNdf();
            if (ndf > 0)
               chi2red.push_back(meta->GetChi2() / ndf);
         }

         const auto &kin = ft->GetKinematics();
         if (kin.kineticEnergy > 0) KE.push_back(kin.kineticEnergy);
         if (kin.theta > 0) theta.push_back(kin.theta * 180.0 / M_PI);

         const auto &vtx = ft->GetVertex();
         vtxR.push_back(std::sqrt(vtx.X() * vtx.X() + vtx.Y() * vtx.Y() + vtx.Z() * vtx.Z()));
      }
   }

   auto stats = [](std::vector<double> &v, const char *label, const char *unit) {
      if (v.empty()) {
         std::cout << "  " << label << ": (empty)\n";
         return;
      }
      std::sort(v.begin(), v.end());
      double sum = 0;
      for (double x : v) sum += x;
      double mean = sum / v.size();
      double med = v[v.size() / 2];
      double p10 = v[(size_t)(0.10 * v.size())];
      double p90 = v[(size_t)(0.90 * v.size())];
      std::cout << "  " << label << " [" << unit << "]"
                << "  n=" << v.size() << "  min=" << v.front() << "  p10=" << p10 << "  med=" << med
                << "  mean=" << mean << "  p90=" << p90 << "  max=" << v.back() << "\n";
   };

   std::cout << "\n=== Fitted-track summary ===\n";
   std::cout << "Total fitted tracks  : " << nFitted << "\n";
   std::cout << "With metadata        : " << nWithMeta << "\n";
   std::cout << "Converged            : " << nConverged << "\n";

   std::cout << "\n=== Best-hypothesis breakdown ===\n";
   for (auto &p : hypoCount) {
      double pct = 100.0 * p.second / std::max(nFitted, 1);
      std::cout << "  " << p.first << " : " << p.second << "  (" << Form("%.1f", pct) << "%)\n";
   }

   std::cout << "\n=== Distributions ===\n";
   stats(chi2red, "chi2/ndf  ", "");
   stats(KE,      "KE        ", "MeV");
   stats(theta,   "theta_lab ", "deg");
   stats(vtxR,    "vtx |R|   ", "mm");
}
