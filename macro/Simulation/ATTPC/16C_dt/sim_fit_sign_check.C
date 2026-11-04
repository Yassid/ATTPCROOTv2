/// @file sim_fit_sign_check.C
/// @brief Compares genfit output on the a1975 D2 simulation with MC truth, one line per file: tracks fitted, median and spread of ln(KE/true) at the vertex,
///        theta - true and (180 - theta) - true, vz - true z. Truth-matched by species on the hits. True KE comes from |p| and the mass argument, since AtMCTrack::GetMass does not know the ions. Used to settle the field sign the sim needs.
///   root -l -b -q 'sim_fit_sign_check.C("fit.root","reco.root","sim.root",1000010030,3,1,2808.921)'
namespace sfc {
double pct(std::vector<double> &v, double f) { if (v.empty()) return 0; std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, size_t(f * v.size()))]; }
void line(const char *n, std::vector<double> &v) { printf("   %-22s median %+8.3f   p16 %+8.3f   p84 %+8.3f\n", n, pct(v, .5), pct(v, .16), pct(v, .84)); }
}
void sim_fit_sign_check(TString fit, TString reco, TString sim, int pdg, int A, int Z, double mass)
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   TFile ff(fit), fr(reco), fs(sim);
   auto *tf = (TTree *)ff.Get("cbmsim"); auto *tr = (TTree *)fr.Get("cbmsim"); auto *ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr, *pa = nullptr, *mc = nullptr;
   tf->SetBranchAddress("AtTrackingEvent", &te); tf->SetBranchAddress("AtPIDEvent", &pe);
   tr->SetBranchStatus("*", 0); tr->SetBranchStatus("AtPatternEvent*", 1); tr->SetBranchAddress("AtPatternEvent", &pa);
   ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc);
   std::vector<double> dk, dt, df, dz, chi;
   long nTrue = 0, nFit = 0, nDir = 0;
   for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
      tf->GetEntry(i);
      if (!pe->GetEntries() || !te->GetEntries()) continue;
      auto *pid = (AtPIDEvent *)pe->At(0);
      if (!pid || pid->GetSpyral().empty()) continue;
      tr->GetEntry(i); ts->GetEntry(i);
      double tke = -1, tth = 0, tz = 0;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (m && m->GetPdgCode() == pdg) { const double pm = m->GetP() * 1000.; tke = std::sqrt(pm * pm + mass * mass) - mass; tth = std::acos(m->GetPz() / m->GetP()) * 57.29578; tz = m->GetStartZ() * 10; break; }
      }
      if (tke < 0 || !pa->GetEntries()) continue;
      std::set<int> truthIds;
      for (auto &t : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         int n = 0, nt = 0;
         for (auto &h : t.GetHitArray()) { auto &ms = h->GetMCSimPointArray(); if (ms.empty()) continue; ++n; if (ms[0].A == A && ms[0].Z == Z) ++nt; }
         if (n && nt > n / 2) truthIds.insert(t.GetTrackID());
      }
      if (truthIds.empty()) continue;
      ++nTrue;
      bool got = false;
      for (auto &ft : ((AtTrackingEvent *)te->At(0))->GetFittedTracks()) {
         if (!ft || !truthIds.count(ft->GetTrackID())) continue;
         auto &k = ft->GetKinematicsXtr();
         const double ndf = ft->GetTrackMetadata()->GetNdf();
         if (k.kineticEnergy <= 0 || ndf <= 0) continue;
         got = true;
         const double th = k.theta * 57.29578;
         dk.push_back(std::log(k.kineticEnergy / tke)); dt.push_back(th - tth); df.push_back(180 - th - tth); dz.push_back(ft->GetVertex().Z() - tz); chi.push_back(ft->GetTrackMetadata()->GetChi2() / ndf);
         nDir += std::abs(th - tth) < std::abs(180 - th - tth);
         break;
      }
      nFit += got;
   }
   printf("%s\n   events with a truth-matched pattern track %ld, fitted %ld (%.1f %%), angle closer to theta than to 180-theta %.1f %%\n", fit.Data(), nTrue, nFit, 100. * nFit / std::max(1L, nTrue), 100. * nDir / std::max(1L, nFit));
   sfc::line("ln(KE/true)", dk); sfc::line("theta - true [deg]", dt); sfc::line("180-theta - true [deg]", df); sfc::line("vz - z_true [mm]", dz); sfc::line("chi2/ndf", chi);
}
