/// @file sim_dump_dp.C
/// @brief Per-event dump of a genfit fit of the (d,p) simulation, tree "r" (ke th vz c2 tke tth tvz), read by the D2_UKF/angdist_bound_dp*.C efficiency macros. Rebuilt 2026-11-03: the
///        original producer was lost with a scratchpad; run on the old fit it is the control against gs_s4001_simres_dump.root. Picks the fitted track nearest to truth, as sim_res_dp.C does.
///        Reads only entries with PID results, so the repeated entries of the genfit trees are not counted twice. dedup = false keeps them, to reproduce a file made without that check.
///   root -l -b -q 'sim_dump_dp.C("fit.root","/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root","dump.root")'
void sim_dump_dp(TString gfFile, TString simFile, TString outFile, bool dedup = true)
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   TFile fg(gfFile), fs(simFile);
   auto *tg = (TTree *)fg.Get("cbmsim"); auto *ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr, *mc = nullptr;
   tg->SetBranchAddress("AtTrackingEvent", &te); tg->SetBranchAddress("AtPIDEvent", &pe);
   ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc);
   TFile fo(outFile, "RECREATE");
   TTree r("r", "per-event sim fit dump");
   double ke, th, vz, c2, tke, tth, tvz;
   r.Branch("ke", &ke); r.Branch("th", &th); r.Branch("vz", &vz); r.Branch("c2", &c2); r.Branch("tke", &tke); r.Branch("tth", &tth); r.Branch("tvz", &tvz);
   for (Long64_t i = 0; i < std::min(tg->GetEntries(), ts->GetEntries()); ++i) {
      tg->GetEntry(i);
      if (!te->GetEntries()) continue;
      if (dedup) {
         auto *pid = pe->GetEntries() ? (AtPIDEvent *)pe->At(0) : nullptr;
         if (!pid || pid->GetSpyral().empty()) continue;
      }
      ts->GetEntry(i);
      tke = -1;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (m && m->GetPdgCode() == 2212) { const double p = m->GetP() * 1000.; tke = std::sqrt(p * p + 938.272 * 938.272) - 938.272; tth = std::acos(m->GetPz() / m->GetP()) * 57.29578; tvz = m->GetStartZ() * 10; break; }
      }
      if (tke < 0) continue;
      double bd = 1e18; bool got = false;
      for (auto &ft : ((AtTrackingEvent *)te->At(0))->GetFittedTracks()) {
         if (!ft) continue;
         auto &k = ft->GetKinematicsXtr();
         if (k.kineticEnergy <= 0) continue;
         const double t = k.theta * 57.29578, d = TMath::Sq(k.kineticEnergy - tke) + TMath::Sq((t - tth) / 5.0);
         if (d >= bd) continue;
         auto &md = ft->GetTrackMetadata();
         const double ndf = md ? md->GetNdf() : -1;
         bd = d; got = true; ke = k.kineticEnergy; th = t; vz = ft->GetVertex().Z(); c2 = ndf > 0 ? md->GetChi2() / ndf : 1e9;
      }
      if (got) r.Fill();
   }
   fo.Write();
   printf("%s: %lld rows\n", outFile.Data(), r.GetEntries());
}
