/// @file dump_dd_sim.C
/// @brief Elastic 16C(d,d) simulation (r280) -> two trees for the scale-test acceptance (2026-11-10).
///   "g": every generated event, truth of the recoil deuteron: ev tke tth tvz (tvz in mm, TRUE frame)
///   "d": reconstructed tracks with EXACTLY the columns and matching rule of dump_dd.C (valid Spyral result matched
///        to a fit, so repeated genfit entries are not counted), plus the event truth: ... tke tth tvz
/// The data selection is applied in the analysis on "d"; the acceptance is selected "d" / "g" per cell.
///   root -l -b -q 'dump_dd_sim.C("/mnt/f/a1975/sim_rebuild/r280/dd_el/el_s6001_genfitter_d.root",
///                               "/mnt/f/a1975/sim_rebuild/r280/dd_el/el_s6001_sim.root","out.root")'
void dump_dd_sim(TString gfFile, TString simFile, TString outFile)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   const double md = 2.01355321 * 931.49401;
   TFile fg(gfFile), fs(simFile);
   auto *tg = (TTree *)fg.Get("cbmsim");
   auto *ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr, *mc = nullptr;
   tg->SetBranchAddress("AtTrackingEvent", &te);
   tg->SetBranchAddress("AtPIDEvent", &pe);
   ts->SetBranchStatus("*", 0);
   ts->SetBranchStatus("MCTrack*", 1);
   ts->SetBranchAddress("MCTrack", &mc);

   TFile fo(outFile, "RECREATE");
   TNtuple g("g", "generated recoil deuterons", "ev:tke:tth:tvz");
   TNtuple d("d", "reconstructed tracks (dump_dd.C columns) + truth",
             "run:ev:ke:th:vz:vr:c2n:brho:sqrtdedx:ncl:ntrk:gate:ic:kefit:thfit:matfx:tke:tth:tvz");
   const Long64_t N = ts->GetEntries();
   std::vector<std::array<float, 3>> truth(N, {-1.f, -1.f, -1.f});
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (!m || m->GetPdgCode() != 1000010020 || m->GetMotherId() != -1)
            continue;
         const double p = m->GetP() * 1000.;
         truth[i] = {(float)(std::sqrt(p * p + md * md) - md), (float)(std::acos(m->GetPz() / m->GetP()) * TMath::RadToDeg()),
                     (float)(m->GetStartZ() * 10.)};
         break;
      }
      g.Fill((float)i, truth[i][0], truth[i][1], truth[i][2]);
   }
   long rows = 0;
   for (Long64_t i = 0; i < std::min(tg->GetEntries(), N); ++i) {
      tg->GetEntry(i);
      if (pe->GetEntries() == 0 || te->GetEntries() == 0)
         continue;
      auto *pidev = (AtPIDEvent *)pe->At(0);
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!pidev || !ev)
         continue;
      std::map<int, AtFittedTrack *> fmap;
      for (auto &ft : ev->GetFittedTracks())
         if (ft)
            fmap[ft->GetTrackID()] = ft.get();
      for (auto &sr : pidev->GetSpyral()) {
         if (!sr.valid)
            continue;
         auto it = fmap.find(sr.trackID);
         if (it == fmap.end())
            continue;
         auto *ft = it->second;
         auto &k = ft->GetKinematicsXtr();
         auto &kf = ft->GetKinematics();
         double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
         auto v = ft->GetVertex();
         float row[19] = {0.f,
                          (float)i,
                          (float)k.kineticEnergy,
                          (float)(k.theta * TMath::RadToDeg()),
                          (float)v.Z(),
                          (float)std::sqrt(v.X() * v.X() + v.Y() * v.Y()),
                          (float)(ndf > 0 ? chi2 / ndf : 1e9),
                          (float)sr.brho,
                          (float)sr.sqrtdEdx,
                          (float)sr.nClusters,
                          (float)fmap.size(),
                          0.f,
                          1150.f, // carbon IC, so the data IC gate passes
                          (float)kf.kineticEnergy,
                          (float)(kf.theta * TMath::RadToDeg()),
                          (float)ft->GetTrackMetadata()->GetMatEffects(),
                          truth[i][0],
                          truth[i][1],
                          truth[i][2]};
         d.Fill(row);
         ++rows;
      }
   }
   fo.cd();
   g.Write();
   d.Write();
   fo.Close();
   printf("%s: generated %lld, reco rows %ld\n", outFile.Data(), N, rows);
}
