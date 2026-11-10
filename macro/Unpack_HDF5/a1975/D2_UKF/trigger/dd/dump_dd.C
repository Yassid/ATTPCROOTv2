/// @file dump_dd.C
/// @brief 16C(d,d) elastic scale test: one row per deuteron-hypothesis genfit track (r280 geometry, CATIMA),
///        from /mnt/f/a1975/gf_dd_r280/<run>_multifit_genfitter_d.root. 2026-11-10.
///
/// Same matching rule as ex_dt_a1975.C: loop over VALID Spyral results and match trackID to the fit map, so the
/// ~22 % duplicated AtTrackingEvent entries (empty AtPIDEvent) are never counted. NO CUTS: the deuteron PID gate
/// is stored as a flag (gate), the IC amplitude is stored (ic, -1 = no IC file), and the kinematics are left to
/// the analysis (Ex and theta_cm depend on Ebeam(vz), which is applied there).
///
/// Columns: run ev ke th vz vr c2n brho sqrtdedx ncl ntrk gate ic kefit thfit matfx
///   root -l -b -q 'dump_dd.C("run_0016","/mnt/f/a1975/caches/dd_r280/run_0016.root")'

void dump_dd(TString run = "run_0016", TString out = "dd_run_0016.root", TString inDir = "/mnt/f/a1975/gf_dd_r280/",
             TString gateFile = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/"
                                "deuteron_band_d2.json",
             TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gErrorIgnoreLevel = kError;
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   int runNo = TString(run(4, 4)).Atoi();

   std::vector<float> icv;
   TFile *fi = nullptr;
   TString icf = icDir + run + "_IC.root";
   for (int att = 0; att < 5 && !fi; ++att) { // drvfs fails intermittently under parallelism
      if (att)
         gSystem->Sleep(300);
      fi = TFile::Open(icf);
      if (fi && fi->IsZombie()) {
         fi->Close();
         fi = nullptr;
      }
   }
   if (fi) {
      auto *ti = (TTree *)fi->Get("ic");
      float icm = 0;
      ti->SetBranchAddress("icmax", &icm);
      for (Long64_t j = 0; j < ti->GetEntries(); ++j) {
         ti->GetEntry(j);
         icv.push_back(icm);
      }
      fi->Close();
   } else
      printf("*** WARNING %s: NO IC FILE -- ic = -1 for every row\n", run.Data());

   TFile *fu = TFile::Open(inDir + run + "_multifit_genfitter_d.root");
   auto *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   tu->SetBranchAddress("AtTrackingEvent", &te);
   tu->SetBranchAddress("AtPIDEvent", &pe);

   TFile fo(out, "RECREATE");
   TNtuple nt("d", "16C(d,d) deuteron-hypothesis tracks",
              "run:ev:ke:th:vz:vr:c2n:brho:sqrtdedx:ncl:ntrk:gate:ic:kefit:thfit:matfx");
   long nval = 0, nmatch = 0;
   for (Long64_t i = 0; i < tu->GetEntries(); ++i) {
      tu->GetEntry(i);
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
         ++nval;
         auto it = fmap.find(sr.trackID);
         if (it == fmap.end())
            continue;
         ++nmatch;
         auto *ft = it->second;
         auto &k = ft->GetKinematicsXtr();
         auto &kf = ft->GetKinematics();
         double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
         auto v = ft->GetVertex();
         float row[16] = {(float)runNo,
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
                          (float)(pid.IsInside(sr.sqrtdEdx, sr.brho) ? 1 : 0),
                          (i < (Long64_t)icv.size()) ? icv[i] : -1.f,
                          (float)kf.kineticEnergy,
                          (float)(kf.theta * TMath::RadToDeg()),
                          (float)ft->GetTrackMetadata()->GetMatEffects()};
         nt.Fill(row);
      }
   }
   fo.cd();
   nt.Write();
   fo.Close();
   printf("%s: entries %lld, IC %zu, valid Spyral %ld, matched %ld, rows %lld\n", run.Data(), tu->GetEntries(),
          icv.size(), nval, nmatch, nmatch);
}
