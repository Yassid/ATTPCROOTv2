/// @file dump_dt_tracks.C  (d,t) tracks of one run with their reco ENTRY index, for joining to the raw FRIB mesh features.
/// Same loop and gate as ex_dt_a1975.C as used for the adopted dt_kin.root (gf_dt_cateloss, triton_d2_dv1104.json, no chi2/KE/theta cut),
/// so the rows reproduce that cache. Columns: run entry trackID ke theta vz vr chi2ndf ic
void dump_dt_tracks(TString run = "run_0016", TString inDir = "/mnt/f/a1975/gf_dt_cateloss/", TString gateFile = "pid/triton_d2_dv1104.json",
                    TString out = "/mnt/f/a1975/sim_rebuild/trig/dt/tracks_0016.txt", TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gErrorIgnoreLevel = kError; gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   TFile *fu = TFile::Open(inDir + run + "_multifit_genfitter_t.root"); TTree *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr; tu->SetBranchAddress("AtTrackingEvent", &te); tu->SetBranchAddress("AtPIDEvent", &pe);
   std::vector<float> icv; TFile *fi = nullptr; for (int a = 0; a < 5 && !fi; ++a) { if (a) gSystem->Sleep(300); fi = TFile::Open(icDir + run + "_IC.root"); if (fi && fi->IsZombie()) { fi->Close(); fi = nullptr; } }
   if (fi) { TTree *ti = (TTree *)fi->Get("ic"); float icm; ti->SetBranchAddress("icmax", &icm); for (Long64_t j = 0; j < ti->GetEntries(); ++j) { ti->GetEntry(j); icv.push_back(icm); } fi->Close(); }
   int runNo = TString(run(run.Length() - 4, 4)).Atoi(); FILE *o = fopen(out, "w"); long n = 0;
   for (Long64_t i = 0; i < tu->GetEntries(); ++i) {
      tu->GetEntry(i); if (pe->GetEntries() == 0 || te->GetEntries() == 0) continue;
      auto *pidev = (AtPIDEvent *)pe->At(0); auto *ev = (AtTrackingEvent *)te->At(0); if (!pidev || !ev) continue;
      std::map<int, AtFittedTrack *> fmap; for (auto &ft : ev->GetFittedTracks()) if (ft) fmap[ft->GetTrackID()] = ft.get();
      for (auto &sr : pidev->GetSpyral()) {
         if (!sr.valid || !pid.IsInside(sr.sqrtdEdx, sr.brho)) continue; auto it = fmap.find(sr.trackID); if (it == fmap.end()) continue;
         auto *ft = it->second; auto &k = ft->GetKinematicsXtr(); double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
         if (k.kineticEnergy <= 0) continue; auto v = ft->GetVertex();
         fprintf(o, "%d %lld %d %.4f %.4f %.2f %.2f %.4f %.1f\n", runNo, i, sr.trackID, k.kineticEnergy, k.theta * TMath::RadToDeg(), v.Z(), std::hypot(v.X(), v.Y()), ndf > 0 ? c2 / ndf : 1e9, i < (Long64_t)icv.size() ? icv[i] : -1.f); ++n;
      }
   }
   fclose(o); printf("done %s %ld rows ic %zu\n", run.Data(), n, icv.size());
}
