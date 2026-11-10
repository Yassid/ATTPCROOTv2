/// @file dump_dd_valid.C
/// @brief Spyral validity per genfit track, data or sim (2026-11-10). Same matching rule as dump_dd.C (loop the
///        Spyral results, match trackID to the fit map: the repeated genfit entries carry no AtPIDEvent) but every
///        Spyral result is kept, with its valid flag. Used to compare the Spyral efficiency of data and sim vs angle.
/// Tree s: every Spyral result, run ev brho polar(deg) valid hasfit ic ncl (fit efficiency on the locus, fit-free).
/// Tree v columns: run ev ke th vz vr c2n valid ic   (kinematics from the genfit fit, which exists also for invalid results;
///   ic = -1 when no IC file is given or the entry has none)
///   root -l -b -q 'dump_dd_valid.C("/mnt/f/a1975/gf_dd_r280/run_0016_multifit_genfitter_d.root","out.root",16,
///                                  "/mnt/f/a1975/ic_d2/run_0016_IC.root")'
void dump_dd_valid(TString in, TString out, int runNo = 0, TString icf = "")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kError;
   std::vector<float> icv;
   if (icf != "") {
      TFile *fi = TFile::Open(icf);
      if (fi && !fi->IsZombie()) {
         auto *ti = (TTree *)fi->Get("ic");
         float icm = 0;
         ti->SetBranchAddress("icmax", &icm);
         for (Long64_t j = 0; j < ti->GetEntries(); ++j) {
            ti->GetEntry(j);
            icv.push_back(icm);
         }
         fi->Close();
      } else
         printf("*** WARNING: no IC file %s\n", icf.Data());
   }
   TFile *fu = TFile::Open(in);
   auto *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   tu->SetBranchAddress("AtTrackingEvent", &te);
   tu->SetBranchAddress("AtPIDEvent", &pe);
   TFile fo(out, "RECREATE");
   TNtuple nt("v", "Spyral validity per fitted track", "run:ev:ke:th:vz:vr:c2n:valid:ic");
   TNtuple ns("s", "every Spyral result: does it have a genfit fit", "run:ev:brho:polar:valid:hasfit:ic:ncl");
   long nall = 0, nval = 0, nfit = 0;
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
         ++nall;
         nval += sr.valid;
         auto it = fmap.find(sr.trackID);
         ns.Fill((float)runNo, (float)i, (float)sr.brho, (float)(sr.polar * TMath::RadToDeg()), (float)(sr.valid ? 1 : 0),
                 (float)(it != fmap.end() ? 1 : 0), (i < (Long64_t)icv.size()) ? icv[i] : -1.f, (float)sr.nClusters);
         if (it == fmap.end())
            continue;
         ++nfit;
         auto *ft = it->second;
         auto &k = ft->GetKinematicsXtr();
         double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
         auto v = ft->GetVertex();
         nt.Fill((float)runNo, (float)i, (float)k.kineticEnergy, (float)(k.theta * TMath::RadToDeg()), (float)v.Z(),
                 (float)std::sqrt(v.X() * v.X() + v.Y() * v.Y()), (float)(ndf > 0 ? chi2 / ndf : 1e9),
                 (float)(sr.valid ? 1 : 0), (i < (Long64_t)icv.size()) ? icv[i] : -1.f);
      }
   }
   fo.cd();
   nt.Write();
   ns.Write();
   fo.Close();
   printf("%s: Spyral results %ld, valid %ld, with a fit %ld\n", out.Data(), nall, nval, nfit);
}
