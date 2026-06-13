/// @file perf_build.C
/// @brief Build the MASTER per-track diagnostic ntuple for the genfit (p,d)
/// performance study. For every PRA candidate track it joins, by event+trackID:
///   - PRA quality   : nClusters, nHits, geoRadius, geoTheta
///   - GenFit outcome: fitted?, KE, chi2/ndf, converged
///   - UKF outcome   : fitted?, KE, chi2/ndf   (deuteron hyp, reco_pd/_ukf_d)
///   - PID (Spyral)  : valid?, brho, sqrtdEdx, in-deuteron-gate?
/// Output -> perf_master.root (ntuple "m"). Every stage macro reads this, instant.
///
///   root -b -q 'pd/perf/perf_build.C("run_0106,run_0107,run_0108,run_0109,run_0110")'

#include <map>
void perf_build(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110",
                TString recoDir = "/mnt/f/a1975/reco/", TString gfDir = "/mnt/f/a1975/reco_gf/",
                TString ukfDir = "/mnt/f/a1975/reco_pd/", TString gateFile = "pid/deuteron_band.json",
                TString out = "pd/perf/perf_master.root")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile *fo = new TFile(out, "RECREATE");
   TNtuple *m = new TNtuple("m", "perf master",
                            "run:evt:tid:nclus:nhits:radius:gtheta:"
                            "gffit:gfke:gfchi2:gfconv:ukfit:ukke:ukchi2:spyvalid:brho:sqdedx:gated");

   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      float rno = TString(run(4, 4)).Atof();
      TString rc = recoDir + run + "_reco.root", gf = gfDir + run + "_genfit.root";
      TString pf = gfDir + run + "_pid.root", uf = ukfDir + run + "_ukf_d.root";
      if (gSystem->AccessPathName(rc) || gSystem->AccessPathName(gf)) {
         printf("  %s: missing reco/genfit, skip\n", run.Data());
         continue;
      }
      bool hasPid = !gSystem->AccessPathName(pf), hasUkf = !gSystem->AccessPathName(uf);
      TFile *fr = TFile::Open(rc), *fg = TFile::Open(gf);
      TFile *fp = hasPid ? TFile::Open(pf) : nullptr, *fu = hasUkf ? TFile::Open(uf) : nullptr;
      TTree *tr = (TTree *)fr->Get("cbmsim"), *tg = (TTree *)fg->Get("cbmsim");
      TTree *tp = hasPid ? (TTree *)fp->Get("cbmsim") : nullptr, *tu = hasUkf ? (TTree *)fu->Get("cbmsim") : nullptr;
      TClonesArray *pe = nullptr, *te = nullptr, *pide = nullptr, *ue = nullptr;
      tr->SetBranchAddress("AtPatternEvent", &pe);
      tg->SetBranchAddress("AtTrackingEvent", &te);
      if (tp) tp->SetBranchAddress("AtPIDEvent", &pide);
      if (tu) tu->SetBranchAddress("AtTrackingEvent", &ue);
      Long64_t N = tr->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         tr->GetEntry(i);
         if (pe->GetEntries() == 0) continue;
         auto *pat = (AtPatternEvent *)pe->At(0);
         if (!pat) continue;
         // genfit fitted tracks by trackID
         std::map<int, AtFittedTrack *> gfm;
         if (i < tg->GetEntries()) { tg->GetEntry(i);
            if (te->GetEntries() > 0) { auto *ev = (AtTrackingEvent *)te->At(0);
               if (ev) for (auto &ft : ev->GetFittedTracks()) if (ft) gfm[ft->GetTrackID()] = ft.get(); } }
         std::map<int, AtFittedTrack *> ukm;
         if (tu && i < tu->GetEntries()) { tu->GetEntry(i);
            if (ue->GetEntries() > 0) { auto *ev = (AtTrackingEvent *)ue->At(0);
               if (ev) for (auto &ft : ev->GetFittedTracks()) if (ft) ukm[ft->GetTrackID()] = ft.get(); } }
         std::map<int, const AtTools::AtSpyralResult *> pim;
         if (tp && i < tp->GetEntries()) { tp->GetEntry(i);
            if (pide->GetEntries() > 0) { auto *pv = (AtPIDEvent *)pide->At(0);
               if (pv) for (auto &sr : pv->GetSpyral()) pim[sr.trackID] = &sr; } }

         for (auto &t2 : pat->GetTrackCand()) {
            int tid = t2.GetTrackID();
            float nclus = t2.GetHitClusterArray()->size();
            float nhits = t2.GetHitArray().size();
            float radius = t2.GetGeoRadius();
            float gtheta = t2.GetGeoTheta();
            float gffit = 0, gfke = -1, gfchi2 = -1, gfconv = -1;
            auto ig = gfm.find(tid);
            if (ig != gfm.end()) { gffit = 1; auto &k = ig->second->GetKinematics();
               double nd = ig->second->GetTrackMetadata()->GetNdf(), c2 = ig->second->GetTrackMetadata()->GetChi2();
               gfke = k.kineticEnergy; gfchi2 = nd > 0 ? c2 / nd : -1;
               gfconv = ig->second->GetTrackMetadata()->GetFitConverged(); }
            float ukfit = 0, ukke = -1, ukchi2 = -1;
            auto iu = ukm.find(tid);
            if (iu != ukm.end()) { ukfit = 1; auto &k = iu->second->GetKinematics();
               double nd = iu->second->GetTrackMetadata()->GetNdf(), c2 = iu->second->GetTrackMetadata()->GetChi2();
               ukke = k.kineticEnergy; ukchi2 = nd > 0 ? c2 / nd : -1; }
            float spyvalid = 0, brho = -1, sqdedx = -1, gated = 0;
            auto ip = pim.find(tid);
            if (ip != pim.end()) { spyvalid = ip->second->valid; brho = ip->second->brho; sqdedx = ip->second->sqrtdEdx;
               gated = (ip->second->valid && pid.IsInside(ip->second->sqrtdEdx, ip->second->brho)) ? 1 : 0; }
            float row[18] = {rno, (float)i, (float)tid, nclus, nhits, radius, gtheta,
                             gffit, gfke, gfchi2, gfconv, ukfit, ukke, ukchi2, spyvalid, brho, sqdedx, gated};
            m->Fill(row);
         }
      }
      fr->Close(); fg->Close(); if (fp) fp->Close(); if (fu) fu->Close();
      printf("  %s done (master rows so far: %lld)\n", run.Data(), m->GetEntries());
   }
   fo->cd(); m->Write(); fo->Close();
   printf("wrote %s\n", out.Data());
}
