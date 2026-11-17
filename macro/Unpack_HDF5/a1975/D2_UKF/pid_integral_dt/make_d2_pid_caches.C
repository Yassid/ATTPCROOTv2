/// @file make_d2_pid_caches.C
/// @brief 2026-11-17 D2 twin of make_h2_pid_caches.C, for the (d,t) triton gate on the HDBSCAN + float-tau D2 re-reco
///        (/mnt/f/a1975/reco_d2_hdbtau/pid_width_NNNN.root). IC from /mnt/f/a1975/ic_d2/run_NNNN_IC.root (tree ic, icmax), which is
///        ALIGNED WITH THE RECO ENTRIES (project_a1975_dp_trigger traps), so ic = icmax[ev]. H2 header follows:
/// Two caches for the buttoned gate drawer (pid/gate_draw_dp.C) from the H2 float-tau re-reco (2026-11-11):
///   <out>_peak.root      sqrtdedx = peak amplitude (sqa)
///   <out>_integral.root  sqrtdedx = trace integral (sqw)
/// Same tree as pid_plane_dt.C ("pts": sqrtdedx brho polar vertexz vertexr ncl run ic), so the drawer reads them
/// unchanged. Source: /mnt/f/a1975/reco_h2_tau/pid_width_NNNN.root (tree w: one row per pattern track, va = Spyral
/// valid, br, pol, sqa, sqw, nh, ev = reco entry). IC = max of the FRIB IC trace over TB 1000-1350 at the same entry,
/// as dump_pd_tracks.C (/mnt/f/a1975/reco/run_NNNN_FRIB.root). vertexz/vertexr are not in pid_width: -1.
///   root -l -b -q 'make_h2_pid_caches.C("/mnt/f/a1975/reco_h2_tau/pid_plane_h2")'
void make_d2_pid_caches(TString out = "/mnt/f/a1975/reco_d2_hdbtau/pid_plane_d2",
                        TString pat = "/mnt/f/a1975/reco_d2_hdbtau/pid_width_", TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kError;
   // every H2 run whose pid_width file is complete (.done), so the cache follows the re-reco as it fills
   std::vector<int> runs;
   for (int r = 16; r <= 103; ++r)
      if (!gSystem->AccessPathName(Form("%s%04d.root.done", pat.Data(), r)))
         runs.push_back(r);
   printf("%zu D2 runs with a finished pid_width\n", runs.size());
   TFile fa(out + "_peak.root", "RECREATE"), fb(out + "_integral.root", "RECREATE");
   float x, y, pol, vz = -1, vr = -1, ncl, ic;
   int run;
   TTree *t[2];
   TFile *fo[2] = {&fa, &fb};
   for (int k = 0; k < 2; ++k) {
      fo[k]->cd();
      t[k] = new TTree("pts", k ? "D2 Spyral PID points, TRACE INTEGRAL (float tau)" : "D2 Spyral PID points, PEAK (float tau)");
      t[k]->Branch("sqrtdedx", &x);
      t[k]->Branch("brho", &y);
      t[k]->Branch("polar", &pol);
      t[k]->Branch("vertexz", &vz);
      t[k]->Branch("vertexr", &vr);
      t[k]->Branch("ncl", &ncl);
      t[k]->Branch("run", &run);
      t[k]->Branch("ic", &ic);
   }
   for (int r : runs) {
      TFile fw(Form("%s%04d.root", pat.Data(), r));
      auto *w = (TTree *)fw.Get("w");
      if (!w) {
         printf("run %d: no pid_width\n", r);
         continue;
      }
      Long64_t ev;
      int va, nh;
      double br, pl, sqa, sqw;
      w->SetBranchAddress("ev", &ev);
      w->SetBranchAddress("va", &va);
      w->SetBranchAddress("nh", &nh);
      w->SetBranchAddress("br", &br);
      w->SetBranchAddress("pol", &pl);
      w->SetBranchAddress("sqa", &sqa);
      w->SetBranchAddress("sqw", &sqw);
      std::vector<float> icv;
      {
         TFile *fic = TFile::Open(icDir + Form("run_%04d_IC.root", r));
         if (fic && !fic->IsZombie()) {
            auto *ti = (TTree *)fic->Get("ic"); float icm = 0; ti->SetBranchAddress("icmax", &icm);
            for (Long64_t j = 0; j < ti->GetEntries(); ++j) { ti->GetEntry(j); icv.push_back(icm); }
            fic->Close();
         } else printf("  *** run %d: NO IC file -> ic = -1 ***\n", r);
      }
      float icEv = -1;
      long n = 0, nic = 0;
      for (Long64_t i = 0; i < w->GetEntries(); ++i) {
         w->GetEntry(i);
         if (va <= 0 || br <= 0.05 || sqa <= 0 || sqw <= 0)
            continue;
         icEv = (ev >= 0 && ev < (Long64_t)icv.size()) ? icv[ev] : -1.f;
         run = r;
         y = br;
         pol = pl;
         ncl = nh;
         ic = icEv;
         nic += ic > 0;
         x = sqa;
         t[0]->Fill();
         x = sqw;
         t[1]->Fill();
         ++n;
      }
      printf("run %d: %ld valid tracks, IC found %ld (IC entries %zu)\n", r, n, nic, icv.size());
   }
   for (int k = 0; k < 2; ++k) {
      fo[k]->cd();
      t[k]->Write();
      fo[k]->Close();
   }
   printf("-> %s_peak.root, %s_integral.root\n", out.Data(), out.Data());
}
