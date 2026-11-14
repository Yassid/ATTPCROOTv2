/// @file make_h2_pid_caches.C
/// @brief Two caches for the buttoned gate drawer (pid/gate_draw_dp.C) from the H2 float-tau re-reco (2026-11-11):
///   <out>_peak.root      sqrtdedx = peak amplitude (sqa)
///   <out>_integral.root  sqrtdedx = trace integral (sqw)
/// Same tree as pid_plane_dt.C ("pts": sqrtdedx brho polar vertexz vertexr ncl run ic), so the drawer reads them
/// unchanged. Source: /mnt/f/a1975/reco_h2_tau/pid_width_NNNN.root (tree w: one row per pattern track, va = Spyral
/// valid, br, pol, sqa, sqw, nh, ev = reco entry). IC = max of the FRIB IC trace over TB 1000-1350 at the same entry,
/// as dump_pd_tracks.C (/mnt/f/a1975/reco/run_NNNN_FRIB.root). vertexz/vertexr are not in pid_width: -1.
///   root -l -b -q 'make_h2_pid_caches.C("/mnt/f/a1975/reco_h2_tau/pid_plane_h2")'
void make_h2_pid_caches(TString out = "/mnt/f/a1975/reco_h2_tau/pid_plane_h2",
                        TString pat = "/mnt/f/a1975/reco_h2_tau/pid_width_", TString fribDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kError;
   // every H2 run whose pid_width file is complete (.done), so the cache follows the re-reco as it fills
   std::vector<int> runs;
   for (int r = 106; r <= 189; ++r)
      if (!gSystem->AccessPathName(Form("%s%04d.root.done", pat.Data(), r)))
         runs.push_back(r);
   printf("%zu H2 runs with a finished pid_width\n", runs.size());
   TFile fa(out + "_peak.root", "RECREATE"), fb(out + "_integral.root", "RECREATE");
   float x, y, pol, vz = -1, vr = -1, ncl, ic;
   int run;
   TTree *t[2];
   TFile *fo[2] = {&fa, &fb};
   for (int k = 0; k < 2; ++k) {
      fo[k]->cd();
      t[k] = new TTree("pts", k ? "H2 Spyral PID points, TRACE INTEGRAL (float tau)" : "H2 Spyral PID points, PEAK (float tau)");
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
      TFile *fc = TFile::Open(fribDir + Form("run_%04d_FRIB.root", r));
      TTree *tc = fc ? (TTree *)fc->Get("cbmsim") : nullptr;
      TClonesArray *re = nullptr;
      if (tc)
         tc->SetBranchAddress("AtRawEvent", &re);
      Long64_t lastEv = -1;
      float icEv = -1;
      long n = 0, nic = 0;
      for (Long64_t i = 0; i < w->GetEntries(); ++i) {
         w->GetEntry(i);
         if (va <= 0 || br <= 0.05 || sqa <= 0 || sqw <= 0)
            continue;
         if (ev != lastEv) {
            icEv = -1;
            if (tc && ev < tc->GetEntries()) {
               tc->GetEntry(ev);
               if (re->GetEntries() > 0) {
                  auto *raw = (AtRawEvent *)re->At(0);
                  if (raw && !raw->GetGenTraces().empty()) {
                     auto &adc = raw->GetGenTraces()[0]->GetADC();
                     double mx = -1e9;
                     for (int b = 1000; b < 1350 && b < (int)adc.size(); ++b)
                        mx = std::max(mx, (double)adc[b]);
                     icEv = mx;
                  }
               }
            }
            lastEv = ev;
         }
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
      printf("run %d: %ld valid tracks, IC found %ld; reco entries (pid_width max ev) vs FRIB %lld\n", r, n, nic,
             tc ? tc->GetEntries() : -1LL);
      if (fc)
         fc->Close();
   }
   for (int k = 0; k < 2; ++k) {
      fo[k]->cd();
      t[k]->Write();
      fo[k]->Close();
   }
   printf("-> %s_peak.root, %s_integral.root\n", out.Data(), out.Data());
}
