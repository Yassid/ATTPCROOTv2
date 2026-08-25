/// @file closure_RT_C14.C
/// @brief Closure of the 2026-08-25 14C fitter change against MC TRUTH, on the RT-gas sample.
///        Both arms are fitted from the SAME reco, so the only difference is the fitter config:
///          _cat  matEffects + native CATIMA dE/dx, ATTPC_H300torr_RT, matFallback OFF
///          _xtr  matEffects OFF + manual CATIMA eloss over the vertex gap  (the ADOPTED data config)
///        Truth matching follows fit_truth_bias_C14.C: primary proton from MCTrack, paired with the
///        fitted track in the same entry closest in theta, |dtheta| < 10 deg.
///   root -b -q 'closure_RT_C14.C'
void closure_RT_C14(TString dir = "diagnostics/RT/", Double_t chi2Cut = 1e9)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   const double u = 931.49401, m_p = 1.007825 * u;
   const char *arms[2] = {"xtr", "cat"};
   const char *lab[2]  = {"matFX OFF + manual gap (ADOPTED)", "matFX ON + CATIMA (NEW)"};
   for (int a = 0; a < 2; ++a) {
      TFile *fs = TFile::Open(dir + "attpcsim_RT.root");
      TFile *ff = TFile::Open(TString::Format("%sRT_genfit_%s.root", dir.Data(), arms[a]));
      if (!fs || fs->IsZombie() || !ff || ff->IsZombie()) { printf("missing files for %s\n", arms[a]); continue; }
      TTree *ts = (TTree *)fs->Get("cbmsim"), *tf = (TTree *)ff->Get("cbmsim");
      if (!ts || !tf) { printf("no tree for %s\n", arms[a]); continue; }
      TClonesArray *mc = nullptr, *te = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);
      Long64_t N = std::min(ts->GetEntries(), tf->GetEntries());
      std::vector<double> rel, relX, dth, c2v; long nTrk = 0, nColl = 0;
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i); tf->GetEntry(i);
         if (!mc || !te || te->GetEntriesFast() == 0) continue;
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *tr = (AtMCTrack *)mc->At(k);
            if (!tr || tr->GetPdgCode() != 2212 || tr->GetMotherId() != -1) continue;
            double px = tr->GetPx()*1000, py = tr->GetPy()*1000, pz = tr->GetPz()*1000;
            double p = std::sqrt(px*px+py*py+pz*pz); if (p <= 0) continue;
            keT = std::sqrt(p*p+m_p*m_p) - m_p; thT = std::acos(pz/p)*TMath::RadToDeg(); break;
         }
         if (keT <= 0) continue;
         auto *evt = (AtTrackingEvent *)te->At(0); if (!evt) continue;
         double best = 1e9, keR = -1, keX = -1, thR = -1, c2 = -1;
         for (auto &ft : evt->GetFittedTracks()) {
            if (!ft) continue;
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double r = (ndf > 0) ? chi2/ndf : 1e9;
            nTrk++; if (r >= 1e8) { nColl++; continue; }
            if (r > chi2Cut) continue;
            double th = ft->GetKinematics().theta * TMath::RadToDeg();
            if (std::fabs(th - thT) < best) { best = std::fabs(th - thT); keR = ft->GetKinematics().kineticEnergy;
                                              keX = ft->GetKinematicsXtr().kineticEnergy; thR = th; c2 = r; }
         }
         if (keR <= 0 || best > 10.0) continue;
         rel.push_back(100.0*(keR-keT)/keT);
         if (keX > 0 && keX < 1000) relX.push_back(100.0*(keX-keT)/keT);
         dth.push_back(thR - thT); c2v.push_back(c2);
      }
      auto med=[](std::vector<double> v){ if(v.empty()) return -999.0; std::sort(v.begin(),v.end()); return v[v.size()/2]; };
      auto iqr=[](std::vector<double> v){ if(v.size()<4) return -1.0; std::sort(v.begin(),v.end());
                                          return v[3*v.size()/4]-v[v.size()/4]; };
      printf("%-34s pairs=%5zu  tracks=%6ld coll=%ld(%.2f%%)  medChi2=%.3f\n"
             "%-34s   KE bias      = %+6.2f %%  (IQR %.2f)\n"
             "%-34s   KE_xtr bias  = %+6.2f %%  (IQR %.2f)\n"
             "%-34s   theta bias   = %+6.3f deg (IQR %.3f)\n",
             lab[a], rel.size(), nTrk, nColl, nTrk?100.*nColl/nTrk:0., med(c2v),
             "", med(rel), iqr(rel), "", med(relX), iqr(relX), "", med(dth), iqr(dth));
      fs->Close(); ff->Close();
   }
}
