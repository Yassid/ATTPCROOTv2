/// @file dir_ab_dp.C
/// @brief A/B of genfit's forward/backward call on a1975 16C(d,p): AtPRA GeoTheta (ctrl) vs the
///        AtSpyralPID estimate (spy, AtGenfitter::SetDirectionFromSpyral), matched track by track.
///
/// Both arms are find805's arguments verbatim (gf_dp_dirAB/run_ab.sh), so the ONLY difference is
/// the direction source. Sanity check built in: where the two calls agree, the fits must match.
/// Per gated Spyral-valid track (ex_dp_dv1104.C's matching rule, which also skips the duplicate
/// genfitter entries): back-extrapolated KE/theta (Xtr slot), Ex at constant Ebeam, chi2/ndf,
/// vertex z and |xy|, plus both direction calls (GeoTheta from the reco file).
///
///   root -l -b -q 'dir_ab_dp.C("run_0016,run_0020")'
#include "ex_dp_dv1104.C"

void dir_ab_dp(TString runsCSV = "run_0016,run_0020", TString abDir = "/mnt/f/a1975/gf_dp_dirAB/",
               TString recoDir = "/mnt/f/a1975/reco_d2_dv1104/", TString gateFile = "pid/proton_d2_dv1104.json",
               double Ebeam = 184.25, TString out = "/mnt/f/a1975/gf_dp_dirAB/dir_ab_dp.root")
{
   gSystem->Load("libAtReconstruction.so");
   gErrorIgnoreLevel = kFatal;
   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;
   const double m_p = 1.00727646688 * u, m_C17 = 17.0225864 * u;
   auto gate = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile fo(out, "RECREATE");
   TTree o("ab", "direction A/B, (d,p)");
   int run, tid, geoBack, spyBack, ncl, F[2];
   Long64_t ev;
   double sbr, ssq, KE[2], TH[2], EX[2], TCM[2], C2[2], VZ[2], VR[2];
   o.Branch("run", &run); o.Branch("ev", &ev); o.Branch("tid", &tid); o.Branch("ncl", &ncl);
   o.Branch("geoBack", &geoBack); o.Branch("spyBack", &spyBack); o.Branch("sbr", &sbr); o.Branch("ssq", &ssq);
   const char *arm[2] = {"c", "s"};
   for (int a = 0; a < 2; ++a) {
      o.Branch(Form("%sfit", arm[a]), &F[a]); o.Branch(Form("%ske", arm[a]), &KE[a]);
      o.Branch(Form("%sth", arm[a]), &TH[a]); o.Branch(Form("%sex", arm[a]), &EX[a]);
      o.Branch(Form("%stcm", arm[a]), &TCM[a]); o.Branch(Form("%schi", arm[a]), &C2[a]);
      o.Branch(Form("%svz", arm[a]), &VZ[a]); o.Branch(Form("%svr", arm[a]), &VR[a]);
   }

   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString rs = ((TObjString *)runs->At(ri))->GetString();
      run = TString(rs(rs.Length() - 4, 4)).Atoi();
      TFile *fr = TFile::Open(recoDir + rs + "_multifit_reco.root");
      TFile *fc = TFile::Open(abDir + "ctrl/" + rs + "_multifit_genfitter_p.root");
      TFile *fs = TFile::Open(abDir + "spy/" + rs + "_multifit_genfitter_p.root");
      if (!fr || !fc || !fs) { printf("skip %s\n", rs.Data()); continue; }
      auto *tr = (TTree *)fr->Get("cbmsim");
      TTree *tg[2] = {(TTree *)fc->Get("cbmsim"), (TTree *)fs->Get("cbmsim")};
      TClonesArray *pa = nullptr, *te[2] = {nullptr, nullptr}, *tp = nullptr;
      tr->SetBranchStatus("*", 0); tr->SetBranchStatus("AtPatternEvent*", 1);
      tr->SetBranchAddress("AtPatternEvent", &pa);
      for (int a = 0; a < 2; ++a) tg[a]->SetBranchAddress("AtTrackingEvent", &te[a]);
      tg[0]->SetBranchAddress("AtPIDEvent", &tp);
      Long64_t N = std::min({tr->GetEntries(), tg[0]->GetEntries(), tg[1]->GetEntries()});
      printf("%s: entries reco %lld ctrl %lld spy %lld\n", rs.Data(), tr->GetEntries(), tg[0]->GetEntries(),
             tg[1]->GetEntries());
      for (ev = 0; ev < N; ++ev) {
         tg[0]->GetEntry(ev);
         auto *pe = (tp && tp->GetEntries()) ? (AtPIDEvent *)tp->At(0) : nullptr;
         if (!pe || pe->GetSpyral().empty()) continue;
         tg[1]->GetEntry(ev);
         tr->GetEntry(ev);
         auto *pat = (pa && pa->GetEntries()) ? (AtPatternEvent *)pa->At(0) : nullptr;
         std::map<int, double> geo;
         if (pat)
            for (auto &t : pat->GetTrackCand()) geo[t.GetTrackID()] = t.GetGeoTheta();
         std::map<int, AtFittedTrack *> fm[2];
         for (int a = 0; a < 2; ++a) {
            auto *e = (te[a] && te[a]->GetEntries()) ? (AtTrackingEvent *)te[a]->At(0) : nullptr;
            if (e)
               for (auto &ft : e->GetFittedTracks())
                  if (ft) fm[a][ft->GetTrackID()] = ft.get();
         }
         for (auto &sr : pe->GetSpyral()) {
            if (!sr.valid || !gate.IsInside(sr.sqrtdEdx, sr.brho)) continue;
            tid = sr.trackID;
            auto g = geo.find(tid);
            if (g == geo.end()) continue;
            geoBack = g->second > M_PI / 2;
            spyBack = (M_PI - sr.polar) > M_PI / 2;
            sbr = sr.brho; ssq = sr.sqrtdEdx; ncl = sr.nClusters;
            for (int a = 0; a < 2; ++a) {
               F[a] = 0;
               KE[a] = TH[a] = EX[a] = TCM[a] = C2[a] = VZ[a] = VR[a] = -999;
               auto it = fm[a].find(tid);
               if (it == fm[a].end()) continue;
               auto *ft = it->second;
               auto &k = ft->GetKinematicsXtr();
               double ndf = ft->GetTrackMetadata()->GetNdf();
               F[a] = 1;
               KE[a] = k.kineticEnergy; TH[a] = k.theta * TMath::RadToDeg();
               C2[a] = ndf > 0 ? ft->GetTrackMetadata()->GetChi2() / ndf : 1e9;
               auto v = ft->GetVertex();
               VZ[a] = v.Z(); VR[a] = std::hypot(v.X(), v.Y());
               if (KE[a] > 0) {
                  auto [ex, tcm] = kine_2b(m_C16, m_d, m_p, m_C17, Ebeam, k.theta, KE[a]);
                  if (!std::isnan(ex)) { EX[a] = ex; TCM[a] = tcm; }
               }
            }
            o.Fill();
         }
      }
      fr->Close(); fc->Close(); fs->Close();
   }
   fo.cd(); o.Write(); fo.Close();
   printf("wrote %s\n", out.Data());
}
