/// @file gf_hyp_read_dp.C
/// @brief Genfit TWO-HYPOTHESIS PID test (2026-11-01): reads the genfit fits of the same simulated reco file under the
///        deuteron and the triton hypothesis (fitGenfitter_a1975_deuterium.C, CATIMA eloss full range + MSC +
///        straggling, no PID gate) and writes one row per track IN THE ROW ORDER OF pid_sbpid_dp.C (sim mode: events
///        with exactly one >= 30-hit candidate and a truth particle of pdg), so it is a friend of pv_*.root.
///        Per hypothesis h = d, t: ok_h, Brho at the vertex from the back-extrapolated KE (bx_h), chi2_h, ndf_h, pval_h.
///        Duplicate-entry trap (reference_genfitter_tree_duplicate_entries): an entry whose AtPIDEvent holds no Spyral
///        result is a copy of the previous one and is ignored.
///   root -l -b -q 'gf_hyp_read_dp.C("<reco>","<sim>","<gf_d>","<gf_t>",1000010020,"out.root")'
void gf_hyp_read_dp(TString reco, TString sim, TString gfd, TString gft, int pdg, TString out)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   TFile fs(sim);
   auto *ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc);
   TFile *fg[2] = {TFile::Open(gfd), TFile::Open(gft)};
   TTree *tg[2];
   TClonesArray *te[2] = {nullptr, nullptr}, *pe[2] = {nullptr, nullptr};
   const double mass[2] = {1875.613, 2808.921};
   for (int h = 0; h < 2; ++h) {
      tg[h] = (TTree *)fg[h]->Get("cbmsim");
      tg[h]->SetBranchAddress("AtTrackingEvent", &te[h]); tg[h]->SetBranchAddress("AtPIDEvent", &pe[h]);
      printf("%s: %lld entries (reco %lld)\n", h ? "t" : "d", tg[h]->GetEntries(), t->GetEntries());
   }
   TFile fo(out, "RECREATE");
   TTree o("g", "genfit d/t hypotheses");
   double tbr, bx[2], b0[2], c2[2], nd[2], pv[2];
   int ok[2];
   o.Branch("tbr", &tbr);
   for (int h = 0; h < 2; ++h) {
      const char *s = h ? "t" : "d";
      o.Branch(Form("ok_%s", s), &ok[h]); o.Branch(Form("bx_%s", s), &bx[h]); o.Branch(Form("b0_%s", s), &b0[h]);
      o.Branch(Form("chi2_%s", s), &c2[h]); o.Branch(Form("ndf_%s", s), &nd[h]); o.Branch(Form("pval_%s", s), &pv[h]);
   }
   long nRow = 0, nFit[2] = {0, 0};
   for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto &cands = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      int nl = 0;
      for (auto &tr : cands) nl += tr.GetHitArray().size() >= 30;
      if (nl != 1) continue;
      ts->GetEntry(ev);
      tbr = -1;
      for (int j = 0; j < mc->GetEntriesFast(); j++) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (!m || m->GetPdgCode() != pdg) continue;
         tbr = std::sqrt(m->GetPx() * m->GetPx() + m->GetPy() * m->GetPy() + m->GetPz() * m->GetPz()) / 0.299792458;
         break;
      }
      if (tbr <= 0) continue;
      std::map<int, AtFittedTrack *> fm[2];
      for (int h = 0; h < 2; ++h) {
         if (ev >= tg[h]->GetEntries()) continue;
         tg[h]->GetEntry(ev);
         if (!pe[h]->GetEntries() || !te[h]->GetEntries()) continue;
         auto *pid = (AtPIDEvent *)pe[h]->At(0);
         auto *tev = (AtTrackingEvent *)te[h]->At(0);
         if (!pid || !tev || pid->GetSpyral().empty()) continue; // duplicate copy
         for (auto &ft : tev->GetFittedTracks()) if (ft) fm[h][ft->GetTrackID()] = ft.get();
      }
      for (auto &tr : cands) {
         if (tr.GetHitArray().size() < 30) continue;
         for (int h = 0; h < 2; ++h) {
            ok[h] = 0; bx[h] = b0[h] = c2[h] = nd[h] = pv[h] = -1;
            auto it = fm[h].find(tr.GetTrackID());
            if (it == fm[h].end()) continue;
            auto *ft = it->second;
            auto kx = ft->GetKinematicsXtr(), k0 = ft->GetKinematics();
            auto P = [&](double ke) { return ke > 0 ? std::sqrt(ke * ke + 2 * ke * mass[h]) / 299.792458 : -1; };
            bx[h] = P(kx.kineticEnergy); b0[h] = P(k0.kineticEnergy);
            c2[h] = ft->GetTrackMetadata()->GetChi2(); nd[h] = ft->GetTrackMetadata()->GetNdf();
            pv[h] = ft->GetTrackMetadata()->GetPValue();
            ok[h] = bx[h] > 0 && nd[h] > 0;
            nFit[h] += ok[h];
         }
         o.Fill(); ++nRow;
      }
   }
   fo.cd(); o.Write();
   printf("rows %ld, fitted d %ld t %ld\n", nRow, nFit[0], nFit[1]);
}
