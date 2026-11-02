/// @file validate_pid_framework_dp.C
/// @brief Track-by-track check that the framework PID options (2026-11-02) reproduce the a1975 analysis macros:
///        AtPIDGainTable::LoadPadGains + AtHelix3D (SetHitRefit(true,true), SetArcAngle) + AtTools::CombinePID, run on
///        one run's reco file in the row order of pid_sbpid_dp.C / range_dp.C (every >= 30-hit track), against
///          bA (pid_sbpid_dp.C)         <-> brhoArc        hs <-> sqrtdEdx (pad-corrected, no run gain)
///          r.midpeak (range_dp.C, RAW charges) <-> midPeak   r.cont && 10 < zfar < 955 <-> contained
///          Spyral veto (plot_pid_spyral_final_dp.C, sqp_g = Spyral x run gain) <-> CombinePID vetoed
///   root -l -b -q 'validate_pid_framework_dp.C("0026",3000)'
void validate_pid_framework_dp(TString run = "0026", Long64_t maxEv = 3000, TString padMap = "pid/padgain_d2_v1.txt")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kWarning;
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   Long64_t off = 0;
   for (int i = 0; i < 4 && run != runs[i]; ++i) {
      TFile fr(Form("/mnt/f/a1975/pid3d_sbpid/pv_%s.root", runs[i]));
      off += ((TTree *)fr.Get("w"))->GetEntries();
   }
   TFile fp("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root");
   auto *w = (TTree *)fp.Get("w");
   w->AddFriend("r=r", "/mnt/f/a1975/pid_range/R_all.root");
   double bA, hs, br, sqp_g, g_sqp, zfar;
   int va, midpeak, cont, rok;
   w->SetBranchAddress("bA", &bA); w->SetBranchAddress("hs", &hs); w->SetBranchAddress("br", &br); w->SetBranchAddress("sqp_g", &sqp_g);
   w->SetBranchAddress("g_sqp", &g_sqp); w->SetBranchAddress("va", &va);
   w->SetBranchAddress("r.midpeak", &midpeak); w->SetBranchAddress("r.cont", &cont); w->SetBranchAddress("r.zfar", &zfar); w->SetBranchAddress("r.ok", &rok);
   auto gains = AtTools::AtPIDGainTable::LoadPadGains(padMap.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   h3.SetBField(2.85);
   h3.SetHitRefit(true, true);
   h3.SetArcAngle(true);
   AtTools::AtPIDCombineOpt opt;
   opt.highSq = 30; opt.addBackScale = 0.965;
   TFile f(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", run.Data()));
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   long n = 0, nB = 0, nBagree = 0, nS = 0, nSagree = 0, nM = 0, nMagree = 0, nC = 0, nCagree = 0, nV = 0, nVagree = 0, nVmac = 0, nVfw = 0, nMpad = 0, nCpad = 0;
   Long64_t row = off;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      for (auto &tr0 : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         if (tr0.GetHitArray().size() < 30) continue;
         w->GetEntry(row++);
         AtTrack tc(tr0);
         for (auto &h : tc.GetHitArray()) { auto it = gains.find(h->GetPadNum()); if (it != gains.end()) h->SetCharge(h->GetCharge() / it->second); }
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr0)); // the macros took Spyral's direction from the raw track
         auto sp = spy.Estimate(tc);
         auto r = h3.Estimate(tc, sr.valid ? sr.direction : -1);
         ++n;
         if (bA > 0 || r.brhoArc > 0) { ++nB; nBagree += bA > 0 && r.brhoArc > 0 && std::abs(r.brhoArc / bA - 1) < 1e-6; }
         if (hs > 0 || r.sqrtdEdx > 0) { ++nS; nSagree += hs > 0 && std::abs(r.sqrtdEdx / hs - 1) < 1e-6; }
         if (rok) {
            // range_dp.C ran on the RAW charges (no pad map): compare like with like
            auto rr = h3.Estimate(tr0, sr.valid ? sr.direction : -1);
            ++nM; nMagree += (bool)midpeak == rr.midPeak; nMpad += (bool)midpeak != r.midPeak;
            ++nC; nCagree += (bool)(cont && zfar > 10 && zfar < 955) == rr.contained; nCpad += (bool)(cont && zfar > 10 && zfar < 955) != r.contained;
         }
         if (va && bA > 0 && br > 0) {
            ++nV;
            bool mac = std::log(br / bA) >= 0.25 && (br < 0.45 || sqp_g > 30);
            AtTools::AtSpyralResult spg = sp;
            spg.brho = sr.brho;            // macro: Brho from the raw track
            spg.sqrtdEdx = sp.sqrtdEdx * g_sqp; // macro: sqp_g = pad-corrected Spyral x run gain
            bool fw = AtTools::CombinePID(spg, r, opt).vetoed;
            nVagree += mac == fw; nVmac += mac; nVfw += fw;
         }
      }
   }
   printf("run %s, %lld events, %ld tracks\n", run.Data(), std::min(maxEv, t->GetEntries()), n);
   printf("  brhoArc == bA        %ld / %ld\n", nBagree, nB);
   printf("  sqrtdEdx == hs       %ld / %ld\n", nSagree, nS);
   printf("  midPeak == midpeak   %ld / %ld (raw charges, as range_dp.C); with the pad map %ld differ\n", nMagree, nM, nMpad);
   printf("  contained == cont    %ld / %ld (raw charges); with the pad map %ld differ\n", nCagree, nC, nCpad);
   printf("  veto agrees          %ld / %ld (macro vetoes %ld, framework %ld)\n", nVagree, nV, nVmac, nVfw);
}
