/// @file pid_padgain_test_dp.C
/// @brief Does the pad-gain map (pad_gain_map_dp.C) improve the PID? Same tracks, Spyral and AtHelix3D
///        computed on raw hit amplitudes and on amplitudes divided by the pad gain. Branch names follow
///        pid_variants_dp.C (sqa/br Spyral raw, s0/b0 3D raw, v3 b3 sq3h th3 sp3 reference) plus
///        sqp (Spyral, pad-corrected) and s1/b1 (3D, pad-corrected), so the usual scorers apply.
///        Run it on runs NOT used to make the map.
///   root -l -b -q 'pid_padgain_test_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2.txt","out.root")'
// col: which gain column of the map to use -- 1 = full gain, 4 = local (ring-normalised) gain
void pid_padgain_test_dp(TString in, TString padMap, TString out, Long64_t maxEv = 1000000000, int col = 1)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   std::map<int, double> gain;
   {
      std::ifstream fm(padMap.Data());
      std::string line;
      while (std::getline(fm, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         int p; double g; long n; double x, y, gl;
         if (!(ss >> p >> g >> n)) continue;
         if (col == 4) { if (ss >> x >> y >> gl) gain[p] = gl; } else gain[p] = g;
      }
   }
   printf("pad map: %zu pads\n", gain.size());
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   TFile fo(out, "RECREATE");
   TTree o("w", "PID raw vs pad-gain corrected");
   int nh, va, vp, v3, ok1;
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s0, b0, s1, b1;
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s0", &s0); o.Branch("b0", &b0); o.Branch("s1", &s1);
   o.Branch("b1", &b1); o.Branch("ok1", &ok1);
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         nh = H.size();
         if (nh < 30) continue;
         // pad-corrected copy of the track
         AtTrack tc(tr);
         for (auto &h : tc.GetHitArray()) {
            auto it = gain.find(h->GetPadNum());
            if (it != gain.end() && it->second > 0) h->SetCharge(h->GetCharge() / it->second);
         }
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
         auto sp = spy.Estimate(tc);
         va = sr.valid; vp = sp.valid;
         br = sr.valid ? sr.brho : -1; sqa = sr.valid ? sr.sqrtdEdx : -1; sqp = sp.valid ? sp.sqrtdEdx : -1;
         pol = sr.valid ? sr.polar * TMath::RadToDeg() : -999;
         auto r0 = h3.Estimate(tr, sr.valid ? sr.direction : -1);
         auto r1 = h3.Estimate(tc, sr.valid ? sr.direction : -1);
         v3 = r0.valid; ok1 = r1.valid;
         s0 = r0.valid ? r0.sqrtdEdx : -1; b0 = r0.valid ? r0.brho : -1;
         s1 = r1.valid ? r1.sqrtdEdx : -1; b1 = r1.valid ? r1.brho : -1;
         b3 = b0; sq3h = s0; th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
