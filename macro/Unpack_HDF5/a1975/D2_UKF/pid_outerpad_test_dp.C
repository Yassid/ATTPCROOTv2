/// @file pid_outerpad_test_dp.C
/// @brief Does the 3D dE/dx improve if it also uses the OUTER (large) pads? Same tracks as
///        pid_padgain_test_dp.C (same branches) plus s2/b2/ok2/dir1/dir2: AtHelix3D, pad-gain corrected, with
///        the small-pad radius raised to rOuter (all pads), and the vertex-end call of both 3D versions.
///        Fast forward tracks leave rho <= 152 mm early, so the inner-only 3D dE/dx is short for them.
///        p1/p2: dE/dx path of the two versions (outer charge per mm = (s2^2 p2 - s1^2 p1)/(p2 - p1)).
///        s3/b3o/ok3: all pads with the outer-pad (rho > 152 mm) charges divided by kOuter (inner/outer scale).
///   root -l -b -q 'pid_outerpad_test_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2_v2.txt","out.root")'
// col: which gain column of the map to use -- 1 = full gain, 4 = local (ring-normalised) gain
void pid_outerpad_test_dp(TString in, TString padMap, TString out, Long64_t maxEv = 1000000000, int col = 1,
                          double rOuter = 1000, double kOuter = 1)
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
   AtTools::AtHelix3D h3, h3o;
   h3o.SetSmallPadRadius(rOuter);
   TFile fo(out, "RECREATE");
   TTree o("w", "PID raw vs pad-gain corrected");
   int nh, va, vp, v3, ok1, ok2, ok3, dir1, dir2;
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s0, b0, s1, b1, s2, b2, s3, b3o, p1, p2;
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s0", &s0); o.Branch("b0", &b0); o.Branch("s1", &s1);
   o.Branch("b1", &b1); o.Branch("ok1", &ok1); o.Branch("s2", &s2); o.Branch("b2", &b2); o.Branch("ok2", &ok2);
   o.Branch("dir1", &dir1); o.Branch("dir2", &dir2); o.Branch("s3", &s3); o.Branch("b3o", &b3o); o.Branch("ok3", &ok3);
   o.Branch("p1", &p1); o.Branch("p2", &p2);
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
         auto r2 = h3o.Estimate(tc, sr.valid ? sr.direction : -1);
         ok2 = r2.valid; s2 = r2.valid ? r2.sqrtdEdx : -1; b2 = r2.valid ? r2.brho : -1;
         p1 = r1.valid ? r1.pidPath : -1; p2 = r2.valid ? r2.pidPath : -1;
         AtTrack to(tc);
         for (auto &h : to.GetHitArray())
            if (std::hypot(h->GetPosition().X(), h->GetPosition().Y()) > 152) h->SetCharge(h->GetCharge() / kOuter);
         auto r3 = h3o.Estimate(to, sr.valid ? sr.direction : -1);
         ok3 = r3.valid; s3 = r3.valid ? r3.sqrtdEdx : -1; b3o = r3.valid ? r3.brho : -1;
         dir1 = r1.valid ? int(r1.backward) : -9; dir2 = r2.valid ? int(r2.backward) : -9;
         b3 = b0; sq3h = s0; th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
