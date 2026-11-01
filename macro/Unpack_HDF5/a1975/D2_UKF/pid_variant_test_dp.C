/// @file pid_variant_test_dp.C
/// @brief AtHelix3D setting variants on data, same tracks and same branches as pid_padgain_test_dp.C (sqa/br
///        Spyral raw, sqp Spyral pad-corrected, s0/b0 3D raw, s1/b1 3D pad-corrected, reference v3 b3 sq3h th3
///        sp3) plus sv<k>/bv<k>/vv<k>: pad-corrected 3D with variant k. variants = "key=value,...;..." with keys
///        hits (SetHitRefit: 0 off, 1 circle, 2 circle + z line), cap (SetEndWindowCap fraction), axis (SetAxisConstraint), guard.
///   root -l -b -q 'pid_padgain_test_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","pid/padgain_d2.txt","out.root")'
// col: which gain column of the map to use -- 1 = full gain, 4 = local (ring-normalised) gain
void pid_variant_test_dp(TString in, TString padMap, TString out, TString variants, Long64_t maxEv = 1000000000,
                         int col = 1)
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
   std::vector<AtTools::AtHelix3D> hv;
   for (auto *o : *variants.Tokenize(";")) {
      hv.emplace_back();
      for (auto *kv : *((TObjString *)o)->GetString().Tokenize(",")) {
         TString k = ((TObjString *)kv)->GetString(), key = k(0, k.Index("=")), val = k(k.Index("=") + 1, k.Length());
         if (key == "hits") hv.back().SetHitRefit(val.Atoi() > 0, val.Atoi() > 1);
         else if (key == "cap") hv.back().SetEndWindowCap(val.Atof());
         else if (key == "axis") hv.back().SetAxisConstraint(val.Atof());
         else if (key == "guard") hv.back().SetRadiusGuard(val.Atof());
      }
   }
   std::vector<double> sv(hv.size()), bv(hv.size());
   std::vector<int> vv(hv.size());
   TFile fo(out, "RECREATE");
   TTree o("w", "PID raw vs pad-gain corrected");
   int nh, va, vp, v3, ok1;
   double br, sqa, sqp, pol, b3, sq3h, th3, sp3, s0, b0, s1, b1;
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vp", &vp); o.Branch("br", &br); o.Branch("sqa", &sqa);
   o.Branch("sqp", &sqp); o.Branch("pol", &pol); o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h);
   o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("s0", &s0); o.Branch("b0", &b0); o.Branch("s1", &s1);
   o.Branch("b1", &b1); o.Branch("ok1", &ok1);
   for (size_t k = 0; k < hv.size(); ++k) {
      o.Branch(Form("sv%zu", k), &sv[k]); o.Branch(Form("bv%zu", k), &bv[k]); o.Branch(Form("vv%zu", k), &vv[k]);
   }
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
         for (size_t k = 0; k < hv.size(); ++k) {
            auto rv = hv[k].Estimate(tc, sr.valid ? sr.direction : -1);
            vv[k] = rv.valid; sv[k] = rv.valid ? rv.sqrtdEdx : -1; bv[k] = rv.valid ? rv.brho : -1;
         }
         b3 = b0; sq3h = s0; th3 = r0.valid ? r0.theta * TMath::RadToDeg() : -1; sp3 = r0.valid ? r0.spacing : 1e9;
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
