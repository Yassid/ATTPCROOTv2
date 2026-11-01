/// @file pid_brho_truth_dp.C
/// @brief Brho of the two PID estimators (AtSpyralPID, AtHelix3D) against MC truth, on a simulation whose
///        ejectile is one known species (tritons: 16C(d,t), pdg 1000010030; protons: (d,p), 2212). Only events
///        with exactly ONE track of >= 30 hits are used, so the track is the ejectile without any matching
///        on the estimators' own output. Truth Brho = |p| / (0.299792458 Z) at the generated vertex.
///        Tree "b": tbr (truth), tth (truth theta lab), br/sq (Spyral), b3/s3/th3/sp3 (3D), nh, va, v3.
///        caps = "variant;variant;..." adds 3D variants as bc<i>/sc<i>/vc<i>; a variant is "key=value,..." with keys
///        cap (SetEndWindowCap fraction), len (cap length mm), axis (SetAxisConstraint dcaMax, 0 = off), guard
///        (SetRadiusGuard, 0 = off), hits (SetHitRefit: 0 off, 1 circle, 2 circle + z line).
///   root -l -b -q 'pid_brho_truth_dp.C("/mnt/f/a1975_C16_dt_sim/gs_s3001_reco.root","/mnt/f/a1975_C16_dt_sim/gs_s3001_sim.root",1000010030,"out.root")'
void pid_brho_truth_dp(TString reco, TString sim, int pdg, TString out, Long64_t maxEv = 1000000000, TString caps = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   TFile fs(sim);
   auto *ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchStatus("*", 0);
   ts->SetBranchStatus("MCTrack*", 1);
   ts->SetBranchAddress("MCTrack", &mc);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   std::vector<AtTools::AtHelix3D> hc;
   for (auto *o : *caps.Tokenize(";")) {
      hc.emplace_back();
      double cap = 1, len = 1e9;
      for (auto *kv : *((TObjString *)o)->GetString().Tokenize(",")) {
         TString k = ((TObjString *)kv)->GetString(), key = k(0, k.Index("=")), val = k(k.Index("=") + 1, k.Length());
         if (key == "cap") cap = val.Atof();
         else if (key == "len") len = val.Atof();
         else if (key == "axis") hc.back().SetAxisConstraint(val.Atof());
         else if (key == "guard") hc.back().SetRadiusGuard(val.Atof());
         else if (key == "hits") hc.back().SetHitRefit(val.Atoi() > 0, val.Atoi() > 1);
      }
      hc.back().SetEndWindowCap(cap, len);
   }
   std::vector<double> bc(hc.size()), sc(hc.size());
   std::vector<int> vc(hc.size());
   TFile fo(out, "RECREATE");
   TTree o("b", "estimator Brho vs MC truth");
   int nh, va, v3;
   double tbr, tth, br, sq, b3, s3, th3, sp3;
   o.Branch("tbr", &tbr); o.Branch("tth", &tth); o.Branch("br", &br); o.Branch("sq", &sq); o.Branch("b3", &b3);
   o.Branch("s3", &s3); o.Branch("th3", &th3); o.Branch("sp3", &sp3); o.Branch("nh", &nh); o.Branch("va", &va);
   o.Branch("v3", &v3);
   for (size_t k = 0; k < hc.size(); ++k) {
      o.Branch(Form("bc%zu", k), &bc[k]); o.Branch(Form("sc%zu", k), &sc[k]); o.Branch(Form("vc%zu", k), &vc[k]);
   }
   long nEv = 0, nOne = 0;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      ++nEv;
      const AtTrack *one = nullptr;
      int nLong = 0;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand())
         if (tr.GetHitArray().size() >= 30) { ++nLong; one = &tr; }
      if (nLong != 1) continue;
      ts->GetEntry(ev);
      tbr = -1;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (!m || m->GetPdgCode() != pdg) continue;
         double px = m->GetPx(), py = m->GetPy(), pz = m->GetPz(), pm = std::sqrt(px * px + py * py + pz * pz);
         tbr = pm / 0.299792458; // GeV/c -> T m, Z = 1
         tth = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (tbr <= 0) continue;
      ++nOne;
      nh = one->GetHitArray().size();
      auto sr = spy.Estimate(const_cast<AtTrack &>(*one));
      auto hr = h3.Estimate(*one, sr.valid ? sr.direction : -1);
      va = sr.valid; v3 = hr.valid;
      br = sr.valid ? sr.brho : -1; sq = sr.valid ? sr.sqrtdEdx : -1;
      b3 = hr.valid ? hr.brho : -1; s3 = hr.valid ? hr.sqrtdEdx : -1;
      for (size_t k = 0; k < hc.size(); ++k) {
         auto rc = hc[k].Estimate(*one, sr.valid ? sr.direction : -1);
         vc[k] = rc.valid; bc[k] = rc.valid ? rc.brho : -1; sc[k] = rc.valid ? rc.sqrtdEdx : -1;
      }
      th3 = hr.valid ? hr.theta * TMath::RadToDeg() : -1; sp3 = hr.valid ? hr.spacing : 1e9;
      o.Fill();
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %ld events with tracks, %ld with exactly one >=30-hit track and a truth ejectile\n", out.Data(), nEv, nOne);
}
