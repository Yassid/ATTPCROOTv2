/// @file embed_dd.C
/// @brief EMBEDDING of simulated 16C(d,d) elastic deuterons into real D2 events (2026-11-10).
/// The elastic scale S = data / (L eps sigma_OM) sits at ~0.76 with the r280 sim acceptance, while Spyral validity
/// (dd_valid_frac.py) and the genfit fit efficiency (dd_fit_eff.py) on the locus are ~0.98-1.0 in data. The remaining
/// stage is finding the track in the real event environment (beam, noise, other tracks), which the sim lacks. Same
/// method as embed_vz_dp.C (which gave .68-.91 for the (d,p) protons): the sim DEUTERON PSA hits (AtEventH, before
/// cleaning, MC A=2 Z=1) are overlaid on random real events of a data run (AtEventH, same stage; data reco has space
/// charge off), and the chain is re-run: AtDirDeDxCleaner -> AtTrackFinderHDBSCAN (data arguments) -> AtSpyralPID.
///   mode 1  deuteron hits alone      -> reference: the Spyral track on the elastic locus (if none, event skipped)
///   mode 2  deuteron hits + data event (nOver per sim event)
/// Success in mode 2 = a valid Spyral track with theta_lab 60-89.5, |theta_cm - ref| < 3 deg, |vz - ref| < 40 mm,
/// |KE/KE_ref - 1| < 0.15. theta_cm = 180 - 2 theta_lab (elastic, non-relativistic; only differences matter).
/// Tree e: sev dev mode found ref(1 if mode 1 found a locus track) tc0 vz0 ke0 tcB vzB keB nhS nhD ntr icD
///   root -l -b -q 'embed_dd.C(16,2,"embed_dd_r16.root","/mnt/f/a1975/sim_rebuild/r280/dd_el/el_s6001_reco.root")'
bool isDeut(const AtHit &h)
{
   for (auto &m : h.GetMCSimPointArray())
      if (m.A == 2 && m.Z == 1)
         return true;
   return false;
}

struct Spy {
   bool ok;
   double th, vz, ke, tc;
};
Spy spyKin(const AtTools::AtSpyralResult &r)
{
   Spy s{false, -999, -999, -999, -999};
   if (!r.valid)
      return s;
   const double md = 1875.613;
   s.th = 180 - r.polar * TMath::RadToDeg();
   s.vz = 1000 - r.vertex.Z();
   double p = 299.792458 * r.brho;
   s.ke = std::sqrt(p * p + md * md) - md;
   s.tc = 180 - 2 * s.th;
   s.ok = s.th > 60 && s.th < 89.5 && s.ke > 0.3 && s.ke < 10;
   return s;
}

void embed_dd(int dataRun = 16, int nOver = 2, TString out = "embed_dd_r16.root",
              TString simReco = "/mnt/f/a1975/sim_rebuild/r280/dd_el/el_s6001_reco.root", int nSimMax = 100000,
              TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   TFile fs(simReco);
   auto ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *hs = nullptr;
   ts->SetBranchStatus("*", 0);
   ts->SetBranchStatus("AtEventH*", 1);
   ts->SetBranchAddress("AtEventH", &hs);
   TFile fd(Form("/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root", dataRun));
   auto td = (TTree *)fd.Get("cbmsim");
   TClonesArray *hd = nullptr;
   td->SetBranchStatus("*", 0);
   td->SetBranchStatus("AtEventH*", 1);
   td->SetBranchAddress("AtEventH", &hd);
   const Long64_t nD = td->GetEntries();
   std::vector<float> icv; // IC of the host event (to split carbon / oxygen hosts in the analysis)
   if (TFile *fi = TFile::Open(icDir + Form("run_%04d_IC.root", dataRun))) {
      if (!fi->IsZombie()) {
         auto *ti = (TTree *)fi->Get("ic");
         float icm = 0;
         ti->SetBranchAddress("icmax", &icm);
         for (Long64_t j = 0; j < ti->GetEntries(); ++j) {
            ti->GetEntry(j);
            icv.push_back(icm);
         }
      }
      fi->Close();
   }

   // --- the chain, configured as embed_vz_dp.C (= unpackReco_multifit.C data arguments: mcs 0, ms 3)
   AtTools::DataCleaning::AtDirDeDxCleaner cleaner;
   AtPATTERN::AtTrackFinderHDBSCAN hdb;
   hdb.SetMinClusterSize(0);
   hdb.SetMinSamples(3);
   hdb.SetClusterSelectionEpsilon(10.0);
   hdb.SetJoinMethod("mover");
   hdb.SetMinClusterSizeJoin(15);
   hdb.SetCircleOverlapRatio(0.25);
   hdb.SetMotionGapTol(40);
   hdb.SetMotionAngleTol(35);
   hdb.SetDiffusionParams(0.0009, 0.0009, 1.10424, 0.16, -1);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TFile fo(out, "RECREATE");
   TTree o("e", "elastic deuteron embedding");
   int sev, dev, mode, found, ref, nhS, nhD, ntr;
   double tc0, vz0, ke0, tcB, vzB, keB;
   float icD;
   o.Branch("sev", &sev);
   o.Branch("dev", &dev);
   o.Branch("mode", &mode);
   o.Branch("found", &found);
   o.Branch("ref", &ref);
   o.Branch("tc0", &tc0);
   o.Branch("vz0", &vz0);
   o.Branch("ke0", &ke0);
   o.Branch("tcB", &tcB);
   o.Branch("vzB", &vzB);
   o.Branch("keB", &keB);
   o.Branch("nhS", &nhS);
   o.Branch("nhD", &nhD);
   o.Branch("ntr", &ntr);
   o.Branch("icD", &icD);

   // returns the best locus track; with a reference, the one consistent with it
   auto run = [&](std::vector<const AtHit *> &in, bool useRef) {
      found = 0;
      ntr = 0;
      tcB = vzB = keB = -999;
      AtEvent ev;
      int id = 0;
      for (auto *h : in) {
         auto &n = ev.AddHit(*h);
         n.SetHitID(id++);
      }
      auto cl = cleaner.CleanData(ev.GetHits());
      AtEvent evc;
      id = 0;
      for (auto &h : cl) {
         h->SetHitID(id++);
         evc.AddHit(std::move(h));
      }
      if (evc.GetHits().size() <= 10 || evc.GetHits().size() >= 5000)
         return;
      std::unique_ptr<AtPatternEvent> pe;
      try {
         pe = hdb.FindTracks(evc);
      } catch (...) {
         return;
      }
      if (!pe)
         return;
      double bestd = 1e9;
      for (auto &tr : pe->GetTrackCand()) {
         ++ntr;
         auto s = spyKin(spy.Estimate(tr));
         if (!s.ok)
            continue;
         int nd = 0;
         for (auto &h : tr.GetHitArray())
            if (isDeut(*h))
               ++nd;
         if (!useRef) { // reference: the locus track holding the most deuteron hits
            if (nd > -bestd) {
               bestd = -nd;
               found = 1;
               tcB = s.tc;
               vzB = s.vz;
               keB = s.ke;
            }
            continue;
         }
         if (std::abs(s.tc - tc0) < 3 && std::abs(s.vz - vz0) < 40 && std::abs(s.ke / ke0 - 1) < 0.15) {
            double d = std::abs(s.tc - tc0) / 3 + std::abs(s.vz - vz0) / 40;
            if (d < bestd) {
               bestd = d;
               found = 1;
               tcB = s.tc;
               vzB = s.vz;
               keB = s.ke;
            }
         }
      }
   };

   TRandom3 rng(4242 + dataRun);
   int k = 0, nref = 0;
   const Long64_t nS = std::min<Long64_t>(ts->GetEntries(), nSimMax);
   for (Long64_t i = 0; i < nS; ++i) {
      ts->GetEntry(i);
      auto *es = (AtEvent *)hs->At(0);
      if (!es)
         continue;
      std::vector<const AtHit *> deut;
      for (auto &h : es->GetHits())
         if (isDeut(*h))
            deut.push_back(h.get());
      if (deut.size() < 5)
         continue; // beam-only event or a deuteron with no hits
      sev = i;
      nhS = deut.size();
      dev = -1;
      nhD = 0;
      icD = -1;
      mode = 1;
      run(deut, false);
      ref = found;
      tc0 = tcB;
      vz0 = vzB;
      ke0 = keB;
      o.Fill();
      ++k;
      if (!ref)
         continue;
      ++nref;
      for (int j = 0; j < nOver; ++j) {
         AtEvent *ed = nullptr;
         for (int tries = 0; tries < 50 && !ed; ++tries) {
            dev = rng.Integer(nD);
            td->GetEntry(dev);
            ed = (AtEvent *)hd->At(0);
            if (ed && ed->GetHits().empty())
               ed = nullptr;
         }
         if (!ed)
            continue;
         icD = dev < (Long64_t)icv.size() ? icv[dev] : -1.f;
         std::vector<const AtHit *> a = deut;
         for (auto &h : ed->GetHits())
            a.push_back(h.get());
         nhD = ed->GetHits().size();
         mode = 2;
         run(a, true);
         o.Fill();
      }
      if (k % 100 == 0) {
         printf("  %d sim deuterons, %d with a reference\n", k, nref);
         fflush(stdout);
      }
   }
   fo.cd();
   o.Write();
   fo.Close();
   printf("done -> %s (%d sim deuterons, %d with a reference)\n", out.Data(), k, nref);
}
