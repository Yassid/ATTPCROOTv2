/// @file embed_dp.C -- EMBEDDING test for the data-only loss of forward-cm 16C(d,p) bound protons (2026-10-12).
/// census_arc_dp.C showed the missing protons never appear as vertex-arc pattern tracks in data, while Spyral passes
/// the ones that do. This asks whether the real DATA EVENT ENVIRONMENT (beam, noise, other tracks) is what kills them:
/// g.s. sim events that the sim chain reconstructs correctly (Spyral bound, theta_cm 10-40, vz 200-400) have their
/// PSA hits (AtEventH, i.e. before cleaning) overlaid on random real data events (run's AtEventH, same stage), and the
/// SAME chain as both productions is re-run: AtDirDeDxCleaner -> AtTrackFinderHDBSCAN (mcs 0, ms 3, mover) -> AtSpyralPID.
/// Data reco ran with space charge OFF (dv_dt.sh), so both clouds are in the same frame and stage.
///   mode 0  sim full event alone           (closure: must reproduce the sim chain)
///   mode 1  sim PROTON hits alone          (closure for the injection: MC A=1,Z=1 on the hits)
///   mode 2  sim proton hits + data event   <- the test
///   mode 3  sim full event + data event
///   mode 4  sim PROTON hits shifted in x,y by a beam-spot offset (offX,offY + Gaus(spotSig)), hits landing inside the
///           20 mm beam hole dropped -- the sim beam is a pencil, the data vertex spread is ~4.5 mm (vertex_cmp_dp.C).
///           nSpot draws per sim event; only when spotSig > 0. Approximate: no pad re-quantisation.
/// Success = a pattern track with Spyral valid, theta_lab 89-161, KE>1.2, Ex in the bound window, |thcm - thcm0| < 5,
/// |vz - vz0| < 40 (thcm0, vz0 = the sim chain's own Spyral values for that event).
/// Usage: root -l -b -q 'embed_dp.C(16,3,"embed_r16.root")'   (data run, overlays per sim event, output)
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
using namespace exsel;
double Eb(double vz) { return 184.25 - 0.011030 * (vz - 500); }

struct Spy {
   bool ok;
   double th, vz, ke, ex, tc;
};
Spy spyKin(const AtTools::AtSpyralResult &r)
{
   Spy s{false, -999, -999, -999, -999, -999};
   if (!r.valid) return s;
   const double mp = 938.272;
   s.th = 180 - r.polar * TMath::RadToDeg(); s.vz = 1000 - r.vertex.Z();
   double p = 299.792458 * r.brho;
   s.ke = (std::sqrt(p * p + mp * mp) - mp) / 1.023;
   if (!kine2b(Eb(s.vz), s.th * TMath::DegToRad(), s.ke, s.ex, s.tc)) return s;
   s.ok = s.th > 89 && s.th < 161 && s.ke > 1.2 && s.ex > -1.5 && s.ex < 0.729;
   return s;
}

void embed_dp(int dataRun = 16, int nOver = 3, TString out = "embed_r16.root", int nSimMax = 100000, double spotSig = 0,
              double offX = 0, double offY = 0, int nSpot = 3)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   // --- the sim events to embed: the sim chain's own Spyral-bound successes
   TFile fc("/mnt/f/a1975/census_arc/sim_gs_s4001.root");
   auto c = (TTree *)fc.Get("c");
   c->Draw("ev:stcm:svz", "sv&&sth>89&&sth<161&&ske>1.2&&sex>-1.5&&sex<0.729&&svz>=200&&svz<400&&stcm>=10&&stcm<40", "goff");
   std::vector<std::array<double, 3>> sel;
   std::set<int> seen;
   for (int i = 0; i < c->GetSelectedRows(); ++i)
      if (seen.insert(int(c->GetV1()[i])).second) sel.push_back({c->GetV1()[i], c->GetV2()[i], c->GetV3()[i]});
   if ((int)sel.size() > nSimMax) sel.resize(nSimMax);
   printf("%zu sim events to embed\n", sel.size());

   TFile fs("/mnt/f/a1975_C16_dp_sim/gs_s4001_reco_g75000.root");
   auto ts = (TTree *)fs.Get("cbmsim");
   TClonesArray *hs = nullptr;
   ts->SetBranchStatus("*", 0); ts->SetBranchStatus("AtEventH*", 1); ts->SetBranchAddress("AtEventH", &hs);
   TFile fd(Form("/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root", dataRun));
   auto td = (TTree *)fd.Get("cbmsim");
   TClonesArray *hd = nullptr;
   td->SetBranchStatus("*", 0); td->SetBranchStatus("AtEventH*", 1); td->SetBranchAddress("AtEventH", &hd);
   const Long64_t nD = td->GetEntries();

   // --- the chain, configured exactly as unpackReco_multifit.C / run_reco_C16dt.C (data args: mcs 0, ms 3)
   AtTools::DataCleaning::AtDirDeDxCleaner cleaner;
   AtPATTERN::AtTrackFinderHDBSCAN hdb;
   hdb.SetMinClusterSize(0); hdb.SetMinSamples(3); hdb.SetClusterSelectionEpsilon(10.0);
   hdb.SetJoinMethod("mover"); hdb.SetMinClusterSizeJoin(15); hdb.SetCircleOverlapRatio(0.25);
   hdb.SetMotionGapTol(40); hdb.SetMotionAngleTol(35);
   hdb.SetDiffusionParams(0.0009, 0.0009, 1.10424, 0.16, -1); // as AtPRAtask logs them from the dv1104 par
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TFile fo(out, "RECREATE");
   TTree o("e", "embedding");
   int sev, dev, mode, found, ntr, nhS, nhP, nhD, nhClean, nhBest;
   double tc0, vz0, tcB, vzB, exB, keB, dx = 0, dy = 0;
   o.Branch("sev", &sev); o.Branch("dev", &dev); o.Branch("mode", &mode); o.Branch("found", &found);
   o.Branch("ntr", &ntr); o.Branch("nhS", &nhS); o.Branch("nhP", &nhP); o.Branch("nhD", &nhD);
   o.Branch("nhClean", &nhClean); o.Branch("nhBest", &nhBest); o.Branch("tc0", &tc0); o.Branch("vz0", &vz0);
   o.Branch("tcB", &tcB); o.Branch("vzB", &vzB); o.Branch("exB", &exB); o.Branch("keB", &keB);
   o.Branch("dx", &dx); o.Branch("dy", &dy);

   auto run = [&](std::vector<const AtHit *> &in) {
      AtEvent ev;
      int id = 0;
      for (auto *h : in) { auto &n = ev.AddHit(*h); n.SetHitID(id++); }
      auto cl = cleaner.CleanData(ev.GetHits());
      nhClean = cl.size();
      AtEvent evc;
      id = 0;
      for (auto &h : cl) { h->SetHitID(id++); evc.AddHit(std::move(h)); }
      found = 0; ntr = 0; nhBest = 0; tcB = vzB = exB = keB = -999;
      if (nhClean <= 10 || nhClean >= 5000) return; // AtPRAtask fMinNumHits / fMaxNumHits
      std::unique_ptr<AtPatternEvent> pe;
      try { pe = hdb.FindTracks(evc); } catch (...) { return; }
      if (!pe) return;
      double bestd = 1e9;
      for (auto &tr : pe->GetTrackCand()) {
         ++ntr;
         auto s = spyKin(spy.Estimate(tr));
         if (!s.ok) continue;
         double d = std::abs(s.tc - tc0) / 5 + std::abs(s.vz - vz0) / 40;
         if (std::abs(s.tc - tc0) < 5 && std::abs(s.vz - vz0) < 40 && d < bestd) {
            bestd = d; found = 1; tcB = s.tc; vzB = s.vz; exB = s.ex; keB = s.ke; nhBest = tr.GetHitArray().size();
         }
      }
   };

   TRandom3 rng(12345);
   int k = 0;
   for (auto &e : sel) {
      sev = e[0]; tc0 = e[1]; vz0 = e[2];
      ts->GetEntry(sev);
      auto *es = (AtEvent *)hs->At(0);
      if (!es) continue;
      std::vector<const AtHit *> full, prot;
      for (auto &h : es->GetHits()) {
         full.push_back(h.get());
         auto &mc = h->GetMCSimPointArray();
         bool isP = false;
         for (auto &m : mc) if (m.A == 1 && m.Z == 1) isP = true;
         if (isP) prot.push_back(h.get());
      }
      nhS = full.size(); nhP = prot.size();
      dev = -1; nhD = 0;
      mode = 0; run(full); o.Fill();
      mode = 1; run(prot); o.Fill();
      for (int j = 0; spotSig > 0 && j < nSpot; ++j) {
         dx = offX + rng.Gaus(0, spotSig); dy = offY + rng.Gaus(0, spotSig);
         std::vector<AtHit> sh;
         sh.reserve(prot.size());
         for (auto *h : prot) {
            auto p = h->GetPosition();
            if (std::hypot(p.X() + dx, p.Y() + dy) < 20) continue;
            sh.push_back(*h);
            sh.back().SetPosition({p.X() + dx, p.Y() + dy, p.Z()});
         }
         std::vector<const AtHit *> sp;
         for (auto &h : sh) sp.push_back(&h);
         mode = 4; run(sp); o.Fill();
      }
      dx = dy = 0;
      for (int j = 0; j < nOver; ++j) {
         AtEvent *ed = nullptr;
         for (int tries = 0; tries < 50 && !ed; ++tries) {
            dev = rng.Integer(nD);
            td->GetEntry(dev);
            ed = (AtEvent *)hd->At(0);
            if (ed && ed->GetHits().empty()) ed = nullptr;
         }
         if (!ed) continue;
         std::vector<const AtHit *> dh;
         for (auto &h : ed->GetHits()) dh.push_back(h.get());
         nhD = dh.size();
         auto a = prot; a.insert(a.end(), dh.begin(), dh.end());
         mode = 2; run(a); o.Fill();
         auto b = full; b.insert(b.end(), dh.begin(), dh.end());
         mode = 3; run(b); o.Fill();
      }
      if (++k % 50 == 0) { printf("  %d / %zu\n", k, sel.size()); fflush(stdout); }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("done -> %s\n", out.Data());
}
