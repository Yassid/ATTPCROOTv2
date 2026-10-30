/// @file dir3d_bench.C
/// @brief Run the 3D direction-free estimator (helix3d.h) on every pattern track, alongside the two
///        production direction calls (AtPRA GeoTheta, AtSpyralPID), into one flat tree per source.
///
/// Same tree for data and simulation so one analysis macro scores both. Simulation adds the MC
/// truth of the generated ejectile (start vertex, KE, theta), which referees WHICH END IS THE
/// VERTEX directly. Nothing about frames is assumed here: hit z is stored in the DIGI frame and
/// both z conventions are scored offline (see helix3d.h).
///
///   data: root -l -b -q 'pid/dir3d_bench.C("/mnt/f/a1975/reco/run_0106_reco.root","d3_h2_0106.root",0,106,8000)'
///   sim : root -l -b -q 'pid/dir3d_bench.C("<reco>","d3_sim.root",1,0,6000,"<sim.root>",2212)'
#include "helix3d.h"

void dir3d_bench(TString in, TString out, int src = 0, int run = 0, Long64_t maxEv = 8000, TString simFile = "",
                 int truthPdg = 2212, double B = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;

   TFile *f = TFile::Open(in);
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *pa = nullptr, *mc = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   TTree *ts = nullptr;
   if (src == 1) {
      auto *fs = TFile::Open(simFile);
      ts = (TTree *)fs->Get("cbmsim");
      ts->SetBranchStatus("*", 0);
      ts->SetBranchStatus("MCTrack*", 1);
      ts->SetBranchAddress("MCTrack", &mc);
   }
   AtTools::AtSpyralPID spy;
   spy.SetBField(B);

   TFile fo(out, "RECREATE");
   TTree o("d3", "3D direction-free estimator bench");
   int S = src, RUN = run, tid, n, nb, ok, spyV, spyD, spyF;
   Long64_t ev;
   double len, rhoMax, geoTh, spyPol, spyBr, spySq, tx, ty, tz, tke, tth;
   // per end k = 0 (backbone A), 1 (B)
   double ex[2], ey[2], ezd[2], dca[2], gap[2], R[2], th[2], br[2], qd[2], qe[2], vz[2], dph[2];
   int sen[2], nu[2];
   o.Branch("src", &S); o.Branch("run", &RUN); o.Branch("ev", &ev); o.Branch("tid", &tid); o.Branch("n", &n);
   o.Branch("nb", &nb); o.Branch("ok", &ok); o.Branch("len", &len); o.Branch("rhoMax", &rhoMax);
   o.Branch("geoTh", &geoTh); o.Branch("spyV", &spyV); o.Branch("spyF", &spyF); o.Branch("spyD", &spyD);
   o.Branch("spyPol", &spyPol); o.Branch("spyBr", &spyBr); o.Branch("spySq", &spySq);
   o.Branch("tx", &tx); o.Branch("ty", &ty); o.Branch("tz", &tz); o.Branch("tke", &tke); o.Branch("tth", &tth);
   o.Branch("ex", ex, "ex[2]/D"); o.Branch("ey", ey, "ey[2]/D"); o.Branch("ezd", ezd, "ezd[2]/D");
   o.Branch("dca", dca, "dca[2]/D"); o.Branch("gap", gap, "gap[2]/D"); o.Branch("R", R, "R[2]/D");
   o.Branch("th", th, "th[2]/D"); o.Branch("br", br, "br[2]/D"); o.Branch("qd", qd, "qd[2]/D");
   o.Branch("qe", qe, "qe[2]/D"); o.Branch("vz", vz, "vz[2]/D"); o.Branch("dph", dph, "dph[2]/D");
   o.Branch("sen", sen, "sen[2]/I"); o.Branch("nu", nu, "nu[2]/I");

   Long64_t N = std::min(maxEv, t->GetEntries());
   for (ev = 0; ev < N; ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || pa->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)pa->At(0);
      if (!pe) continue;
      tx = ty = tz = tke = tth = -999;
      if (ts) {
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); ++j) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (!m || m->GetPdgCode() != truthPdg) continue;
            double px = m->GetPx(), py = m->GetPy(), pz = m->GetPz(), pm = std::sqrt(px * px + py * py + pz * pz);
            double M = m->GetMass() * 1000.0;
            tke = std::sqrt(pm * pm * 1e6 + M * M) - M;
            tth = std::acos(pz / pm) * TMath::RadToDeg();
            tx = m->GetStartX() * 10; ty = m->GetStartY() * 10; tz = m->GetStartZ() * 10;
            break;
         }
      }
      for (auto &tr : pe->GetTrackCand()) {
         auto &hits = tr.GetHitArray();
         n = hits.size();
         if (n < 10) continue;
         tid = tr.GetTrackID();
         geoTh = tr.GetGeoTheta() * TMath::RadToDeg();
         std::vector<h3d::P> h;
         h.reserve(n);
         rhoMax = 0;
         for (auto &x : hits) {
            auto p = x->GetPosition();
            h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()});
            rhoMax = std::max(rhoMax, std::hypot(p.X(), p.Y()));
         }
         auto r = spy.Estimate(tr);
         spyV = r.valid; spyF = r.failCode; spyD = r.direction;
         spyPol = r.valid ? r.polar * TMath::RadToDeg() : -999;
         spyBr = r.valid ? r.brho : -999; spySq = r.valid ? r.sqrtdEdx : -999;

         // digi-frame z (zSign +1, zOff 0) so the stored geometry is raw; theta/vz then refer to
         // +z_digi and are converted offline under whichever convention truth selects
         auto res = h3d::estimate(h, B, +1, 0);
         ok = res.ok; nb = res.nBackbone; len = res.length;
         for (int k = 0; k < 2; ++k) {
            auto &E = res.e[k];
            int idx = res.order.empty() ? 0 : (k == 0 ? res.order.front() : res.order.back());
            ex[k] = h[idx].x; ey[k] = h[idx].y; ezd[k] = h[idx].z;
            dca[k] = E.dca; gap[k] = E.gap; R[k] = E.R; th[k] = E.theta * TMath::RadToDeg(); br[k] = E.brho;
            qd[k] = E.qdens; qe[k] = E.qEdge; vz[k] = E.vz; dph[k] = E.dphi; sen[k] = E.sense; nu[k] = E.nUsed;
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
