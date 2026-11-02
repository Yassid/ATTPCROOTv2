/// @file range_dp.C
/// @brief RANGE per track for a range-vs-momentum PID test (2026-11-01). One row per track in the row order of
///        pid_sbpid_dp.C (sim mode: events with exactly one >= 30-hit candidate and a truth particle of pdg; data mode
///        when simFile is empty: every >= 30-hit track), so it is a friend of pv_*.root. Columns:
///          len   AtHelix3D backbone length of the whole track [mm] (3D centroid path)
///          lhel  the same path measured along the vertex-end helix where possible: sqrt(R^2+(dz/dphi)^2)*|dphi| summed
///                over consecutive centroids (helix-path is less noisy than the centroid polyline, cf. the dE/dx study)
///          cont  far (stopping) end inside the gas: r < rCont and 5 < z < 995 mm
///          rfar, zfar  far-end position; tbr truth Brho (sim)
///          fpk, tail, midpeak  charge profile along the track (10 mm bins, 3-bin smoothed, hits ordered from the vertex
///                end): peak position / track length, mean charge after peak+20 mm / peak, and the MIDPEAK flag
///                (fpk < 0.75, peak > 30 mm before the end, tail < 0.35) -- the particle stopped and the candidate
///                carries extra hits (classify_short_t_dp.C, same definition)
///        Vertex end from the same 2-of-3 vote as pid_sbpid_dp.C.
///   root -l -b -q 'range_dp.C("<reco>","<sim>","out.root",1000010030)'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
void range_dp(TString in, TString simFile, TString out, int pdg, double rCont = 145, Long64_t firstEv = 0,
              Long64_t maxEv = 1000000000)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   namespace HX = AtTools::Helix3D;
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   TFile *fs = simFile.Length() ? TFile::Open(simFile) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   TClonesArray *mc = nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   HX::EndOpt eo;
   eo.hitRefit = true; eo.hitRefitZ = true;
   TFile fo(out, "RECREATE");
   TTree o("r", "track range");
   double len, lhel, rfar, zfar, tbr, fpk, tail;
   int cont, ok, midpeak;
   o.Branch("ok", &ok); o.Branch("len", &len); o.Branch("lhel", &lhel); o.Branch("cont", &cont);
   o.Branch("rfar", &rfar); o.Branch("zfar", &zfar); o.Branch("tbr", &tbr);
   o.Branch("fpk", &fpk); o.Branch("tail", &tail); o.Branch("midpeak", &midpeak);
   for (Long64_t ev = firstEv; ev < std::min(firstEv + maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto &cands = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      tbr = -1;
      if (ts) {
         int nl = 0;
         for (auto &tr : cands) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (!m || m->GetPdgCode() != pdg) continue;
            tbr = std::sqrt(m->GetPx() * m->GetPx() + m->GetPy() * m->GetPy() + m->GetPz() * m->GetPz()) / 0.299792458;
            break;
         }
         if (tbr <= 0) continue;
      }
      for (auto &tr : cands) {
         if (tr.GetHitArray().size() < 30) continue;
         ok = cont = midpeak = 0; len = lhel = rfar = zfar = fpk = tail = -1;
         auto sr = spy.Estimate(tr);
         std::vector<HX::P> h;
         for (auto &x : tr.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
         if (r.ok) {
            int spyEnd = -1;
            if (sr.valid && (sr.direction == 0 || sr.direction == 1))
               spyEnd = (sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1;
            auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
            const int mC = r.cen.size();
            if (mC >= 2) {
               const auto &cf = r.cen[V.end == 0 ? mC - 1 : 0];
               rfar = std::hypot(cf.x, cf.y); zfar = cf.z;
               cont = rfar < rCont && cf.z > 5 && cf.z < 995;
               len = r.length;
               // path along a LOCAL helix: each centroid step projected on xy arc of a 3-point circle is overkill; use
               // the xy chord plus dz between consecutive centroids, smoothed over a 3-centroid window
               lhel = 0;
               for (int j = 2; j < mC; j += 2) {
                  const auto &a = r.cen[j - 2], &b = r.cen[j];
                  lhel += std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
               }
               // charge profile from the vertex end
               {
                  const int n = r.order.size();
                  const int nb = std::max(1, (int)std::ceil(r.length / 10.0));
                  std::vector<double> prof(nb, 0), sm(nb, 0);
                  for (int i = 0; i < n; ++i) {
                     double sv = V.end == 0 ? r.s[i] : r.length - r.s[i];
                     prof[std::min(nb - 1, std::max(0, (int)(sv / 10.0)))] += h[r.order[i]].q;
                  }
                  for (int b = 0; b < nb; ++b) {
                     double sum = 0; int c = 0;
                     for (int d = -1; d <= 1; ++d) if (b + d >= 0 && b + d < nb) { sum += prof[b + d]; ++c; }
                     sm[b] = sum / c;
                  }
                  int bp = std::max_element(sm.begin(), sm.end()) - sm.begin();
                  const double sPk = (bp + 0.5) * 10;
                  fpk = sPk / std::max(1.0, r.length);
                  double ta = 0; int tc = 0;
                  for (int b = bp + 2; b < nb; ++b) { ta += sm[b]; ++tc; }
                  tail = tc ? ta / tc / std::max(1e-9, sm[bp]) : 1;
                  midpeak = fpk < 0.75 && r.length - sPk > 30 && tail < 0.35;
               }
               ok = 1;
            }
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write();
   printf("done: %lld rows\n", o.GetEntries());
}
