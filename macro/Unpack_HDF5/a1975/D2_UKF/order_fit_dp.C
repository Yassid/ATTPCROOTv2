/// @file order_fit_dp.C
/// @brief Z SCRAMBLING, tested against truth (2026-11-02): does the order the clusters reach genfit in matter, and is the
///        3D backbone (AtHelix3D, no z sort, no DAQ order) better than the z sort? One row per track with >= 30 hits (sim:
///        events with exactly one such track and a truth particle of pdg; rows align with pid_sbpid_dp.C), and per mode the
///        fit of the SAME track with the SAME fitter settings as gf_hyp_run_dp.sh (CATIMA eloss full range + MSC + straggling):
///          mode 0  z sort, direction from GeoTheta                 the production
///          mode 1  z sort, direction from the 3D vote              separates direction from order
///          mode 2  stored cluster order, direction from GeoTheta   SetUseClusterOrder
///          mode 3  clusters sorted along the 3D backbone, 3D vote
///          mode 4  HITS sorted along the backbone, RE-CLUSTERED (ClusterizeSmooth3D), 3D vote: the test named on
///                  2026-09-06 and never run -- are the clusters already broken when they are built?
///          mode 5  as 4 but direction from GeoTheta (separates the new clusters from the 3D direction)
///          mode 6  as 5 but the fitter z-sorts the new clusters (does the z sort still hurt with good clusters?)
///          mode 7  as 5 but keeping the hits the backbone drops as blobs, placed by their nearest centroid
///          mode 8  ONLY the blob removal: the hits AtHelix3D drops are removed, the rest keep their stored order, and clustering,
///                  z sort and direction are the production's (no ordering change at all: the single-variable test)
///          mode 9  as 8 with the LARGE-PAD GUARD: dropped hits beyond 152 mm from the axis are kept (the larger pad spacing opens gaps
///                  that split the backbone, so on sim the outer part of real tracks was being dropped)
///        modeMask selects the modes to run (bit m = mode m).
///        Scrambling metrics per track: path of the clusters in z / stored / backbone order, net-vs-summed phase.
///   root -l -b -q 'order_fit_dp.C("<reco>","<sim>","out.root",1000010030,3.01550072,"t")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/AtTools/AtPID/AtHelix3DCore.h"
using XYZPoint = ROOT::Math::XYZPoint;
namespace ordfit {
double path(const std::vector<XYZPoint> &v)
{
   double s = 0;
   for (size_t i = 1; i < v.size(); ++i) s += std::sqrt((v[i] - v[i - 1]).Mag2());
   return s;
}
/// backbone coordinate of a position (digi frame), from the vertex end
double sOf(const AtTools::Helix3D::Result &r, int vEnd, const XYZPoint &q)
{
   double bd = 1e30;
   int bj = 0;
   for (size_t j = 0; j < r.cen.size(); ++j) {
      const auto &c = r.cen[j];
      double d = (q.X() - c.x) * (q.X() - c.x) + (q.Y() - c.y) * (q.Y() - c.y) + (1000 - q.Z() - c.z) * (1000 - q.Z() - c.z);
      if (d < bd) { bd = d; bj = j; }
   }
   return vEnd == 0 ? r.sc[bj] : r.length - r.sc[bj];
}
bool zLess(const XYZPoint &a, const XYZPoint &b) { return a.Z() < b.Z(); }
void loadEvTid(TString file, std::set<std::pair<Long64_t, int>> &s)
{
   std::ifstream il(file.Data());
   Long64_t e;
   int t;
   while (il >> e >> t) s.insert({e, t});
}
void loadRows(TString file, std::set<Long64_t> &rows)
{
   std::ifstream il(file.Data());
   Long64_t v;
   while (il >> v) rows.insert(v);
}
} // namespace ordfit
void order_fit_dp(TString in, TString simFile, TString out, int pdg, double massAmu, TString tag, Long64_t firstEv = 0,
                  Long64_t maxEv = 1000000000, double measSigma = 4.0, double maxTrueBrho = 1e9, int modeMask = 31,
                  TString rowList = "", Long64_t rowBase = 0, int rowMod = 1, int rowRem = 0, TString evTidList = "")
{
   // evTidList: text file of "event trackID" pairs to fit: survives a change of pattern recognition, where row numbers shift
   std::set<std::pair<Long64_t, int>> wantedET;
   if (evTidList.Length()) ordfit::loadEvTid(evTidList, wantedET);
   // rowList: text file of pid_sbpid_dp.C row numbers to fit (e.g. the proton band); rowBase = merged-tree row of the first track of
   // this run. With it the whole file is scanned from event 0 to count rows, and only listed rows with row % rowMod == rowRem are fitted.
   std::set<Long64_t> wanted;
   if (rowList.Length()) {
      ordfit::loadRows(rowList, wanted);
      firstEv = 0;
      maxEv = 1000000000;
   }
   Long64_t rowCount = rowBase;
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   namespace HX = AtTools::Helix3D;
   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TFile gf(dir + "/geometry/ATTPC_D300torr_v2_geomanager.root");
   gf.Get("FAIRGeom");
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   TFile *fs = simFile.Length() ? TFile::Open(simFile) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   TClonesArray *mc = nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   auto fitter = std::make_unique<EventFit::AtGenfitter>(-2.85, pdg, massAmu, 1, "", kFALSE, 2, 5);
   fitter->SetCatimaMaterial(kTRUE, kTRUE);
   fitter->SetCatimaELoss(kTRUE, kTRUE);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetMatEffectsFallback(kFALSE);
   fitter->SetThetaWindow(10.0, 170.0);
   fitter->SetBackwardSeedFix(kTRUE);
   fitter->SetBackExtrapToAxis(kTRUE);
   fitter->Init();
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtTrackTransformer tf;
   const int NM = 10;
   TFile fo(out, "RECREATE");
   TTree o("o", "order test");
   double tbr, tth, pz, ps, pb, lenB, th3, b3;
   int nh, ncl, ncl4, h3ok, bw3, bwG, skip, nDrop, tid;
   Long64_t rowId, evId;
   double vzFit[10];
   int ok[NM];
   double ke[NM], th[NM], c2[NM], nd[NM], vr[NM];
   o.Branch("tbr", &tbr); o.Branch("tth", &tth); o.Branch("nh", &nh); o.Branch("ncl", &ncl); o.Branch("ncl4", &ncl4);
   o.Branch("pz", &pz); o.Branch("ps", &ps); o.Branch("pb", &pb); o.Branch("lenB", &lenB); o.Branch("th3", &th3); o.Branch("b3", &b3);
   o.Branch("h3ok", &h3ok); o.Branch("bw3", &bw3); o.Branch("bwG", &bwG); o.Branch("skip", &skip); o.Branch("nDrop", &nDrop);
   o.Branch("row", &rowId); o.Branch("ev", &evId); o.Branch("tid", &tid);
   for (int m = 0; m < NM; ++m) {
      o.Branch(Form("ok%d", m), &ok[m]); o.Branch(Form("ke%d", m), &ke[m]); o.Branch(Form("th%d", m), &th[m]);
      o.Branch(Form("c2%d", m), &c2[m]); o.Branch(Form("nd%d", m), &nd[m]); o.Branch(Form("vr%d", m), &vr[m]); o.Branch(Form("vz%d", m), &vzFit[m]);
   }
   long nDone = 0;
   for (Long64_t ev = firstEv; ev < std::min(firstEv + maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto cands = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      tbr = tth = -1;
      if (ts) {
         int nl = 0;
         for (auto &tr : cands) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (!m || m->GetPdgCode() != pdg) continue;
            const double p = std::sqrt(m->GetPx() * m->GetPx() + m->GetPy() * m->GetPy() + m->GetPz() * m->GetPz());
            tbr = p / 0.299792458;
            tth = std::acos(m->GetPz() / p);
            break;
         }
         if (tbr <= 0) continue;
      }
      for (auto &tr0 : cands) {
         if (tr0.GetHitArray().size() < 30) continue;
         rowId = rowCount++;
         if (!wanted.empty() && (!wanted.count(rowId) || rowId % rowMod != rowRem)) continue;
         evId = ev; tid = tr0.GetTrackID();
         if (!wantedET.empty() && !wantedET.count({evId, tid})) continue;
         nh = tr0.GetHitArray().size();
         for (int m = 0; m < NM; ++m) { ok[m] = 0; ke[m] = th[m] = c2[m] = nd[m] = vr[m] = -1; vzFit[m] = -9999; }
         pz = ps = pb = lenB = th3 = b3 = -1; h3ok = bw3 = 0; ncl = ncl4 = 0; skip = 0; nDrop = -1;
         bwG = tr0.GetGeoTheta() > M_PI / 2;
         if (tbr > maxTrueBrho) { skip = 1; o.Fill(); continue; }
         // the 3D backbone and vote (default options, as AtHelix3D::Estimate)
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr0));
         std::vector<HX::P> h;
         for (auto &x : tr0.GetHitArray()) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         HX::EndOpt eo;
         auto r = HX::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, eo);
         int vEnd = -1;
         if (r.ok && r.cen.size() >= 2) {
            int spyEnd = -1;
            if (sr.valid && (sr.direction == 0 || sr.direction == 1))
               spyEnd = (sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1;
            auto V = HX::majorityEnd(r, spyEnd, 1, 0.2, 152.0, 150.0);
            vEnd = V.end;
            nDrop = (int)h.size() - (int)r.order.size();
            h3ok = !V.ambiguous;
            const auto &E = r.e[V.end];
            th3 = E.theta; b3 = E.brho; lenB = r.length;
            const double zV = 1000 - h[V.end == 0 ? r.order.front() : r.order.back()].z, zF = 1000 - h[V.end == 0 ? r.order.back() : r.order.front()].z;
            bw3 = zV > zF;
         }
         auto *hc0 = const_cast<AtTrack &>(tr0).GetHitClusterArray();
         if (!hc0 || hc0->size() < 3) { o.Fill(); continue; }
         ncl = hc0->size();
         std::vector<XYZPoint> vs, vz;
         for (auto &c : *hc0) vs.push_back(c.GetPosition());
         vz = vs;
         std::sort(vz.begin(), vz.end(), ordfit::zLess);
         ps = ordfit::path(vs); pz = ordfit::path(vz);
         for (int m = 0; m < NM; ++m) {
            if (!((modeMask >> m) & 1)) continue;
            if (m != 0 && m != 2 && vEnd < 0) continue;
            AtTrack tc(tr0);
            if (m == 1 || m == 3 || m == 4) tc.SetGeoTheta(bw3 ? 3.0 : 0.5); // only its side of 90 deg is used by the fitter
            auto *hc = tc.GetHitClusterArray();
            if (m >= 3) {
               std::vector<std::pair<double, int>> idx;
               if (m >= 4) {
                  // hits along the backbone (AtHelix3D's order, from the vertex end), blob hits dropped
                  std::vector<std::pair<double, int>> oh;
                  if (m == 8) {
                     for (size_t i = 0; i < r.order.size(); ++i) oh.push_back({(double)r.order[i], r.order[i]}); // stored order
                  } else if (m == 9) {
                     std::set<int> kept(r.order.begin(), r.order.end());
                     for (size_t i = 0; i < h.size(); ++i)
                        if (kept.count(i) || std::hypot(h[i].x, h[i].y) > 152.0) oh.push_back({(double)i, (int)i});
                  } else if (m == 7) {
                     for (size_t i = 0; i < h.size(); ++i)
                        oh.push_back({ordfit::sOf(r, vEnd, tr0.GetHitArray()[i]->GetPosition()), (int)i});
                  } else {
                     for (size_t i = 0; i < r.order.size(); ++i) oh.push_back({vEnd == 0 ? r.s[i] : r.length - r.s[i], r.order[i]});
                  }
                  std::sort(oh.begin(), oh.end());
                  AtTrack t4;
                  t4.SetTrackID(tr0.GetTrackID()); t4.SetGeoTheta(tc.GetGeoTheta()); t4.SetGeoRadius(tc.GetGeoRadius());
                  t4.SetGeoCenter(tc.GetGeoCenter()); t4.SetGeoPhi(tc.GetGeoPhi());
                  auto &H0 = tr0.GetHitArray();
                  for (auto &pr : oh) t4.AddHit(*H0[pr.second]);
                  tf.ClusterizeSmooth3D(t4, 15.0, 7.5);
                  tc = t4;
                  hc = tc.GetHitClusterArray();
                  ncl4 = hc->size();
                  if (ncl4 < 3) continue;
               }
               for (size_t i = 0; i < hc->size(); ++i) idx.push_back({m >= 8 ? (double)i : ordfit::sOf(r, vEnd, hc->at(i).GetPosition()), (int)i});
               std::sort(idx.begin(), idx.end());
               // the fitter's array runs low z_lab -> high z_lab for a forward seed and is reversed for a backward one:
               // front = vertex when forward, back = vertex when backward
               if (bw3 && m < 8) std::reverse(idx.begin(), idx.end());
               // modes 5-7 take the direction from GeoTheta, so the array must simply run low z_lab -> high z_lab
               if (m >= 5 && m < 8 && hc->at(idx.front().second).GetPosition().Z() < hc->at(idx.back().second).GetPosition().Z())
                  std::reverse(idx.begin(), idx.end()); // z_digi falls as z_lab rises
               std::vector<AtHitCluster> re;
               for (auto &pr : idx) re.push_back(hc->at(pr.second));
               if (m == 3) { std::vector<XYZPoint> vb; for (auto &c : re) vb.push_back(c.GetPosition()); pb = ordfit::path(vb); }
               *hc = re;
            }
            fitter->SetUseClusterOrder(m >= 2 && m != 6 && m < 8);
            AtPatternEvent pe;
            pe.AddTrack(tc);
            AtTrackingEvent tev;
            fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
            if (tev.GetFittedTracks().empty()) continue;
            auto *ft = tev.GetFittedTracks().front().get();
            auto &k = ft->GetKinematicsXtr();
            ke[m] = k.kineticEnergy; th[m] = k.theta;
            c2[m] = ft->GetTrackMetadata()->GetChi2(); nd[m] = ft->GetTrackMetadata()->GetNdf();
            auto v = ft->GetVertex();
            vr[m] = std::hypot(v.X(), v.Y()); vzFit[m] = v.Z();
            ok[m] = ke[m] > 0 && nd[m] > 0;
         }
         o.Fill();
         if (++nDone % 200 == 0) { printf("%ld tracks, ev %lld\n", nDone, ev); fflush(stdout); }
      }
   }
   fo.cd(); o.Write();
   printf("done: %lld rows\n", o.GetEntries());
}
