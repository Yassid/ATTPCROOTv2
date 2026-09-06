/// @file fit_one_dp.C
/// @brief Refit ONE track with genfit, standalone, so fitter options can be tried in seconds
///        instead of re-running a 400 MB run or a whole reco.
///
/// WHY THIS IS POSSIBLE AT ALL: the genfitter output files keep the PATTERN tracks
/// (AtTrackingEvent::GetTrackArray) with their hit clusters, and AtGenfitter::Init() needs only a
/// loaded TGeoManager -- no FairRunAna, no reco file, no par file (cluster positions are already
/// in the file and ZPadPlane is set explicitly). So one track can be refitted in about a second.
///
/// THE QUESTION IT WAS BUILT FOR. AtGenfitter orders clusters by the drift coordinate z
/// (AtGenfitter.cxx:230). That is correct only while the helix advances in z faster than the
/// z resolution. On run_0020 entry 9789 the advance is 1.42 mm per cluster against 1.84 mm for one
/// time bucket (dv 1.10424 cm/us, 6 MHz), so the sort orders NOISE: the z-sorted sequence is
/// 2.15x longer than the stored order and 4.3x longer than the 1518 mm helix the proton actually
/// flew. The Kalman filter then transports the state and integrates CATIMA energy loss along a
/// path the particle never took, and converges on the wrong helix (R 53.6 mm instead of 78 mm at
/// the vertex, KE 0.385 MeV instead of ~2.6). Measured over 150 tracks: where z-sorting preserves
/// the path (inflation < 1.05) the median chi2/ndf is 0.01; where it scrambles it (> 1.5) it is 22.
///
///   orderMode 0  z-sort, exactly what the production does          (the control)
///             1  keep the stored cluster order (SetUseClusterOrder)
///             2  arc walk: greedy nearest-neighbour from an end     (no re-reco needed)
///
/// Defaults reproduce prod_dp_catima.sh exactly, so orderMode 0 must return the value already in
/// the cache. If it does not, something else changed and the comparison is void -- check that
/// first, before believing any improvement from mode 1 or 2.
///
///   root -b -q 'fit_one_dp.C("run_0020_multifit",9789,0,0)'    // control
///   root -b -q 'fit_one_dp.C("run_0020_multifit",9789,0,2)'    // arc walk

#include <algorithm>
#include <vector>

/// Greedy nearest-neighbour walk. Starts from the cluster farthest from the one farthest from the
/// centroid -- a cheap, deterministic way to pick a genuine END of the track rather than a point
/// in the middle, which is what makes the walk follow the trajectory instead of cutting across it.
static std::vector<int> dpArcWalk(const std::vector<AtHitCluster> &hc)
{
   const int n = hc.size();
   std::vector<ROOT::Math::XYZPoint> p(n);
   for (int i = 0; i < n; ++i) p[i] = hc[i].GetPosition();
   double cx = 0, cy = 0, cz = 0;
   for (auto &q : p) { cx += q.X(); cy += q.Y(); cz += q.Z(); }
   cx /= n; cy /= n; cz /= n;
   auto d2 = [&](int a, int b) {
      return (p[a].X()-p[b].X())*(p[a].X()-p[b].X()) + (p[a].Y()-p[b].Y())*(p[a].Y()-p[b].Y())
           + (p[a].Z()-p[b].Z())*(p[a].Z()-p[b].Z());
   };
   int far = 0; double best = -1;
   for (int i = 0; i < n; ++i) {
      double d = (p[i].X()-cx)*(p[i].X()-cx)+(p[i].Y()-cy)*(p[i].Y()-cy)+(p[i].Z()-cz)*(p[i].Z()-cz);
      if (d > best) { best = d; far = i; }
   }
   int start = far; best = -1;
   for (int i = 0; i < n; ++i) if (d2(far, i) > best) { best = d2(far, i); start = i; }
   std::vector<char> used(n, 0);
   std::vector<int> out;
   int cur = start;
   used[cur] = 1; out.push_back(cur);
   for (int k = 1; k < n; ++k) {
      int nxt = -1; double bd = 1e30;
      for (int i = 0; i < n; ++i) if (!used[i] && d2(cur, i) < bd) { bd = d2(cur, i); nxt = i; }
      if (nxt < 0) break;
      used[nxt] = 1; out.push_back(nxt); cur = nxt;
   }
   return out;
}

static double dpPath(const std::vector<AtHitCluster> &hc)
{
   double s = 0;
   for (size_t i = 0; i + 1 < hc.size(); ++i) {
      auto a = hc[i].GetPosition(), b = hc[i + 1].GetPosition();
      s += std::sqrt((a.X()-b.X())*(a.X()-b.X())+(a.Y()-b.Y())*(a.Y()-b.Y())+(a.Z()-b.Z())*(a.Z()-b.Z()));
   }
   return s;
}

void fit_one_dp(TString runTag = "run_0020_multifit", Long64_t entry = 9789, int trackID = 0,
                int orderMode = 0, double measSigma = 4.0, bool matEffects = true,
                bool catimaOn = true, bool seedFromSpyral = false, bool backwardSeedFix = true,
                int minIter = 2, int maxIter = 5,
                TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                TString geoName = "ATTPC_D300torr_v2_geomanager.root",
                TString eLossTable = "proton_D2_300torr.txt", TString jsonOut = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf = TFile::Open(dir + "/geometry/" + geoName); gf->Get("FAIRGeom"); }

   TString fn = gfDir + runTag + "_genfitter_p.root";
   TFile *f = TFile::Open(fn);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", fn.Data()); return; }
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   t->GetEntry(entry);
   auto *ev = (AtTrackingEvent *)te->At(0);
   if (!ev) { printf("no AtTrackingEvent at entry %lld\n", entry); return; }
   auto trks = ev->GetTrackArray();
   AtTrack *src = nullptr;
   for (auto &tr : trks) if (tr.GetTrackID() == trackID) src = &tr;
   if (!src) { printf("no trackID %d at entry %lld\n", trackID, entry); return; }

   AtTrack track = *src;                       // copy: the cluster order is edited in place
   auto *hc = track.GetHitClusterArray();
   const double pathStored = dpPath(*hc);
   std::vector<AtHitCluster> orig = *hc;

   const char *modeName[] = {"z-sort (production)", "stored cluster order", "arc walk (kNN)"};
   if (orderMode == 2) {
      auto ord = dpArcWalk(orig);
      std::vector<AtHitCluster> re; re.reserve(ord.size());
      for (int i : ord) re.push_back(orig[i]);
      *hc = re;
   }
   // what z-sorting WOULD do, reported so the ordering is visible rather than assumed
   std::vector<AtHitCluster> zs = orig;
   std::sort(zs.begin(), zs.end(), [](const AtHitCluster &a, const AtHitCluster &b) {
      return (1000.0 - a.GetPosition().Z()) < (1000.0 - b.GetPosition().Z()); });

   printf("\n\033[1m=== %s entry %lld track %d : %s ===\033[0m\n", runTag.Data(), entry, trackID,
          modeName[orderMode < 0 || orderMode > 2 ? 0 : orderMode]);
   printf("  clusters %zu   PRA GeoTheta %.2f deg  GeoRadius %.1f mm\n",
          hc->size(), track.GetGeoTheta() * TMath::RadToDeg(), track.GetGeoRadius());
   printf("  path along  stored %.0f mm | z-sorted %.0f mm | this order %.0f mm\n",
          pathStored, dpPath(zs), dpPath(*hc));

   TString elossPath = eLossTable.Length() ? dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable : TString("");
   auto fitter = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 1.00782503207, 1,
                                                         elossPath.Data(), !matEffects, minIter, maxIter);
   if (catimaOn) { fitter->SetCatimaMaterial(kTRUE, kTRUE); fitter->SetCatimaELoss(kTRUE, kFALSE); }
   // MUST match fitGenfitter_a1975_deuterium.C: whenever an eloss table is passed the production
   // also turns on the HYBRID mode (table below beta*gamma = 0.05, Bethe-Bloch above). Leaving it
   // out was not a subtle difference -- the control came back at KE 4.517 MeV against the cache's
   // 0.385, with 131 fitted points instead of 289. Any setter missing here silently makes this
   // tool a DIFFERENT fitter, which is why orderMode 0 is checked against the cache every time.
   if (eLossTable.Length()) fitter->SetELossHybrid(kTRUE, 6.61e-5);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(measSigma);
   fitter->SetSeedFromSpyral(seedFromSpyral);
   fitter->SetMatEffectsFallback(kFALSE);
   fitter->SetThetaWindow(10.0, 170.0);
   fitter->SetBackwardSeedFix(backwardSeedFix);
   fitter->SetBackExtrapToAxis(kTRUE);
   fitter->SetUseClusterOrder(orderMode != 0);   // 0 lets AtGenfitter do its own z sort
   fitter->Init();

   // GetFittedTrack is protected; FitEvent is the public entry point, so the single track is
   // handed over as a one-track pattern event. That also keeps this tool on EXACTLY the code path
   // the production uses, rather than a parallel one that could diverge from it silently.
   AtPatternEvent pe;
   pe.AddTrack(track);
   AtTrackingEvent tev;
   fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
   if (tev.GetFittedTracks().empty()) {
      printf("  \033[1;31mno fitted track returned (gated out, or the fit threw)\033[0m\n");
      return;
   }
   auto *ft = tev.GetFittedTracks().front().get();
   auto &k = ft->GetKinematicsXtr();
   auto &kf = ft->GetKinematics();
   double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
   auto v = ft->GetVertex();
   auto &sp = ft->GetSmoothedPositions();
   // hit-to-fit rms, with the fit z brought back to the cluster frame (z_fit + z_cluster = 1000)
   double s2 = 0;
   for (auto &c : *hc) {
      auto q = c.GetPosition(); double bd = 1e30;
      for (auto &m : sp) {
         double dz = q.Z() - (1000.0 - m.Z());
         double d = (q.X()-m.X())*(q.X()-m.X()) + (q.Y()-m.Y())*(q.Y()-m.Y()) + dz*dz;
         if (d < bd) bd = d;
      }
      s2 += bd;
   }
   printf("  \033[1mKE_xtr %8.3f MeV   theta %6.2f deg   chi2/ndf %8.2f   ndf %.0f\033[0m\n",
          k.kineticEnergy, k.theta * TMath::RadToDeg(), ndf > 0 ? c2 / ndf : -1, ndf);
   printf("  KE_fit %8.3f MeV   theta %6.2f deg   vertex (%.1f, %.1f, %.1f)\n",
          kf.kineticEnergy, kf.theta * TMath::RadToDeg(), v.X(), v.Y(), v.Z());
   printf("  fit points %zu / %zu clusters    hit-fit rms %.1f mm\n",
          sp.size(), hc->size(), std::sqrt(s2 / hc->size()));
   if (jsonOut.Length()) {
      std::ofstream o(jsonOut.Data());
      o << "{\"mode\":\"" << modeName[orderMode] << "\",\"ke\":" << k.kineticEnergy
        << ",\"theta\":" << k.theta * TMath::RadToDeg() << ",\"chi2ndf\":" << (ndf > 0 ? c2 / ndf : -1)
        << ",\"fit\":[";
      for (size_t i = 0; i < sp.size(); ++i)
        o << (i ? "," : "") << "[" << sp[i].X() << "," << sp[i].Y() << "," << 1000.0 - sp[i].Z() << "]";
      o << "],\"hits\":[";
      for (size_t i = 0; i < hc->size(); ++i) { auto q = (*hc)[i].GetPosition();
        o << (i ? "," : "") << "[" << q.X() << "," << q.Y() << "," << q.Z() << "," << (*hc)[i].GetCharge() << "]"; }
      o << "]}";
      printf("  wrote %s\n", jsonOut.Data());
   }
}
