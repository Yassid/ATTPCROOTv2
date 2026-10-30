/// @file pid_width_dp.C
/// @brief Does the fitted pulse WIDTH improve the Spyral PID? Same tracks, same estimator, only the
///        per-hit charge differs.
///
/// AtSpyralPID sums AtHit::GetCharge() = the PEAK AMPLITUDE. For a fixed deposited charge a wider
/// pulse has a lower peak, and the width depends on how the track crosses the pad (dip angle,
/// diffusion), not only on dE/dx. AtPSAMultiFit's stored "integral" is A*int(response)*tau_NOMINAL,
/// i.e. amplitude times a constant, so it carries no width either. With SetFloatPeakingTime the
/// fitted width is stashed in AtHit::GetTimeStampCorrInter() (TB); this macro needs a reco made
/// that way (unpackReco_multifit.C with fitChi2 < 0, e.g. /mnt/f/a1975/reco_d2_tau/).
///
/// Per track, AtSpyralPID::Estimate on a copy of the track whose hit charges are:
///   amp  the peak amplitude (production)
///   wid  amplitude * tau_fit / tau_nominal   (proportional to the pulse integral)
/// and the 3D estimator (UKF/pid/helix3d.h): MST order, local helix at each end, vertex end from the
/// rotation sense (measured +1 from the vertex end in sim AND data), else the helix closest
/// approach, else the end nearer the axis. Brho from the vertex-end helix, dE/dx = charge per mm of
/// 3D helix path over that half-turn window, again with amp and wid charges.
/// Output tree "w": brho, sqrt(dE/dx) for all four, polar, hits, median tau, backbone spacing.
///
///   root -l -b -q 'pid_width_dp.C("/mnt/f/a1975/reco_d2_tau/run_0016_multifit_reco.root","/mnt/f/a1975/reco_d2_tau/pid_width_0016.root")'
#include "../UKF/pid/helix3d.h"

void pid_width_dp(TString in, TString out, Long64_t maxEv = 1000000000, double tauNomTB = 4.5)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TFile fo(out, "RECREATE");
   TTree o("w", "PID with peak amplitude vs width-corrected charge");
   int tid, nh, va, vw;
   Long64_t ev;
   double br, bw, sqa, sqw, pol, tauMed, tauRms, b3, sq3a, sq3w, th3, sp3;
   int v3, end3, spyD, how3, agreeEnd;
   double sWin3, sIn3, pathIn3, path3, R3, dca3, dcaO, R3o, zA, zB, sq3c, sq3cw, pathC, sq3f, sq3fw,
          spyDE, spyArc, pathF, hybB, hybBw, sq3t, sq3tw, pathH, sq3h, sq3hw, hybC, hybCw;
   o.Branch("ev", &ev); o.Branch("tid", &tid); o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("vw", &vw);
   o.Branch("br", &br); o.Branch("bw", &bw); o.Branch("sqa", &sqa); o.Branch("sqw", &sqw); o.Branch("pol", &pol);
   o.Branch("tauMed", &tauMed); o.Branch("tauRms", &tauRms);
   o.Branch("v3", &v3); o.Branch("end3", &end3); o.Branch("b3", &b3); o.Branch("sq3a", &sq3a);
   o.Branch("sq3w", &sq3w); o.Branch("th3", &th3); o.Branch("sp3", &sp3);
   o.Branch("spyD", &spyD); o.Branch("how3", &how3); o.Branch("agreeEnd", &agreeEnd); o.Branch("sWin3", &sWin3);
   o.Branch("sIn3", &sIn3); o.Branch("pathIn3", &pathIn3); o.Branch("path3", &path3); o.Branch("R3", &R3);
   o.Branch("sq3f", &sq3f); o.Branch("sq3fw", &sq3fw);
   o.Branch("spyDE", &spyDE); o.Branch("spyArc", &spyArc); o.Branch("pathF", &pathF); o.Branch("hybB", &hybB);
   o.Branch("pathH", &pathH); o.Branch("sq3h", &sq3h); o.Branch("sq3hw", &sq3hw); o.Branch("hybC", &hybC);
   o.Branch("hybCw", &hybCw);
   o.Branch("hybBw", &hybBw); o.Branch("sq3t", &sq3t); o.Branch("sq3tw", &sq3tw);
   o.Branch("sq3c", &sq3c); o.Branch("sq3cw", &sq3cw); o.Branch("pathC", &pathC);
   o.Branch("dca3", &dca3); o.Branch("dcaO", &dcaO); o.Branch("R3o", &R3o); o.Branch("zA", &zA); o.Branch("zB", &zB);
   for (ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         nh = H.size();
         if (nh < 10) continue;
         tid = tr.GetTrackID();
         std::vector<double> taus;
         AtTrack tw;
         tw.SetTrackID(tid);
         for (auto &h : H) {
            double tau = h->GetTimeStampCorrInter();
            taus.push_back(tau);
            AtHit c(*h);
            c.SetCharge(h->GetCharge() * tau / tauNomTB);
            tw.AddHit(c);
         }
         std::sort(taus.begin(), taus.end());
         tauMed = taus[taus.size() / 2];
         tauRms = 0.5 * (taus[taus.size() * 84 / 100] - taus[taus.size() * 16 / 100]);
         auto ra = spy.Estimate(tr);
         auto rw = spy.Estimate(tw);
         va = ra.valid; vw = rw.valid;
         br = ra.valid ? ra.brho : -1; bw = rw.valid ? rw.brho : -1;
         sqa = ra.valid ? ra.sqrtdEdx : -1; sqw = rw.valid ? rw.sqrtdEdx : -1;
         pol = ra.valid ? ra.polar * TMath::RadToDeg() : -999;
         spyD = ra.valid ? ra.direction : -1;
         spyDE = ra.valid ? ra.dE : -1; spyArc = ra.valid ? ra.arclength : -1;
         double spyDEw = rw.valid ? rw.dE : -1;

         // ---- 3D ----
         std::vector<h3d::P> h3;
         std::vector<double> qw;
         for (size_t i = 0; i < H.size(); ++i) {
            auto p = H[i]->GetPosition();
            h3.push_back({p.X(), p.Y(), p.Z(), H[i]->GetCharge()});
            qw.push_back(H[i]->GetCharge() * H[i]->GetTimeStampCorrInter() / tauNomTB);
         }
         auto r3 = h3d::estimate(h3, 2.85, -1, 1000);
         v3 = r3.ok; end3 = -1; b3 = sq3a = sq3w = th3 = -1;
         sq3c = sq3cw = pathC = sq3f = sq3fw = pathF = hybB = hybBw = sq3t = sq3tw = -1;
         pathH = sq3h = sq3hw = hybC = hybCw = -1;
         how3 = -1; agreeEnd = -1; sWin3 = sIn3 = pathIn3 = path3 = R3 = dca3 = dcaO = R3o = zA = zB = -1;
         sp3 = r3.nBackbone > 0 ? r3.length / r3.nBackbone : 1e9;
         if (r3.ok) {
            auto &E0 = r3.e[0], &E1 = r3.e[1];
            // vertex end by the 2-of-3 majority (sense, Spyral, Bragg). h3 is in the DIGI frame, so
            // Spyral forward (vertex at the lowest z_digi) maps to the end with the lower h3 z.
            int spyEnd = -1;
            if (spyD >= 0) {
               double zA0 = h3[r3.order.front()].z, zB0 = h3[r3.order.back()].z;
               spyEnd = (spyD == 0) == (zA0 < zB0) ? 0 : 1;
            }
            auto V = h3d::majorityEnd(r3, spyEnd);
            end3 = V.end;
            how3 = V.ambiguous ? 9 : V.nFor; // 3, 2, 1 voters agreeing; 9 = ambiguous
            (void)E0; (void)E1;
            auto &E = r3.e[end3];
            b3 = E.brho; th3 = E.theta * TMath::RadToDeg();
            zA = h3[r3.order.front()].z; zB = h3[r3.order.back()].z; // DIGI frame: h3 is not transformed, only the estimator copy is
            // Spyral's vertex end: forward (0) = lowest z_DIGI = highest z_lab
            if (spyD >= 0) { int spyEnd = (spyD == 0) == (zA < zB) ? 0 : 1; agreeEnd = spyEnd == end3; } // zA,zB digi
            sWin3 = E.sWin; path3 = E.path; R3 = E.R; dca3 = E.dca;
            dcaO = r3.e[1 - end3].dca; R3o = r3.e[1 - end3].R;
            double qa = 0, qww = 0;
            // hits inside the vertex-end window AND on the small inner pads (rho < 152 mm, as Spyral's
            // fSmallPadRadius): the outer pads are bigger and collect more charge per hit, so a window
            // straddling the boundary mixes two charge scales (seen as a split deuteron peak).
            double sIn = 0;
            for (size_t i = 0; i < r3.order.size(); ++i) {
               double sFromEnd = end3 == 0 ? r3.s[i] : r3.length - r3.s[i];
               auto &p = h3[r3.order[i]];
               if (sFromEnd > E.sWin || std::hypot(p.x, p.y) > 152.0) continue;
               qa += p.q; qww += qw[r3.order[i]];
               sIn = std::max(sIn, sFromEnd);
            }
            double pathIn = E.sWin > 0 ? sIn * E.path / E.sWin : 0; // backbone length -> helix path
            sIn3 = sIn; pathIn3 = pathIn;
            // CENTRELINE dE/dx: charge and 3D length of the centroid polyline over the SAME segment
            // (vertex-end window, small pads). Numerator and denominator from one definition.
            {
               const int m = r3.cen.size();
               // width-corrected charge per centroid: rebuild helix3d's bins (sequential, 4 mm in s)
               std::vector<double> qwc;
               for (size_t i = 0; i < r3.order.size();) {
                  size_t j = i;
                  double q = 0;
                  while (j < r3.order.size() && r3.s[j] - r3.s[i] < 4.0) { q += qw[r3.order[j]]; ++j; }
                  qwc.push_back(q);
                  i = j;
               }
               double qc = 0, lc = 0, qcw = 0;
               int prev = -1;
               for (int jj = 0; jj < m; ++jj) {
                  int j = end3 == 0 ? jj : m - 1 - jj;
                  double sFromEnd = end3 == 0 ? r3.sc[j] : r3.length - r3.sc[j];
                  auto &p = r3.cen[j];
                  if (sFromEnd > E.sWin) break;
                  if (std::hypot(p.x, p.y) > 152.0) { prev = -1; continue; }
                  qc += p.q;
                  if ((int)qwc.size() == m) qcw += qwc[j];
                  if (prev >= 0) {
                     auto &u = r3.cen[prev];
                     lc += std::sqrt((u.x - p.x) * (u.x - p.x) + (u.y - p.y) * (u.y - p.y) + (u.z - p.z) * (u.z - p.z));
                  }
                  prev = j;
               }
               pathC = lc;
               // FIRST ARC only (Spyral's choice): from the vertex end to the centroid FARTHEST from it
               // in x-y, so a stopping track's Bragg rise beyond the turning point is not included.
               {
                  // INLIERS of the vertex-end helix only: stray clusters attached to the track's end
                  // must not enter the charge, the path, or the "farthest point" reference.
                  auto inl = [&](const h3d::P &p) {
                     return std::abs(std::hypot(p.x - E.xc, p.y - E.yc) - E.R) < E.resCut;
                  };
                  int j0 = end3 == 0 ? 0 : m - 1, jFar = 0;
                  for (int jj = 0; jj < m; ++jj) {
                     int j = end3 == 0 ? jj : m - 1 - jj;
                     if (inl(r3.cen[j])) { j0 = j; break; }
                  }
                  double dFar = -1;
                  for (int jj = 0; jj < m; ++jj) {
                     int j = end3 == 0 ? jj : m - 1 - jj;
                     if (!inl(r3.cen[j])) continue;
                     double d = std::hypot(r3.cen[j].x - r3.cen[j0].x, r3.cen[j].y - r3.cen[j0].y);
                     if (d > dFar) { dFar = d; jFar = jj; }
                  }
                  double qf = 0, qfw = 0, lf = 0;
                  int pv = -1;
                  for (int jj = 0; jj <= jFar; ++jj) {
                     int j = end3 == 0 ? jj : m - 1 - jj;
                     auto &p = r3.cen[j];
                     if (!inl(p)) continue;
                     if (std::hypot(p.x, p.y) > 152.0) { pv = -1; continue; }
                     qf += p.q;
                     if ((int)qwc.size() == m) qfw += qwc[j];
                     if (pv >= 0) {
                        auto &u = r3.cen[pv];
                        lf += std::sqrt((u.x - p.x) * (u.x - p.x) + (u.y - p.y) * (u.y - p.y) + (u.z - p.z) * (u.z - p.z));
                     }
                     pv = j;
                  }
                  if (lf > 5) { sq3f = std::sqrt(qf / lf); if ((int)qwc.size() == m) sq3fw = std::sqrt(qfw / lf); }
                  pathF = lf;
                  // HELIX PATH of the same first-arc / small-pad segment: unwrap the centroid phases
                  // round the fitted vertex-end circle and take sqrt(R^2 + (dz/dphi)^2) * |dphi| from
                  // the first to the last inner centroid -- smooth, unlike the centroid polyline.
                  {
                     double phPrev = 0, phAcc = 0, phFirst = 0, phLast = 0;
                     bool started = false, stopped = false;
                     for (int jj = 0; jj <= jFar && !stopped; ++jj) {
                        int j = end3 == 0 ? jj : m - 1 - jj;
                        auto &p = r3.cen[j];
                        double a = std::atan2(p.y - E.yc, p.x - E.xc);
                        if (jj == 0) { phPrev = a; phAcc = 0; }
                        else {
                           double d = a - phPrev;
                           while (d > M_PI) d -= 2 * M_PI;
                           while (d < -M_PI) d += 2 * M_PI;
                           phAcc += d; phPrev = a;
                        }
                        bool inner = std::hypot(p.x, p.y) <= 152.0;
                        if (!inl(p)) continue;
                        if (inner && !started) { phFirst = phAcc; started = true; }
                        if (inner && started) phLast = phAcc;
                        if (!inner && started) stopped = true; // Spyral stops at the first exit too
                     }
                     if (started)
                        pathH = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi) * std::abs(phLast - phFirst);
                     if (pathH > 5) {
                        sq3h = std::sqrt(qf / pathH);
                        if ((int)qwc.size() == m) sq3hw = std::sqrt(qfw / pathH);
                        if (spyDE > 0) hybC = std::sqrt(spyDE / pathH);
                        if (spyDEw > 0) hybCw = std::sqrt(spyDEw / pathH);
                     }
                  }
                  // HYBRID B: Spyral's integrated (spline-smoothed) charge over ITS first arc, divided by
                  // the 3D first-arc path (Spyral's own arc length comes from x(z),y(z) splines, which
                  // collapse near 90 deg). Assumes both first arcs cover the same stretch of track.
                  if (lf > 5 && spyDE > 0) hybB = std::sqrt(spyDE / lf);
                  if (lf > 5 && spyDEw > 0) hybBw = std::sqrt(spyDEw / lf);
                  // TRUNCATED MEAN: local dE/dx per centreline step over the first arc (small pads),
                  // mean of the lowest 70 % -- the standard TPC answer to Landau fluctuations.
                  std::vector<double> loc, locw;
                  int pv2 = -1;
                  for (int jj = 0; jj <= jFar; ++jj) {
                     int j = end3 == 0 ? jj : m - 1 - jj;
                     auto &p = r3.cen[j];
                     if (std::hypot(p.x, p.y) > 152.0) { pv2 = -1; continue; }
                     if (pv2 >= 0) {
                        auto &u = r3.cen[pv2];
                        double ds = std::sqrt((u.x - p.x) * (u.x - p.x) + (u.y - p.y) * (u.y - p.y) + (u.z - p.z) * (u.z - p.z));
                        if (ds > 0.5) {
                           loc.push_back(p.q / ds);
                           if ((int)qwc.size() == m) locw.push_back(qwc[j] / ds);
                        }
                     }
                     pv2 = j;
                  }
                  auto tmean = [](std::vector<double> v) {
                     if (v.size() < 4) return -1.0;
                     std::sort(v.begin(), v.end());
                     size_t k = std::max<size_t>(2, (size_t)(0.7 * v.size()));
                     double sum = 0;
                     for (size_t i = 0; i < k; ++i) sum += v[i];
                     return sum / k;
                  };
                  double t1 = tmean(loc), t2 = tmean(locw);
                  if (t1 > 0) sq3t = std::sqrt(t1);
                  if (t2 > 0) sq3tw = std::sqrt(t2);
               }
               if (lc > 5) { sq3c = std::sqrt(qc / lc); if ((int)qwc.size() == m) sq3cw = std::sqrt(qcw / lc); }
            }
            if (pathIn > 5) { sq3a = std::sqrt(qa / pathIn); sq3w = std::sqrt(qww / pathIn); }
            else v3 = 0;
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
