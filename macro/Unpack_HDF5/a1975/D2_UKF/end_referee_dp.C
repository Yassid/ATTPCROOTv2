/// @file end_referee_dp.C
/// @brief END REFEREE for a1975 D2 data: measure every track from BOTH ends and keep the end whose
///        result is physically consistent.
///
/// For each end hypothesis k (helix3d.h, 3D, no z sort): Brho and theta_lab from the end-k helix,
/// sqrt(dE/dx) = charge / HELIX path over the first arc on the small pads (the estimator that won
/// the deuteron test, 2026-10-30), species from that end's own ln m band, then Ex for the matching
/// reaction on 16C + d:  p -> (d,p)17C,  d -> (d,d)16C,  t -> (d,t)15C.
/// An end is CONSISTENT if its species is assigned and |Ex| < exWin. Starting from the wrong end
/// measures the stopping part of the track (small radius, reversed direction, Bragg dE/dx), which
/// should almost never land on a reaction locus. Where EXACTLY ONE end is consistent, it is the
/// vertex end -- the referee -- and the two production-style calls are graded against it:
///   sense   3D rotation sense (+1 from the vertex end, measured in sim and H2 data)
///   dca     end whose helix passes closer to the beam axis
///   spyral  AtSpyralPID direction (forward = vertex at LOWEST z_digi)
///
///   root -l -b -q 'end_referee_dp.C("/mnt/f/a1975/reco_d2_tau/run_0016_multifit_reco.root")'
#include "../UKF/pid/helix3d.h"
#include <functional>

namespace eref {
const double u = 931.49410242, me = 0.51099895;
const double M16 = 16.014701 * u - 6 * me, Md = 2.01410178 * u - me;
const double Mp = 1.00782503 * u - me, Mt = 3.01604928 * u - me;
const double M17 = 17.022586 * u - 6 * me, M15 = 15.0105993 * u - 6 * me;
// Ex of the residual for ejectile mass m3 at (theta, Brho), beam 16C on d at Eb
double ex(double Eb, double m3, double M4, double th, double brho)
{
   double p = 299.792458 * brho, E3 = std::sqrt(p * p + m3 * m3);
   double Et = Eb + M16, pb = std::sqrt(Et * Et - M16 * M16), E4 = Et + Md - E3;
   double px = -p * std::sin(th), pz = pb - p * std::cos(th), m2 = E4 * E4 - px * px - pz * pz;
   return m2 > 0 ? std::sqrt(m2) - M4 : -999;
}
struct EndMeas {
   bool ok{false};
   double brho{0}, th{0}, sq{0}, dca{1e9};
};
// first-arc, small-pad charge over the helix path, measured from end `end`
EndMeas measure(const h3d::Result &r, int end)
{
   EndMeas M;
   auto &E = r.e[end];
   if (!E.ok) return M;
   const int m = r.cen.size();
   int j0 = end == 0 ? 0 : m - 1, jFar = 0;
   double dFar = -1;
   for (int jj = 0; jj < m; ++jj) {
      int j = end == 0 ? jj : m - 1 - jj;
      double d = std::hypot(r.cen[j].x - r.cen[j0].x, r.cen[j].y - r.cen[j0].y);
      if (d > dFar) { dFar = d; jFar = jj; }
   }
   double q = 0, phPrev = 0, phAcc = 0, phFirst = 0, phLast = 0;
   bool started = false;
   for (int jj = 0; jj <= jFar; ++jj) {
      int j = end == 0 ? jj : m - 1 - jj;
      auto &p = r.cen[j];
      double a = std::atan2(p.y - E.yc, p.x - E.xc);
      if (jj == 0) phPrev = a;
      else {
         double d = a - phPrev;
         while (d > M_PI) d -= 2 * M_PI;
         while (d < -M_PI) d += 2 * M_PI;
         phAcc += d; phPrev = a;
      }
      bool inner = std::hypot(p.x, p.y) <= 152.0;
      if (!inner && started) break;
      if (inner) {
         if (!started) { phFirst = phAcc; started = true; }
         phLast = phAcc;
         q += p.q;
      }
   }
   double path = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi) * std::abs(phLast - phFirst);
   if (!started || path < 5 || E.brho <= 0) return M;
   M.ok = true; M.brho = E.brho; M.th = E.theta; M.sq = std::sqrt(q / path); M.dca = E.dca;
   return M;
}
} // namespace eref

// minBrhoRatio: the referee only rules where the two ends measure DIFFERENT momenta (energy loss
// visible). With equal Brho the two hypotheses are mirror images (theta <-> 180-theta) and the locus
// test just picks the image closer to a low-Ex curve -- measured 2026-10-30: 41 of 56 "overrules"
// were that, near 90 deg, Brho ratio 1.00. 1.0 reproduces the unrestricted referee.
// BRAGG vote: charge per centreline bin (4 mm of track) on the small pads only (outer pads collect a
// different charge per hit); mean of the last third over the first third. Rising toward one end =
// that end is where the particle STOPS, so the vertex is the other one. Abstains within +-braggThr.
// simFile: if given, the MC truth vertex (first MCTrack of pdg truthPdg) referees the end instead.
void end_referee_dp(TString in, Long64_t maxEv = 1000000000, double Eb = 184.25, double exWin = 1.5,
                    double minBrhoRatio = 1.0, double braggThr = 0.2, TString simFile = "", int truthPdg = 2212)
{
   using namespace eref;
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

   TTree *ts = nullptr;
   TClonesArray *mc = nullptr;
   if (simFile.Length()) {
      auto *fs = TFile::Open(simFile);
      ts = (TTree *)fs->Get("cbmsim");
      ts->SetBranchStatus("*", 0);
      ts->SetBranchStatus("MCTrack*", 1);
      ts->SetBranchAddress("MCTrack", &mc);
   }
   struct Row {
      EndMeas m[2];
      int sense, dca, spy, bragg, truth, maj, majN;
      bool majAmb;
      double thHint, bRatio, rhoMax;
   };
   std::vector<Row> rows;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      double tx = 0, ty = 0, tz = -1e9;
      if (ts) {
         ts->GetEntry(ev);
         for (int j = 0; j < mc->GetEntriesFast(); ++j) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (m && m->GetPdgCode() == truthPdg) { tx = m->GetStartX() * 10; ty = m->GetStartY() * 10; tz = m->GetStartZ() * 10; break; }
         }
         if (tz < -1e8) continue;
      }
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         if (H.size() < 30) continue;
         std::vector<h3d::P> h;
         for (auto &x : H) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto r = h3d::estimate(h, 2.85, -1, 1000); // lab frame: z_lab = 1000 - z_digi
         if (!r.ok || r.length / r.nBackbone >= 6) continue; // coherent tracks only
         Row R;
         R.m[0] = measure(r, 0); R.m[1] = measure(r, 1);
         R.sense = r.e[0].sense == r.e[1].sense ? -1 : (r.e[0].sense == 1 ? 0 : 1);
         R.dca = std::abs(r.e[0].dca - r.e[1].dca) < 5 ? -1 : (r.e[0].dca < r.e[1].dca ? 0 : 1);
         auto sr = spy.Estimate(tr);
         // Spyral forward = vertex at LOWEST z_digi = HIGHEST z_lab. h is in the digi frame.
         double zA = h[r.order.front()].z, zB = h[r.order.back()].z;
         R.spy = !sr.valid ? -1 : ((sr.direction == 0) == (zA < zB) ? 0 : 1);
         R.thHint = r.e[0].theta;
         // Bragg profile on the small pads, along the centreline A -> B
         {
            std::vector<double> qi;
            R.rhoMax = 0;
            for (auto &c : r.cen) {
               R.rhoMax = std::max(R.rhoMax, std::hypot(c.x, c.y));
               if (std::hypot(c.x, c.y) <= 152.0) qi.push_back(c.q);
            }
            R.bragg = -1; R.bRatio = 0;
            if (qi.size() >= 6) {
               size_t k = qi.size() / 3;
               double a = 0, b = 0;
               for (size_t i = 0; i < k; ++i) { a += qi[i]; b += qi[qi.size() - 1 - i]; }
               R.bRatio = a > 0 ? b / a : 0;
               if (R.bRatio > 1 + braggThr) R.bragg = 0;            // rises toward B: stops at B, vertex A
               else if (R.bRatio < 1 / (1 + braggThr)) R.bragg = 1; // rises toward A
            }
         }
         {
            auto V = h3d::majorityEnd(r, R.spy, 1, braggThr);
            R.maj = V.end; R.majN = V.nFor; R.majAmb = V.ambiguous;
         }
         R.truth = -1;
         if (ts) { // the centroid end nearest the generated vertex (lab frame; z_lab = 1000 - z_digi measured)
            auto &a = r.cen.front(), &b = r.cen.back();
            double da = std::sqrt(std::pow(a.x - tx, 2) + std::pow(a.y - ty, 2) + std::pow(a.z - tz, 2));
            double db = std::sqrt(std::pow(b.x - tx, 2) + std::pow(b.y - ty, 2) + std::pow(b.z - tz, 2));
            if (std::min(da, db) < 40 && std::max(da, db) > 60) R.truth = da < db ? 0 : 1;
         }
         rows.push_back(R);
      }
   }
   // normalise sqrt(dE/dx) by its median over all end measurements, as in the other macros
   std::vector<double> all;
   for (auto &R : rows)
      for (int k = 0; k < 2; ++k) if (R.m[k].ok) all.push_back(R.m[k].sq);
   std::sort(all.begin(), all.end());
   double med = all[all.size() / 2];
   // species bands in ln(sqrt(dE/dx)/med * Brho), from the deuteron referee (d at -1.47)
   auto species = [&](const EndMeas &M) {
      double l = std::log(M.sq / med * M.brho);
      if (l > -2.35 && l < -1.75) return 1;  // p
      if (l > -1.70 && l < -1.25) return 2;  // d
      if (l > -1.20 && l < -0.80) return 3;  // t
      return 0;
   };
   auto consistent = [&](const EndMeas &M, int &sp, double &x) {
      if (!M.ok) return false;
      sp = species(M);
      if (sp == 1) x = ex(Eb, Mp, M17, M.th, M.brho);
      else if (sp == 2) x = ex(Eb, Md, M16, M.th, M.brho);
      else if (sp == 3) x = ex(Eb, Mt, M15, M.th, M.brho);
      else return false;
      return std::abs(x) < exWin;
   };

   long nRef = 0, nBoth = 0, nNone = 0, byS[4] = {0};
   long good[3] = {0}, dec[3] = {0}, goodDis[3] = {0}, decDis[3] = {0}, nDis = 0;
   auto decisive = [&](const Row &R) {
      if (!R.m[0].ok || !R.m[1].ok) return false;
      double q = R.m[0].brho / R.m[1].brho;
      return std::max(q, 1 / q) >= minBrhoRatio;
   };
   long nIndecisive = 0;
   for (auto &R : rows) {
      if (!decisive(R)) { ++nIndecisive; continue; }
      int s0, s1;
      double x0, x1;
      bool c0 = consistent(R.m[0], s0, x0), c1 = consistent(R.m[1], s1, x1);
      if (c0 && c1) { ++nBoth; continue; }
      if (!c0 && !c1) { ++nNone; continue; }
      int v = c0 ? 0 : 1;
      ++nRef;
      ++byS[c0 ? s0 : s1];
      int call[3] = {R.sense, R.dca, R.spy};
      bool dis = R.sense >= 0 && R.spy >= 0 && R.sense != R.spy;
      nDis += dis;
      for (int k = 0; k < 3; ++k) {
         if (call[k] < 0) continue;
         ++dec[k]; good[k] += call[k] == v;
         if (dis) { ++decDis[k]; goodDis[k] += call[k] == v; }
      }
   }
   printf("ends with Brho ratio < %.2f (mirror images, not refereed): %ld\n", minBrhoRatio, nIndecisive);
   printf("coherent tracks %zu: referee decides %ld (p %ld, d %ld, t %ld), both ends consistent %ld, neither %ld\n",
          rows.size(), nRef, byS[1], byS[2], byS[3], nBoth, nNone);
   const char *nm[3] = {"3D rotation sense", "3D closest approach", "Spyral direction"};
   printf("\n%-22s %18s %26s\n", "call", "all refereed", "where sense != Spyral");
   for (int k = 0; k < 3; ++k)
      printf("%-22s %6.1f %% (n=%5ld)   %6.1f %% (n=%4ld)\n", nm[k], dec[k] ? 100.0 * good[k] / dec[k] : 0, dec[k],
             decDis[k] ? 100.0 * goodDis[k] / decDis[k] : 0, decDis[k]);
   printf("(sense and Spyral disagree on %ld refereed tracks)\n", nDis);

   // ---- BRAGG ----
   auto grade = [&](const char *lab, std::function<int(const Row &)> ref, std::function<bool(const Row &)> use) {
      long n = 0, dB = 0, gB = 0, dS = 0, gS = 0, dP = 0, gP = 0;
      for (auto &R : rows) {
         if (!use(R)) continue;
         int v = ref(R);
         if (v < 0) continue;
         ++n;
         if (R.bragg >= 0) { ++dB; gB += R.bragg == v; }
         if (R.sense >= 0) { ++dS; gS += R.sense == v; }
         if (R.spy >= 0) { ++dP; gP += R.spy == v; }
      }
      long dM = 0, gM = 0, amb = 0;
      for (auto &R : rows) {
         if (!use(R)) continue;
         int v = ref(R);
         if (v < 0) continue;
         ++dM; gM += R.maj == v; amb += R.majAmb;
      }
      long dU = 0, gU = 0;
      for (auto &R : rows) {
         if (!use(R) || R.majAmb) continue;
         int v = ref(R);
         if (v < 0) continue;
         ++dU; gU += R.maj == v;
      }
      printf("  %-44s MAJORITY %5.1f %% (ambiguous %4.1f %%; on unambiguous %5.1f %%)\n", "", dM ? 100.0 * gM / dM : 0,
             dM ? 100.0 * amb / dM : 0, dU ? 100.0 * gU / dU : 0);
      printf("  %-44s n=%5ld | Bragg %5.1f %% (decides %3.0f %%) | sense %5.1f %% (%3.0f %%) | Spyral %5.1f %% (%3.0f %%)\n", lab,
             n, dB ? 100.0 * gB / dB : 0, n ? 100.0 * dB / n : 0, dS ? 100.0 * gS / dS : 0, n ? 100.0 * dS / n : 0,
             dP ? 100.0 * gP / dP : 0, n ? 100.0 * dP / n : 0);
   };
   auto kinRef = [&](const Row &R) {
      int s0, s1;
      double x0, x1;
      bool c0 = consistent(R.m[0], s0, x0), c1 = consistent(R.m[1], s1, x1);
      return c0 == c1 ? -1 : (c0 ? 0 : 1);
   };
   auto brRatio = [](const Row &R) {
      if (!R.m[0].ok || !R.m[1].ok) return 1.0;
      double q = R.m[0].brho / R.m[1].brho;
      return std::max(q, 1 / q);
   };
   printf("\nBRAGG vote (threshold +-%.0f %%):\n", 100 * braggThr);
   if (ts) {
      grade("SIM TRUTH, all", [](const Row &R) { return R.truth; }, [](const Row &) { return true; });
      grade("SIM TRUTH, short (end Brho ratio < 1.15)", [](const Row &R) { return R.truth; },
            [&](const Row &R) { return brRatio(R) < 1.15; });
   }
   grade("KINEMATIC referee (end Brho ratio >= 1.15)", kinRef, [&](const Row &R) { return brRatio(R) >= 1.15; });
   // short tracks: no referee; report how often Bragg agrees with the sense+Spyral consensus
   auto cons = [](const Row &R) { return (R.sense >= 0 && R.sense == R.spy) ? R.sense : -1; };
   grade("SHORT, vs sense+Spyral CONSENSUS (no referee)", cons, [&](const Row &R) { return brRatio(R) < 1.15; });
   {
      long n = 0, amb = 0, n3 = 0, n2 = 0, n1 = 0;
      for (auto &R : rows) { ++n; amb += R.majAmb; if (!R.majAmb) { n3 += R.majN == 3; n2 += R.majN == 2; n1 += R.majN == 1; } }
      printf("\nMAJORITY over all %ld coherent tracks: 3/3 %.1f %%, 2/3 %.1f %%, single voter %.1f %%, AMBIGUOUS %.1f %%\n", n,
             100.0 * n3 / n, 100.0 * n2 / n, 100.0 * n1 / n, 100.0 * amb / n);
   }
   // CONTAINED = the whole track stays on the small pads, so the profile sees where it stops
   grade("  short AND contained (rho_max < 150 mm)", cons,
         [&](const Row &R) { return brRatio(R) < 1.15 && R.rhoMax < 150; });
   grade("  short, NOT contained", cons, [&](const Row &R) { return brRatio(R) < 1.15 && R.rhoMax >= 150; });
   if (ts) {
      grade("SIM TRUTH, contained", [](const Row &R) { return R.truth; }, [&](const Row &R) { return R.rhoMax < 150; });
      grade("SIM TRUTH, not contained", [](const Row &R) { return R.truth; }, [&](const Row &R) { return R.rhoMax >= 150; });
   }
   grade("KINEMATIC referee, contained", kinRef, [&](const Row &R) { return brRatio(R) >= 1.15 && R.rhoMax < 150; });

   // where the referee OVERRULES sense and Spyral together: what do the two end hypotheses look like?
   printf("\nreferee overrules sense AND Spyral (they agree with each other):\n");
   printf("  %5s %6s | %-28s | %-28s\n", "", "", "end THEY pick (rejected)", "end REFEREE picks");
   printf("  %5s %6s | %6s %6s %6s %6s | %6s %6s %6s %6s\n", "sp", "th", "Brho", "lnm", "Ex", "dca", "Brho", "lnm", "Ex", "dca");
   int shown = 0;
   long nOver = 0, nOverNear90 = 0;
   std::vector<double> dBr;
   for (auto &R : rows) {
      int s0, s1;
      double x0 = -999, x1 = -999;
      bool c0 = consistent(R.m[0], s0, x0), c1 = consistent(R.m[1], s1, x1);
      if (!decisive(R) || c0 == c1 || R.sense < 0 || R.spy < 0 || R.sense != R.spy) continue;
      int v = c0 ? 0 : 1;
      if (R.sense == v) continue;
      ++nOver;
      int o = 1 - v;
      double thv = R.m[v].th * TMath::RadToDeg();
      if (std::abs(thv - 90) < 15) ++nOverNear90;
      dBr.push_back(R.m[o].brho / R.m[v].brho);
      // Ex of the rejected end under the referee end's species, for comparison
      int sv = species(R.m[v]);
      double xo = sv == 1 ? ex(Eb, Mp, M17, R.m[o].th, R.m[o].brho) : sv == 2 ? ex(Eb, Md, M16, R.m[o].th, R.m[o].brho)
                                                                             : ex(Eb, Mt, M15, R.m[o].th, R.m[o].brho);
      if (shown++ < 15)
         printf("  %5d %6.1f | %6.3f %6.2f %6.1f %6.1f | %6.3f %6.2f %6.1f %6.1f\n", sv, thv, R.m[o].brho,
                std::log(R.m[o].sq / med * R.m[o].brho), xo, R.m[o].dca, R.m[v].brho,
                std::log(R.m[v].sq / med * R.m[v].brho), c0 ? x0 : x1, R.m[v].dca);
   }
   std::sort(dBr.begin(), dBr.end());
   int nb = dBr.size();
   if (nb > 4)
      printf("  %ld overrules, %ld within 15 deg of 90; Brho(their end)/Brho(referee end) p16 %.2f med %.2f p84 %.2f\n",
             nOver, nOverNear90, dBr[nb * 16 / 100], dBr[nb / 2], dBr[nb * 84 / 100]);
}
