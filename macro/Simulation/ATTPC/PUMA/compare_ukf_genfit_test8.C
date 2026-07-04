/// @file compare_ukf_genfit_test8.C
/// @brief Compare UKF vs GENFIT on the PUMA branch-8 sample (back-to-back pi+/pi-,
///        |p| = 374.9 MeV/c, theta = 90 deg, vertex (0,0,75) mm).
///
/// Reads ./data/output_digi_both8.root (AtTrackingEventUKF + AtTrackingEventGenfit
/// + AtPatternEvent) and ./data/attpcsim.root (MCTrack + AtTpcPoint, correlated by
/// event index) and reports, per fitter:
///   - momentum resolution   (p_fit - 374.9)/374.9
///   - polar angle           theta_fit vs 90 deg
///   - charge-sign accuracy  vs truth (matched via (x,y) hits)
///   - vertex residual       (CAVEAT: UKF = back-extrap POCA; genfit = first-point state)
///
/// Run: root -b -q compare_ukf_genfit_test8.C
double medianOf(std::vector<double> v)
{
   if (v.empty()) return 0;
   std::sort(v.begin(), v.end());
   size_t n = v.size();
   return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
double iqrSigma(std::vector<double> v)
{
   if (v.size() < 4) return 0;
   std::sort(v.begin(), v.end());
   double q1 = v[v.size() / 4], q3 = v[3 * v.size() / 4];
   return (q3 - q1) / 1.349; // Gaussian-equivalent sigma
}

struct Stats {
   std::vector<double> dpFrac;  // (p-p0)/p0
   std::vector<double> theta;   // deg
   std::vector<double> dz, dr;  // vertex residual mm
   int nFit = 0, chargeRight = 0, chargeTot = 0;
};

void reportStats(const char *name, Stats &s, double p0)
{
   printf("\n\033[1;36m=== %s ===\033[0m\n", name);
   printf("  fits                : %d\n", s.nFit);
   if (s.nFit == 0) return;
   printf("  p resolution        : median %+.2f%%  sigma_IQR %.2f%%   (truth |p|=%.1f MeV/c)\n",
          100 * medianOf(s.dpFrac), 100 * iqrSigma(s.dpFrac), p0);
   printf("  theta [deg]         : median %.2f   sigma_IQR %.2f      (truth 90)\n", medianOf(s.theta),
          iqrSigma(s.theta));
   printf("  charge-sign accuracy: %.1f%%  (%d/%d)\n", s.chargeTot ? 100.0 * s.chargeRight / s.chargeTot : 0.,
          s.chargeRight, s.chargeTot);
   printf("  vertex dz [mm]      : median %+.2f  sigma_IQR %.2f\n", medianOf(s.dz), iqrSigma(s.dz));
   printf("  vertex dr(xy) [mm]  : median %.2f   sigma_IQR %.2f\n", medianOf(s.dr), iqrSigma(s.dr));
}

void compare_ukf_genfit_test8(TString digiFile = "./data/output_digi_both8.root",
                              TString simFile = "./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   const double m_pi = 139.57039;                                     // MeV
   const double p0 = std::sqrt(0.4 * 0.4 - (m_pi / 1000) * (m_pi / 1000)) * 1000; // 374.9 MeV/c
   const double vx0 = 0, vy0 = 0, vz0 = 75.0;                         // mm truth vertex
   const double kMatchTol = 3.0;                                      // mm (x,y) hit->MC match

   TFile fD(digiFile);
   TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);
   TTree *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("\033[1;31mmissing tree\033[0m\n"); return; }

   TClonesArray *ukfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *gfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   tD->SetBranchAddress("AtPatternEvent", &patArr);

   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   Stats sU, sG;
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());

   // Analyse one fitter's fits for the current event (arrays already loaded).
   auto analyse = [&](TClonesArray *teArr, Stats &s, std::vector<AtTrack> &trackArr, int nMC, std::vector<double> &mcX,
                      std::vector<double> &mcY, std::vector<int> &mcPdg) {
      if (teArr->GetEntries() == 0) return;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      const auto &fits = te->GetFittedTracks();
      for (const auto &ft : fits) {
         const auto &kin = ft->GetKinematics(0);
         double KE = kin.kineticEnergy;
         if (!(KE > 0)) continue;
         double p = std::sqrt(KE * KE + 2 * KE * m_pi);
         s.nFit++;
         s.dpFrac.push_back((p - p0) / p0);
         s.theta.push_back(kin.theta * 180.0 / M_PI);

         // vertex: UKF back-extrap (initialPositionXtr) if set, else genfit fVertex
         const auto &props = ft->GetTrackPropertiesStruct();
         double vx = props.initialPositionXtr.X(), vy = props.initialPositionXtr.Y(), vz = props.initialPositionXtr.Z();
         if (std::abs(vx) < 1e-9 && std::abs(vy) < 1e-9 && std::abs(vz) < 1e-9) {
            const auto &vtx = ft->GetVertex(0);
            vx = vtx.X(); vy = vtx.Y(); vz = vtx.Z();
         }
         s.dz.push_back(vz - vz0);
         s.dr.push_back(std::hypot(vx - vx0, vy - vy0));

         // fitted charge sign: genfit -> signed ParticleInfo.charge; UKF -> idPDG name "pi+"/"pi-"
         const auto &pinfo = ft->GetParticleInfo(0);
         int fitSign = 0;
         if (pinfo.charge != 0)
            fitSign = (pinfo.charge > 0) ? +1 : -1;
         else {
            TString id = pinfo.idPDG;
            if (id.Contains("+")) fitSign = +1;
            else if (id.Contains("-")) fitSign = -1;
         }

         // truth charge via (x,y) hit majority vote to MC pion PDG
         int tid = ft->GetTrackID();
         if (fitSign != 0 && tid >= 0 && tid < (int)trackArr.size()) {
            std::map<int, int> votes;
            for (const auto &h : trackArr[tid].GetHitArray()) {
               const auto &pp = h->GetPosition();
               double best = kMatchTol * kMatchTol;
               int bestPdg = 0;
               for (int k = 0; k < nMC; ++k) {
                  double d2 = (pp.X() - mcX[k]) * (pp.X() - mcX[k]) + (pp.Y() - mcY[k]) * (pp.Y() - mcY[k]);
                  if (d2 < best) { best = d2; bestPdg = mcPdg[k]; }
               }
               if (bestPdg != 0) votes[bestPdg]++;
            }
            int truthPdg = 0, bestV = 0;
            for (auto &kv : votes) if (kv.second > bestV) { bestV = kv.second; truthPdg = kv.first; }
            if (truthPdg != 0) {
               int truthSign = (truthPdg > 0) ? +1 : -1; // +211 -> +, -211 -> -
               s.chargeTot++;
               if (truthSign == fitSign) s.chargeRight++;
            }
         }
      }
   };

   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e);
      tS->GetEntry(e);
      if (patArr->GetEntries() == 0) continue;
      auto *pat = (AtPatternEvent *)patArr->At(0);
      auto trackArr = pat->GetTrackCand(); // vector<AtTrack>

      int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC);
      std::vector<int> mcPdg(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10.0;
         mcY[k] = mp->GetY() * 10.0;
         int tid = mp->GetTrackID();
         auto *mt = (tid >= 0 && tid < mcTrks->GetEntries()) ? (AtMCTrack *)mcTrks->At(tid) : nullptr;
         mcPdg[k] = mt ? mt->GetPdgCode() : 0;
      }

      analyse(ukfArr, sU, trackArr, nMC, mcX, mcY, mcPdg);
      analyse(gfArr, sG, trackArr, nMC, mcX, mcY, mcPdg);
   }

   printf("\n\033[1;33m########  UKF vs GENFIT  —  PUMA branch 8 (pi+/pi-)  ########\033[0m\n");
   printf("events analysed: %lld   truth: |p|=%.1f MeV/c, theta=90 deg, vertex (0,0,75) mm\n", nE, p0);
   reportStats("UKF  (AtFitterUKF, CATIMA)", sU, p0);
   reportStats("GENFIT (KalmanRefTrack, Bethe-Bloch)", sG, p0);
   printf("\n\033[2mNote: vertex for UKF is the back-extrapolated POCA (initialPositionXtr);\n"
          "      for genfit it is the fitted first-point state — not yet extrapolated to the\n"
          "      beam axis, so genfit vertex dr is expected larger by construction.\033[0m\n");
}
