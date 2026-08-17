#include "AtMCQMinimization.h"

#include "AtDigiPar.h"
#include "AtELossModel.h"
#include "AtHit.h"
#include "AtMap.h"

#include <FairLogger.h>
#include <FairRun.h>
#include <FairRuntimeDb.h>

#include <Math/Point2D.h>
#include <TMath.h>
#include <TRandom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

/* Ported from the FairRootv18.00 branch (AtMCQMinimization, Y. Ayyad and W. Mittig). The
 * simulation, the objective function and the random walk are the ones of MinimizeGen(); the
 * following changes were made:
 *
 *  - The pad plane is taken from AtMap (GetPadNum() and CalcPadCenter()) instead of being passed
 *    as a TH2Poly plus a boost::multi_array of pad corners. The histogram is no longer filled, so
 *    the minimizer no longer has a side effect on the caller's pad plane.
 *  - The 8 element `Double_t parameter[]` array became AtMinimization::TrackSeed and the 19
 *    reference arguments of MCvar() became MCState, split in SeedState() and VaryState().
 *  - The tilt, Lorentz, rotation and pad angles and the maximum range are no longer in AtDigiPar.
 *    They default to an aligned detector and have setters.
 *  - The energy loss and range to energy `std::function` arrays (one per particle, only the first
 *    ever used) became a single function each. The recommended way to set them is SetELossModel(),
 *    which takes the stopping power from CATIMA (AtTools::AtELossCATIMA) and inverts its range;
 *    the parametrized functions of the old macros are still accepted by the setters.
 *  - The state of a fit is cleared at the start of every Minimize(). The original appended the
 *    experimental and simulated tracks of every event to the same vectors.
 *  - Hits whose pad number is out of range are skipped instead of writing out of bounds.
 *  - Pad 0 now takes part in the fit: the original kept `bin > 0` after converting a TH2Poly bin
 *    into a pad, which dropped it.
 *  - Time buckets that the simulated track never reached hold a sentinel instead of the content of
 *    the previous simulation, so their contribution to the position chi2 is the (clamped) maximum
 *    penalty instead of an undefined value.
 *  - Chi2Range() uses the last point of the simulated trajectory. The original used the last entry
 *    of a vector that also held the pad centroids, i.e. the pad with the highest number.
 *  - Dead code was dropped: the deprecated MinimizeOptMapAmp(), the four functions that returned
 *    nothing, the gain calibration plots, and the buffers that were filled but never read.
 *
 * The random walk draws the same number of variates in the same order as the original, so the two
 * implementations explore the same sequence of trial tracks.
 */

namespace MCMinimization {

namespace {
/// Charge threshold of the simulated track, in number of primary electrons.
constexpr double kSimThreshold = 5.;
/// Charge a pad needs, in both the simulation and the data, to enter the charge chi2 [ADC].
constexpr double kChi2Threshold = 50.;
/// Squared distance between a pad and the track above which the two are unrelated, in units of
/// the width of the charge weighted position.
constexpr double kMaxPullZ2 = 20.;
/// Number of pads below which the charge comparison is not meaningful, and its penalty.
constexpr int kMinChi2Pads = 5;
constexpr double kChi2NoPads = 10000.;
/// Error on the center of gravity of the hits of a time bucket [mm2].
constexpr double kPosSigma2 = 36.;
/// Largest contribution a single time bucket can make to the position chi2, once the steps of
/// the random walk have shrunk and while they are still large.
constexpr double kPosChi2Cap = 10.;
constexpr double kPosChi2CapCoarse = 100.;
/// Factor the steps of the random walk shrink by at every coarse iteration.
constexpr double kStepShrink = 1.4;
/// Atomic mass unit [MeV], as used by the original implementation.
constexpr double kAMU = 931.49432;
/// Speed of light [cm/ns].
constexpr double kSpeedOfLight = 29.9792;
} // namespace

bool AtMCQMinimization::Init()
{
   if (fMap == nullptr) {
      LOG(error) << "AtMCQMinimization: no pad plane map, call SetMap() before minimizing!";
      return false;
   }

   auto *run = FairRun::Instance();
   auto *db = run == nullptr ? nullptr : run->GetRuntimeDb();
   auto *par = db == nullptr ? nullptr : dynamic_cast<AtDigiPar *>(db->getContainer("AtDigiPar")); // NOLINT
   if (par == nullptr) {
      LOG(error) << "AtMCQMinimization: AtDigiPar not found, cannot initialize!";
      return false;
   }

   fDriftVelocity = par->GetDriftVelocity();
   fTBTime = par->GetTBTime();
   fBField = par->GetBField();
   fZPadPlane = par->GetZPadPlane();
   fDensity = par->GetDensity();
   fPressure = par->GetGasPressure();
   if (!fEntTBSet)
      fEntTB = par->GetTBEntrance();

   // Charge response, identical to the one of the digitization chain
   fEIonize = par->GetEIonize() / 1.e6;   // [MeV], as AtClusterize does
   fCoefT = par->GetCoefDiffusionTrans(); // [cm2/us]
   fCoefL = par->GetCoefDiffusionLong();  // [cm2/us]
   if (!fGainSet)
      fGain = GainFromParameters(*par); // ADC per primary electron, as AtPulse does

   if (fEIonize <= 0) {
      LOG(error) << "AtMCQMinimization: the ionization energy of the parameter file is " << fEIonize << " MeV!";
      return false;
   }

   fDzStep = fDriftVelocity * fTBTime / 1000.; // [cm]
   fBFieldGauss = fBField * 1.e4 * fBFieldSign;

   const auto numPads = fMap->GetNumPads();
   fQExp.assign(numPads, 0.);
   fZExp.assign(numPads, 0.);
   fQSim.assign(numPads, 0.);
   fZSim.assign(numPads, 0.);
   fTBTrack.assign(kMaxIntegrationSteps, XYZPoint(kNoPoint, kNoPoint, kNoPoint));

   fStepX.assign(kMaxIntegrationSteps, 0.);
   fStepY.assign(kMaxIntegrationSteps, 0.);
   fStepZ.assign(kMaxIntegrationSteps, 0.);
   fStepQ.assign(kMaxIntegrationSteps, 0.);

   LOG(info) << "AtMCQMinimization: drift velocity " << fDriftVelocity << " cm/us, time bucket " << fTBTime << " ns, B "
             << fBField * fBFieldSign << " T, density " << fDensity << ", pressure " << fPressure << " torr";
   LOG(info) << "AtMCQMinimization: ionization energy " << fEIonize * 1.e6 << " eV, transverse diffusion " << fCoefT
             << " cm2/us, " << fGain << " ADC per primary electron";

   fIsInit = true;
   return true;
}

double AtMCQMinimization::GainFromParameters(const AtDigiPar &par)
{
   // Same conversion AtPulse applies: the charge of one electron amplified by the micromegas,
   // in units of the full scale of the GET electronics, whose gain is given in fC.
   constexpr double kElectronCharge = 1.602e-19; // [C]
   constexpr double kADCFullScale = 4096.;

   const double getGain = par.GetGETGain(); // [fC]
   if (getGain <= 0) {
      LOG(error) << "AtMCQMinimization: the GET gain of the parameter file is " << getGain << " fC!";
      return 0.;
   }

   return par.GetGain() * kElectronCharge * kADCFullScale / (getGain * 1.e-15);
}

void AtMCQMinimization::SetParticle(int A, int Z)
{
   fA = A;
   fZ = Z;
   fRestMass = fA * kAMU;
   // Charge over mass in [e/m electron cm2/(V ns2)]
   fChargeToMass = fZ * 1.75879e-3 * 0.510998918 / fRestMass;
}

void AtMCQMinimization::SetELossModel(ELossModelPtr model)
{
   fELossModel = std::move(model);
   if (fELossModel == nullptr) {
      LOG(error) << "AtMCQMinimization: the energy loss model is a null pointer!";
      return;
   }

   // The model gives MeV/mm at the density it was built for; the simulation works in MeV/cm
   fStoppingPower = [this](double energy) { return fELossModel->GetdEdx(energy) * 10.; };
   fRangeToEnergy = [this](double range) { return RangeToEnergyFromModel(range); };
}

double AtMCQMinimization::RangeToEnergyFromModel(double range) const
{
   // The range grows monotonically with the energy, so a bisection is enough to invert it
   constexpr double kMinEnergy = 1.e-3; // [MeV]
   constexpr double kMaxEnergy = 1.e3;  // [MeV]
   constexpr double kTolerance = 1.e-4; // Relative width of the bracket
   constexpr int kMaxIter = 60;

   if (range <= 0.)
      return 0.;
   if (fELossModel->GetRange(kMaxEnergy) < range) {
      LOG(warn) << "AtMCQMinimization: a range of " << range << " mm is longer than the range of a " << kMaxEnergy
                << " MeV particle, the energy is clipped!";
      return kMaxEnergy;
   }

   double low = kMinEnergy;
   double high = kMaxEnergy;
   for (int i = 0; i < kMaxIter && (high - low) > kTolerance * high; i++) {
      const double mid = 0.5 * (low + high);
      if (fELossModel->GetRange(mid) < range)
         low = mid;
      else
         high = mid;
   }

   return 0.5 * (low + high);
}

void AtMCQMinimization::SetStoppingPower(ParametrisedFunc func, std::vector<double> par)
{
   // The parametrizations of the FairRootv18.00 macros give the mass stopping power in the units
   // of the density of the parameter file, i.e. MeV/(mg/cm2) for the AT-TPC parameter files
   fStoppingPower = [this, func = std::move(func), par = std::move(par)](double energy) mutable {
      return func(energy, par) * fDensity;
   };
}

void AtMCQMinimization::SetRangeToEnergy(ParametrisedFunc func, std::vector<double> par)
{
   // The parametrizations of the FairRootv18.00 macros are defined at 760 torr
   fRangeToEnergy = [this, func = std::move(func), par = std::move(par)](double range) mutable {
      return func(range * fPressure / 760., par);
   };
}

void AtMCQMinimization::SetStepParameters(const std::vector<double> &par)
{
   if (par.size() < fStepPar.size()) {
      LOG(error) << "AtMCQMinimization: expected at least " << fStepPar.size() << " step parameters, got " << par.size()
                 << ", keeping the current ones!";
      return;
   }
   std::copy_n(par.begin(), fStepPar.size(), fStepPar.begin());
}

void AtMCQMinimization::SetNumIterations(int coarse, int fine)
{
   fNumCoarseIter = coarse;
   fNumFineIter = fine;
}

void AtMCQMinimization::ResetParameters()
{
   AtMinimization::ResetParameters();
   ClearTrack();
}

void AtMCQMinimization::ClearTrack()
{
   std::fill(fQExp.begin(), fQExp.end(), 0.);
   std::fill(fZExp.begin(), fZExp.end(), 0.);
   std::fill(fQSim.begin(), fQSim.end(), 0.);
   std::fill(fZSim.begin(), fZSim.end(), 0.);
   std::fill(fTBTrack.begin(), fTBTrack.end(), XYZPoint(kNoPoint, kNoPoint, kNoPoint));

   fExpTrack.clear();
   fExpTB.clear();
   fHitsByTB.clear();
   fSimTrack.clear();
   fSimPads.clear();
   fSimCharge.clear();
   fBackTrack.clear();

   fHasExpEndPoint = false;
   fSimEnergy = 0;
   fBeamRange = 0;
   fChi2Points = 0;
   fLastTB = 0;
}

bool AtMCQMinimization::Minimize(const TrackSeed &seed, const HitVector &hits)
{
   if (!fIsInit && !Init())
      return false;

   if (fA == 0 || fZ == 0) {
      LOG(error) << "AtMCQMinimization: the particle was not set, call SetParticle()!";
      return false;
   }
   if (!fStoppingPower) {
      LOG(error) << "AtMCQMinimization: no stopping power function, call SetStoppingPower()!";
      return false;
   }
   if (fBFieldGauss == 0 && !fRangeToEnergy) {
      LOG(error) << "AtMCQMinimization: without a magnetic field the energy is taken from the range, "
                 << "call SetRangeToEnergy()!";
      return false;
   }
   if (hits.empty()) {
      LOG(error) << "AtMCQMinimization: nothing to fit, the hit array is empty!";
      return false;
   }

   ResetParameters();
   FillExperimentalTrack(hits, seed.fVertexTB);

   auto best = SeedState(seed);
   LOG(info) << "AtMCQMinimization: fitting A " << fA << " Z " << fZ << " seeded at (" << best.fX << ", " << best.fY
             << ", " << best.fZ << ") cm, theta " << best.fTheta * TMath::RadToDeg() << " deg, phi "
             << best.fPhi * TMath::RadToDeg() << " deg, radius " << best.fRadius << " mm, with " << seed.fNumExpPoints
             << " experimental points";

   double chi2Best = 1.e7;
   double chi2QBest = 0;
   double chi2RangeBest = 0;

   for (int iCoarse = 0; iCoarse < fNumCoarseIter; iCoarse++) {
      for (int iFine = 0; iFine < fNumFineIter; iFine++) {

         auto trial = VaryState(best, iCoarse);
         SimulateTrack(trial);

         auto chi2Q = Chi2Q();
         double chi2Pos = fUsePosChi2 ? Chi2Pos(iCoarse, seed.fNumExpPoints) : 0.;
         double chi2Range = fUseRangeChi2 ? Chi2Range() : 0.;
         double chi2 = fUseRangeChi2 ? (chi2Q + chi2Range) / 2. : (chi2Q + chi2Pos) / 2.;

         if (chi2 < chi2Best) {
            chi2Best = chi2;
            chi2QBest = chi2Q;
            chi2RangeBest = chi2Range;
            best = trial;
         }
      }
   }

   // Regenerate the track of the best set of parameters
   SimulateTrack(best);

   fFitPar.fTheta = best.fTheta;
   fFitPar.fPhi = best.fPhi;
   fFitPar.fEnergy = fSimEnergy;
   fFitPar.fPos = XYZPoint(best.fX, best.fY, best.fZ);
   fFitPar.fRadius = best.fRadius;
   fFitPar.fB = best.fB;
   fFitPar.fChi2 = chi2Best;
   fFitPar.fChi2Q = chi2QBest;
   fFitPar.fChi2Range = chi2RangeBest;
   fFitPar.fNumMCPoint = fSimTrack.size();

   if (fBackwardProp) {
      BackwardExtrapolation(best);
   } else {
      fFitPar.fVertexPos = fFitPar.fPos;
      if (fRangeToEnergy)
         fFitPar.fVertexEner = fRangeToEnergy(fMaxRange - fBeamRange);
   }

   if (fVerbose)
      PrintResult();

   return true;
}

void AtMCQMinimization::FillExperimentalTrack(const HitVector &hits, int vertexTB)
{
   const auto numPads = static_cast<int>(fQExp.size());

   for (const auto &hitPtr : hits) {
      const auto &hit = *hitPtr;
      const auto pad = hit.GetPadNum();
      if (pad < 0 || pad >= numPads) {
         LOG(debug) << "AtMCQMinimization: skipping hit " << hit.GetHitID() << " with pad number " << pad;
         continue;
      }

      // A pad can hold several hits. The simulated side accumulates every electron that lands on
      // a pad and takes the charge weighted position, so the data is summed the same way.
      fQExp[pad] += hit.GetCharge();
      fZExp[pad] += hit.GetCharge() * hit.GetPosition().Z(); // Normalized below

      fExpTrack.push_back(hit.GetPosition());
      fExpTB.push_back(hit.GetTimeStamp());

      // The simulated track is propagated towards the pad plane, so it ends at the hit of lowest z
      if (!fHasExpEndPoint || hit.GetPosition().Z() < fExpEndPoint.Z()) {
         fExpEndPoint = hit.GetPosition();
         fHasExpEndPoint = true;
      }
   }

   for (int pad = 0; pad < numPads; pad++)
      if (fQExp[pad] > 0)
         fZExp[pad] /= fQExp[pad];

   if (!fUsePosChi2)
      return;

   // The hits of the track sorted by descending time bucket, starting at the one of the vertex,
   // which is the order the simulated track is generated in
   fHitsByTB.assign(std::max(vertexTB, 0), {});
   for (const auto &hitPtr : hits) {
      const int index = vertexTB - hitPtr->GetTimeStamp();
      if (index >= 0 && index < vertexTB)
         fHitsByTB[index].push_back(hitPtr.get());
   }
}

AtMCQMinimization::MCState AtMCQMinimization::SeedState(const TrackSeed &seed) const
{
   MCState state;

   state.fX = seed.fVertex.X() / 10.;
   state.fY = seed.fVertex.Y() / 10.;

   if (fUseGeoVertex) {
      // Vertex found by the pattern recognition, relative to the calibrated entrance position
      state.fZ = seed.fVertex.Z() / 10. + (fZPadPlane - fEntZ0) / 10.;
      state.fPhi = seed.fPhi;
   } else {
      // Vertex given by its time bucket
      state.fZ = fZPadPlane / 10. - (fEntTB - seed.fVertexTB) * fDzStep;
      state.fPhi = TMath::Pi() - seed.fPhi - fThetaPad;
   }

   state.fTheta = seed.fTheta;
   state.fRadius = seed.fRadius;
   state.fB = fBFieldGauss;
   state.fDensityScale = 1.;

   return state;
}

AtMCQMinimization::MCState AtMCQMinimization::VaryState(const MCState &state, int iteration) const
{
   // The step is divided by 1.4 every time the coarse iteration number goes up
   const double factStep = std::pow(kStepShrink, -iteration);
   const auto step = [this, factStep](size_t i) { return fStepPar[i] * factStep * (0.5 - gRandom->Rndm()); };

   MCState varied = state;
   varied.fTheta = state.fTheta + step(0) * TMath::DegToRad();
   varied.fPhi = state.fPhi + step(1) * TMath::DegToRad();
   varied.fRadius = state.fRadius * (1. + step(2));
   varied.fX = state.fX + step(3);
   varied.fY = state.fY + step(4);
   varied.fZ = state.fZ + step(5);
   varied.fB = state.fB * (1. + step(6));
   varied.fDensityScale = state.fDensityScale * (1. + step(7));

   // The original implementation drew a ninth variate and never used it. It is kept so that both
   // implementations walk through the same sequence of trial tracks.
   gRandom->Rndm();

   return varied;
}

void AtMCQMinimization::SimulateTrack(const MCState &state)
{
   std::fill(fQSim.begin(), fQSim.end(), 0.);
   std::fill(fZSim.begin(), fZSim.end(), 0.);
   std::fill(fTBTrack.begin(), fTBTrack.end(), XYZPoint(kNoPoint, kNoPoint, kNoPoint));
   fSimTrack.clear();
   fSimPads.clear();
   fSimCharge.clear();

   /*** Energy of the particle entering the gas ***/
   double energy = 0;      // Per nucleon when the track is seeded from its curvature
   double totalEnergy = 0; // Kinetic energy followed along the track [MeV]

   if (fBFieldGauss != 0) {
      // Magnetic rigidity corrected for the angle [Tm], signed like the field
      const double bRhoTheta = state.fB * state.fRadius * 1.e-7 / TMath::Sin(state.fTheta);
      energy = GetEnergy(fA, fZ, bRhoTheta);
      totalEnergy = energy * fA;
   } else {
      // Without a magnetic field the radius of the seed holds the range of the track
      energy = fRangeToEnergy(state.fRadius);
      totalEnergy = energy;
   }
   fSimEnergy = energy;

   /*** Initial conditions, in the laboratory frame ***/
   const auto posLab = TransformIniPos(XYZPoint(state.fX, state.fY, state.fZ));
   double x = posLab.X();
   double y = posLab.Y();
   double z = posLab.Z();
   const double zIni = z;

   double ekin = totalEnergy;
   const double v0 = TMath::Sqrt(2. * totalEnergy / (fA * 931.49)) * kSpeedOfLight; // [cm/ns]
   double dt = fDzStep / (v0 * TMath::Cos(state.fTheta)) / fIntegrationSteps;       // [ns]

   double dxdt = v0 * TMath::Sin(state.fTheta) * TMath::Cos(state.fPhi);
   double dydt = v0 * TMath::Sin(state.fTheta) * TMath::Sin(state.fPhi);
   double dzdt = v0 * TMath::Cos(state.fTheta);

   /*** Propagation ***/
   int numSteps = 0;
   int tbIni = 0;

   for (int k = 0; k < kMaxIntegrationSteps; k++) {
      numSteps++;

      // Laboratory frame to pad plane frame, in mm. The direction of propagation is reversed
      // because the track is simulated from the pad plane towards the entrance window.
      const double xcmm = x * 10.;
      const double ycmm = y * 10.;
      const double zcmm = -z * 10. + 2. * zIni * 10.;

      const double xsol = xcmm - zcmm * TMath::Sin(fThetaLorentz) * TMath::Sin(fThetaRot);
      const double ysol = ycmm + zcmm * TMath::Sin(fThetaLorentz) * TMath::Cos(fThetaRot);
      const double zsol = zcmm;

      const double xdet = xsol;
      const double ydet = -(fZPadPlane - zsol) * TMath::Sin(fTiltAng) + ysol * TMath::Cos(fTiltAng);
      const double zdet = zsol * TMath::Cos(fTiltAng) - ysol * TMath::Sin(fTiltAng);

      const double xpad = xdet * TMath::Cos(fThetaPad) - ydet * TMath::Sin(fThetaPad);
      const double ypad = xdet * TMath::Sin(fThetaPad) + ydet * TMath::Cos(fThetaPad);
      double zpad = zdet;

      if (fUseGeoVertex)
         zpad -= fZPadPlane - fEntZ0;
      if (k == 0 && !fBackwardProp)
         fBeamRange = fEntZ0 - zpad;

      fStepX[k] = xpad;
      fStepY[k] = ypad;
      fStepZ[k] = zpad;
      fStepQ[k] = 0.;

      // Time bucket of this point, relative to the one the track started at
      const int iterCorr = static_cast<int>(zpad / (fDzStep * 10.) + 0.5);
      if (k == 0)
         tbIni = iterCorr;

      const int tb = tbIni - iterCorr;
      fLastTB = tb;
      if (tb < 0)
         break;

      if (tb < static_cast<int>(fTBTrack.size()))
         fTBTrack[tb] = XYZPoint(xpad, ypad, zpad);
      fSimTrack.emplace_back(xpad, ypad, zpad);

      /*** One step of the equation of motion ***/
      const double ddxddt = fChargeToMass * state.fB * 10. * dydt;
      const double ddyddt = -fChargeToMass * state.fB * 10. * dxdt;

      x += dxdt * dt + 0.5 * ddxddt * dt * dt;
      y += dydt * dt + 0.5 * ddyddt * dt * dt;
      z += dzdt * dt;

      dxdt += ddxddt * dt;
      dydt += ddyddt * dt;

      // NB: as in the original implementation, the length of the step is evaluated with the
      // velocity at the end of the step and not with the one that was used to move the particle.
      const double step = TMath::Sqrt(TMath::Power(dxdt * dt + 0.5 * ddxddt * dt * dt, 2) +
                                      TMath::Power(dydt * dt + 0.5 * ddyddt * dt * dt, 2) + TMath::Power(dzdt * dt, 2));

      /*** Energy loss ***/
      const double eLoss = fStoppingPower(ekin) * state.fDensityScale * step; // [MeV]
      fStepQ[k] = eLoss / fEIonize;                                           // Number of primary electrons

      const double beta = TMath::Sqrt(dxdt * dxdt + dydt * dydt + dzdt * dzdt) / kSpeedOfLight;
      const double ekinNoLoss = fA * 931.494 * 0.5 * beta * beta; // Non relativistic
      ekin -= eLoss;

      const double slowDown = TMath::Sqrt(ekin / ekinNoLoss);
      dxdt *= slowDown;
      dydt *= slowDown;
      dzdt *= slowDown;

      dt = fDzStep / dzdt / fIntegrationSteps;

      if (zpad < 0. || ekin < 0.01 || std::isnan(ekin))
         break;
   }

   /*** Spread the primary electrons over the pad plane ***/
   // Radii of the 4 rings holding a quarter of the transverse straggling each, in units of sigma
   constexpr std::array<double, 4> kStragglingRadii = {0.37925, 0.968, 1.421, 6.69};
   constexpr int kNumStragglingPoints = 8;
   const auto numPads = static_cast<int>(fQSim.size());

   for (int k = 0; k < numSteps; k++) {
      // Transverse diffusion of the electrons of this step, the same expression AtClusterize
      // uses. NB: the longitudinal one is not simulated, only the center of gravity is used.
      const double driftTime = std::max(0., fStepZ[k]) / (10. * fDriftVelocity); // [us]
      const double sigmaTrans = 10. * TMath::Sqrt(2. * fCoefT * driftTime);      // [mm]
      const double charge = fStepQ[k] / (kStragglingRadii.size() * kNumStragglingPoints);

      for (auto radius : kStragglingRadii) {
         for (int iPhi = 0; iPhi < kNumStragglingPoints; iPhi++) {
            const double phi = TMath::TwoPi() / kNumStragglingPoints * iPhi;
            const ROOT::Math::XYPoint point(fStepX[k] + TMath::Cos(phi) * radius * sigmaTrans,
                                            fStepY[k] + TMath::Sin(phi) * radius * sigmaTrans);

            const auto pad = fMap->GetPadNum(point);
            if (pad < 0 || pad >= numPads)
               continue;

            fQSim[pad] += charge;
            fZSim[pad] += charge * fStepZ[k]; // Normalized below
         }
      }
   }

   /*** Normalize the simulated track and bring its charge to the scale of the data ***/
   for (int pad = 0; pad < numPads; pad++) {
      const double electrons = fQSim[pad];
      if (electrons <= 0)
         continue;

      fZSim[pad] /= electrons; // Charge weighted position of the pad
      fQSim[pad] = electrons * fGain;

      if (electrons <= kSimThreshold)
         continue;

      const auto center = fMap->CalcPadCenter(pad);
      fSimPads.emplace_back(center.X(), center.Y(), fZSim[pad]);
      fSimCharge.push_back(electrons);
   }
}

void AtMCQMinimization::BackwardExtrapolation(const MCState &state)
{
   constexpr int kMaxSteps = 200;

   // NB: fSimEnergy is the energy per nucleon only when the track was seeded from its curvature.
   double ekin = fSimEnergy * fA;
   const double v0 = TMath::Sqrt(2. * ekin / (fA * 931.49)) * kSpeedOfLight; // [cm/ns]
   double dt = -fDzStep / (v0 * TMath::Cos(state.fTheta)) / fIntegrationSteps;

   const auto posLab = TransformIniPos(XYZPoint(state.fX, state.fY, state.fZ));
   double x = posLab.X();
   double y = posLab.Y();
   double z = posLab.Z();
   const double zIni = z;

   double dxdt = v0 * TMath::Sin(state.fTheta) * TMath::Cos(state.fPhi);
   double dydt = v0 * TMath::Sin(state.fTheta) * TMath::Sin(state.fPhi);
   double dzdt = v0 * TMath::Cos(state.fTheta);

   double minDist = 1.e10;
   int tbIni = 0;
   bool foundVertex = false;

   for (int k = 0; k < kMaxSteps; k++) {

      // Distance to the beam axis, which is assumed to go through the center of the pad plane
      const double dist = TMath::Sqrt(x * x + y * y);
      if (dist < minDist) {
         minDist = dist;
      } else {
         // The track is going away from the beam axis: the previous point was the vertex
         fFitPar.fVertexPos = XYZPoint(x, y, 2. * zIni - z);
         fFitPar.fVertexEner = ekin;
         fFitPar.fMinDistAppr = dist;
         foundVertex = true;
         break;
      }

      // Laboratory frame to pad plane frame, in mm
      const double xcmm = x * 10.;
      const double ycmm = y * 10.;
      const double zcmm = -z * 10. + 2. * zIni * 10.;

      const double xsol = xcmm - zcmm * TMath::Sin(fThetaLorentz) * TMath::Sin(fThetaRot);
      const double ysol = ycmm + zcmm * TMath::Sin(fThetaLorentz) * TMath::Cos(fThetaRot);
      const double zsol = zcmm;

      const double xdet = xsol;
      const double ydet = -(fZPadPlane - zsol) * TMath::Sin(fTiltAng) + ysol * TMath::Cos(fTiltAng);
      const double zdet = zsol * TMath::Cos(fTiltAng) - ysol * TMath::Sin(fTiltAng);

      const double xpad = xdet * TMath::Cos(fThetaPad) - ydet * TMath::Sin(fThetaPad);
      const double ypad = xdet * TMath::Sin(fThetaPad) + ydet * TMath::Cos(fThetaPad);
      const double zpad = zdet;

      const int iterCorr = static_cast<int>(zpad / (fDzStep * 10.) + 0.5);
      if (k == 0)
         tbIni = iterCorr;
      if (iterCorr - tbIni >= 0)
         fBackTrack.emplace_back(xpad, ypad, zpad);

      /*** One step of the equation of motion, backwards ***/
      // NB: as in the original implementation the field is the one of the parameter file and not
      // the (by default not varied) field of the fit.
      const double ddxddt = fChargeToMass * fBFieldGauss * 10. * dydt;
      const double ddyddt = -fChargeToMass * fBFieldGauss * 10. * dxdt;

      x += dxdt * dt + 0.5 * ddxddt * dt * dt;
      y += dydt * dt + 0.5 * ddyddt * dt * dt;
      z += dzdt * dt;

      dxdt += ddxddt * dt;
      dydt += ddyddt * dt;

      const double step = TMath::Sqrt(TMath::Power(dxdt * dt + 0.5 * ddxddt * dt * dt, 2) +
                                      TMath::Power(dydt * dt + 0.5 * ddyddt * dt * dt, 2) + TMath::Power(dzdt * dt, 2));

      /*** The particle gains back the energy it lost ***/
      const double eLoss = fStoppingPower(ekin) * state.fDensityScale * step;

      const double beta = TMath::Sqrt(dxdt * dxdt + dydt * dydt + dzdt * dzdt) / kSpeedOfLight;
      const double ekinNoLoss = fA * 931.494 * 0.5 * beta * beta;
      ekin += eLoss;

      const double speedUp = TMath::Sqrt(ekin / ekinNoLoss);
      dxdt *= speedUp;
      dydt *= speedUp;
      dzdt *= speedUp;

      dt = -fDzStep / dzdt / fIntegrationSteps;
   }

   if (!foundVertex)
      LOG(warn) << "AtMCQMinimization: the backward extrapolation did not find a minimum distance of approach "
                << "to the beam axis in " << kMaxSteps << " steps!";
}

double AtMCQMinimization::Chi2Q()
{
   double chi2Charge = 0.;
   double chi2Z = 0.;
   int numPoints = 0;

   for (size_t pad = 0; pad < fQSim.size(); pad++) {

      if (fQSim[pad] <= kChi2Threshold || fQExp[pad] <= kChi2Threshold)
         continue;

      const double deltaQ = fQSim[pad] - fQExp[pad];
      const double sumQ = fQSim[pad] + fQExp[pad];
      const double deltaZ = fZSim[pad] - fZExp[pad];

      const double pullQ = deltaQ / (sumQ * fSigmaQ);
      const double pullZ = deltaZ / fSigmaZ;

      // Pads whose center of gravity is too far away are not part of the same track
      if (pullZ * pullZ >= kMaxPullZ2)
         continue;

      numPoints++;
      chi2Charge += pullQ * pullQ;
      chi2Z += pullZ * pullZ;
   }

   fChi2Points = numPoints;

   // Too few pads for the comparison to be meaningful
   if (numPoints < kMinChi2Pads)
      return kChi2NoPads;

   chi2Charge /= numPoints;
   chi2Z /= numPoints;

   // NB: the normalization of the original implementation, which favours the trials reproducing
   // the largest number of pads (the objective function goes as 1/n^3).
   return (chi2Charge + chi2Z) / (4. * numPoints * numPoints);
}

double AtMCQMinimization::Chi2Pos(int iteration, int numExpPoints) const
{
   const int numTB = std::min(
      {std::max(fLastTB, numExpPoints), static_cast<int>(fHitsByTB.size()), static_cast<int>(fTBTrack.size())});

   double chi2 = 0.;
   int numPoints = 0;

   for (int iTB = 0; iTB < numTB; iTB++) {

      const auto &hitsAtTB = fHitsByTB[iTB];
      if (hitsAtTB.empty())
         continue;

      // Center of gravity of the hits of this time bucket
      double cmsX = 0.;
      double cmsY = 0.;
      double totalCharge = 0.;
      for (const auto *hit : hitsAtTB) {
         cmsX += hit->GetPosition().X() * hit->GetCharge();
         cmsY += hit->GetPosition().Y() * hit->GetCharge();
         totalCharge += hit->GetCharge();
      }
      if (totalCharge == 0.)
         continue;

      const double diffX = cmsX / totalCharge - fTBTrack[iTB].X();
      const double diffY = cmsY / totalCharge - fTBTrack[iTB].Y();
      const double chi2Point = (diffX * diffX + diffY * diffY) / kPosSigma2;

      // The contribution of a point is clamped so that a single outlier, or a time bucket the
      // simulated track never reached, cannot drive the fit. The clamp is loosened while the
      // steps of the random walk are still large.
      if (iteration > 2)
         chi2 += std::min(chi2Point, kPosChi2Cap);
      else
         chi2 += std::min(chi2Point, kPosChi2CapCoarse);

      numPoints++;
   }

   if (numPoints == 0)
      return 1.e10;

   // NB: same 1/n^3 normalization as the charge term
   return chi2 / (TMath::Power(numPoints, 3) * 2.0);
}

double AtMCQMinimization::Chi2Range() const
{
   if (!fHasExpEndPoint || fSimTrack.empty() || fChi2Points == 0)
      return 1.e10;

   /* NB: the original implementation took the first hit of the array as the end of the track,
    * which held when the hits came sorted by descending time bucket. The hits of an AtTrack are
    * in the order the pattern recognition found them, so the end point is looked up instead. */
   const auto &endPoint = fExpEndPoint;
   const auto &simEndPoint = fSimTrack.back();

   const double dist =
      TMath::Sqrt(TMath::Power(endPoint.X() - simEndPoint.X(), 2) + TMath::Power(endPoint.Y() - simEndPoint.Y(), 2) +
                  TMath::Power(endPoint.Z() - simEndPoint.Z(), 2));

   return dist / (fSigmaZ * fSigmaZ * fChi2Points * fChi2Points);
}

double AtMCQMinimization::GetEnergy(double A, double Z, double bRho)
{
   const double amu = 931.5;
   const double x = bRho / 0.1439 * Z / A;
   return TMath::Sqrt(2. * amu * x * x + amu * amu) - amu;
}

AtMCQMinimization::XYZPoint AtMCQMinimization::TransformIniPos(const XYZPoint &pos) const
{
   const double xDet = pos.X() * TMath::Cos(fThetaPad) + pos.Y() * TMath::Sin(fThetaPad);
   const double yDet = -pos.X() * TMath::Sin(fThetaPad) + pos.Y() * TMath::Cos(fThetaPad);
   const double zDet = pos.Z();

   const double xSol = xDet;
   const double zSol = zDet * TMath::Cos(fTiltAng) + yDet * TMath::Sin(fTiltAng) +
                       fZPadPlane / 10. * TMath::Power(TMath::Sin(fTiltAng), 2);
   const double ySol = (yDet + (fZPadPlane / 10. - zSol) * TMath::Sin(fTiltAng)) / TMath::Cos(fTiltAng);

   const double zCmm = zSol;
   const double xCmm = xSol + zCmm * TMath::Sin(fThetaLorentz) * TMath::Sin(fThetaRot);
   const double yCmm = ySol - zCmm * TMath::Sin(fThetaLorentz) * TMath::Cos(fThetaRot);

   return {xCmm, yCmm, zCmm};
}

AtMCQMinimization::XYZPoint AtMCQMinimization::InvTransIniPos(const XYZPoint &pos) const
{
   const double xSol = pos.X() - pos.Z() * TMath::Sin(fThetaLorentz) * TMath::Sin(fThetaRot);
   const double ySol = pos.Y() + pos.Z() * TMath::Sin(fThetaLorentz) * TMath::Cos(fThetaRot);
   const double zSol = pos.Z();

   const double xDet = xSol;
   const double yDet = -(fZPadPlane / 10. - zSol) * TMath::Sin(fTiltAng) + ySol * TMath::Cos(fTiltAng);
   const double zDet = zSol * TMath::Cos(fTiltAng) - ySol * TMath::Sin(fTiltAng);

   const double xPad = xDet * TMath::Cos(fThetaPad) - yDet * TMath::Sin(fThetaPad);
   const double yPad = xDet * TMath::Sin(fThetaPad) + yDet * TMath::Cos(fThetaPad);

   return {xPad, yPad, zDet};
}

void AtMCQMinimization::PrintResult() const
{
   LOG(info) << "AtMCQMinimization: minimum at theta " << fFitPar.fTheta * TMath::RadToDeg() << " deg, phi "
             << fFitPar.fPhi * TMath::RadToDeg() << " deg, radius " << fFitPar.fRadius << " mm, B " << fFitPar.fB
             << " G, energy " << fFitPar.fEnergy << " MeV/u";
   LOG(info) << "AtMCQMinimization: track starts at (" << fFitPar.fPos.X() << ", " << fFitPar.fPos.Y() << ", "
             << fFitPar.fPos.Z() << ") cm, vertex at (" << fFitPar.fVertexPos.X() << ", " << fFitPar.fVertexPos.Y()
             << ", " << fFitPar.fVertexPos.Z() << ") cm with " << fFitPar.fVertexEner << " MeV";
   LOG(info) << "AtMCQMinimization: chi2 " << fFitPar.fChi2 << " (charge " << fFitPar.fChi2Q << ", range "
             << fFitPar.fChi2Range << ") over " << fChi2Points << " pads";
}

} // namespace MCMinimization
