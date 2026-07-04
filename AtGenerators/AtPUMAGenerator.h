/// @file AtPUMAGenerator.h
/// @brief FairGenerator port of PUMA's ExN03PrimaryGeneratorAction::GeneratePrimaries.
///
/// Reproduces the physics of the PUMA simulation
/// (https://github.com/.../puma-tpc-simulation) without pulling in any Geant4
/// dependency. Default channel matches PUMA's `ReactionPath = 1`, `branch = 9`:
///
///     ^3He + p-bar  ->  K+ + K+ + pi- + hexaquark
///
/// Phase-space sampled with `TGenPhaseSpace`. Vertex geometry: Gaussian (x,y)
/// constrained to a TrapRadius cylinder, uniform z in [-DeltaZ, +DeltaZ].
/// Beam (p-bar) momentum: Gaussian magnitude * isotropic direction.
///
/// Hexaquark (PDG 1600000, m=2.0 GeV) is a PUMA-custom particle. It will not
/// be transported by Geant4 unless the user registers it; the generator still
/// pushes it into the stack so analysis sees it as a primary.

#ifndef ATPUMAGENERATOR_H
#define ATPUMAGENERATOR_H

#include <FairGenerator.h>
#include <Rtypes.h>

class FairPrimaryGenerator;

class AtPUMAGenerator : public FairGenerator {
public:
   AtPUMAGenerator();
   ~AtPUMAGenerator() override = default;

   /// FairGenerator interface.
   Bool_t ReadEvent(FairPrimaryGenerator *primGen) override;

   /// Beam vertex distribution (mm).
   void SetVertexXY(Double_t meanX, Double_t sigmaX, Double_t meanY, Double_t sigmaY)
   {
      fMeanX = meanX;
      fSigmaX = sigmaX;
      fMeanY = meanY;
      fSigmaY = sigmaY;
   }
   void SetVertexZHalfRange(Double_t deltaZ_mm) { fDeltaZ = deltaZ_mm; }
   /// Centre of the uniform vertex-z distribution (mm). Vertex z drawn from
   /// [fMeanZ - fDeltaZ, fMeanZ + fDeltaZ]. Defaults to 0 (PUMA standalone),
   /// set to the geometric drift-volume centre when running inside ATTPCROOT.
   void SetVertexZOffset(Double_t meanZ_mm) { fMeanZ = meanZ_mm; }
   void SetTrapRadius(Double_t r_mm) { fTrapRadius = r_mm; }

   /// p-bar momentum distribution (Gaussian magnitude in GeV/c, isotropic direction).
   void SetPbarMomentum(Double_t mean_GeV, Double_t sigma_GeV)
   {
      fMeanMomentum = mean_GeV;
      fSigmaMomentum = sigma_GeV;
   }

   /// Override target-mass (default He-3, 2.80835 GeV/c^2).
   void SetTargetMass(Double_t mass_GeV) { fTargetMass = mass_GeV; }
   /// Override beam (anti-proton) mass.
   void SetBeamMass(Double_t mass_GeV) { fBeamMass = mass_GeV; }

   /// Hexaquark mass / PDG (PUMA defaults: m=2.0 GeV, PDG=1600000).
   void SetHexaquarkMass(Double_t mass_GeV) { fHexMass = mass_GeV; }
   void SetHexaquarkPdg(Int_t pdg) { fHexPdg = pdg; }

   /// Select the reaction channel (upstream PUMA "branch" index).
   ///   9  (default) : ^3He + p-bar -> K+ K+ pi- hexaquark (phase space).
   ///   8            : test channel — a single pi+/pi- pair emitted
   ///                  back-to-back, each with fixed total energy
   ///                  fTestEnergy (default 0.4 GeV), momentum in the
   ///                  xy-plane (pz=0, i.e. perpendicular to a z-directed
   ///                  B field).
   ///   10           : test channel — a K+/K+ pair (same charge)
   ///                  back-to-back, each with fixed total energy
   ///                  fTestEnergy (set 0.777 GeV -> p=600 MeV/c, or
   ///                  0.578 -> 300 MeV/c), pz=0. Heavier + same-charge, so
   ///                  the two kaons trace mirror circles (both curve the
   ///                  same way) and sit in the same +Brho band at higher
   ///                  dE/dx than pions.
   void SetChannel(Int_t branch) { fChannel = branch; }
   /// Per-particle total energy (GeV) for the branch-8/10 test channels.
   void SetTestEnergy(Double_t E_GeV) { fTestEnergy = E_GeV; }

private:
   /// Branch-8 back-to-back pi+/pi- test emitter. Vertex already computed
   /// (cm) by ReadEvent. Returns kFALSE and skips the event if fTestEnergy
   /// is below the pion mass.
   Bool_t GenerateBranch8(FairPrimaryGenerator *primGen, Double_t vx_cm, Double_t vy_cm, Double_t vz_cm);

   /// Branch-10 back-to-back K+/K+ test emitter (same charge). Skips the
   /// event if fTestEnergy is at/below the kaon mass.
   Bool_t GenerateBranch10(FairPrimaryGenerator *primGen, Double_t vx_cm, Double_t vy_cm, Double_t vz_cm);

   // Vertex
   Double_t fMeanX{0.}, fSigmaX{4.}, fMeanY{0.}, fSigmaY{4.};
   Double_t fDeltaZ{22.5};
   Double_t fMeanZ{0.};
   Double_t fTrapRadius{10.};

   // Beam
   Double_t fMeanMomentum{0.};   // GeV/c
   Double_t fSigmaMomentum{0.};  // GeV/c
   Double_t fTargetMass{2.80835}; // ^3He
   Double_t fBeamMass{0.93826};   // p-bar

   // Hexaquark
   Double_t fHexMass{2.0};
   Int_t fHexPdg{1600000};

   // Channel selection
   Int_t fChannel{9};          // upstream PUMA "branch" index
   Double_t fTestEnergy{0.4};  // GeV per particle, branch-8 test channel

   ClassDefOverride(AtPUMAGenerator, 2);
};

#endif
