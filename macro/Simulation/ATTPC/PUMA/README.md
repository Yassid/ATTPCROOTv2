# PUMA-channel sandbox

ATTPCROOT port of the PUMA simulation primary generator
(`puma-tpc-simulation/src/ExN03PrimaryGeneratorAction.cc`,
default `ReactionPath = 1`, `branch = 9`):

```
^3He + p-bar  ->  K+  +  K+  +  pi-  +  hexaquark
```

The generator (`AtPUMAGenerator`) is a `FairGenerator` subclass — no Geant4
dependency, uses ROOT `TGenPhaseSpace` + `TRandom`. Vertex sampling
matches PUMA's `Inputs/InputAll.txt` defaults: Gaussian (x,y) with sigma 4 mm
clipped to TrapRadius=10 mm, uniform z in [-22.5, +22.5] mm. p-bar momentum
is Gaussian magnitude * isotropic direction (defaults to at rest).

## Run order

```bash
# 1) Geant4 transport
root -b -q 'PUMA_sim.C(1000)'

# 2) Digi (square 2x2 mm pads, P10 gas — same as pi_TPC sandbox)
root -b -q 'run_digi_attpc.C(1000)'

# 3) UKF — needs per-track configuration (different particles per primary).
#    The current run_ukf_only.C in pi_TPC/ assumes a single particle type.
#    For PUMA we have 3 charged particles per event (2 K+, 1 pi-) plus the
#    invisible hexaquark; pick the channel of interest and configure mass/charge
#    accordingly (or extend AtFitterTask to handle per-track particle hypotheses).
```

## Caveat — hexaquark

Hexaquark uses PDG code 1600000 (custom, mass = 2.0 GeV) defined in PUMA's
`MinosModularList.cc`. It is not registered in Geant4's particle table by
default. The generator still pushes it into the primary stack (so analysis
sees it), but Geant4 will warn / skip transport for unknown PDG codes.
To transport it, register it via a custom physics list before `run->Init()`.

## Tunables

```cpp
puma->SetVertexXY(meanX, sigmaX, meanY, sigmaY); // mm
puma->SetVertexZHalfRange(deltaZ_mm);
puma->SetTrapRadius(r_mm);                        // mm
puma->SetPbarMomentum(mean_GeV, sigma_GeV);
puma->SetTargetMass(mass_GeV);                    // default 2.80835 (^3He)
puma->SetBeamMass(mass_GeV);                      // default 0.93826 (p-bar)
puma->SetHexaquarkMass(mass_GeV);                 // default 2.0
puma->SetHexaquarkPdg(pdg);                       // default 1600000
```
