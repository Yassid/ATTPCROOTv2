# Pions in the AT-TPC

Single-particle pion sandbox for UKF development. No reaction generator —
charged pions are spawned uniformly inside the drift volume and shot
isotropically.

## Run order

```bash
# 1. Geant4 transport (default: pi+, KE 5–50 MeV, 1000 events)
root -b -q 'pi_TPC_sim.C(1000, +1)'
#  pi-:  root -b -q 'pi_TPC_sim.C(1000, -1)'

# 2. Digitization + PSA + PRA
root -b -q run_digi_attpc.C

# 3. UKF fit
root -b -q 'run_ukf_only.C(1000, +1)'
```

## Particle / energy-loss notes

- PDG ±211; rest mass 139.57039 MeV/c²; cτ ≈ 7.8 m → decay-in-flight in
  the TPC volume is rare but non-zero (Geant4 handles it).
- `AtFitterUKF` is mass-configurable; just pass `mass_pi_MeV` and the
  signed charge.
- `AtELossCATIMA` is nuclear-projectile-oriented. We pass A=Z=1 (unit
  charge, correct Bethe-Bloch Z² scaling) and `mass_amu = m_pi/u` so the
  internal `T/m` argument yields the correct β. Shell/density corrections
  are not pion-tuned — adequate for first-pass UKF kinematics, revisit if
  few-percent dE/dx accuracy matters.

## Vertex / direction

- Vertex: uniform in box (|x|,|y| < 15 cm, z ∈ [10, 90] cm) — inscribed in
  the R=25 cm, L=100 cm drift cylinder with margin from the walls.
- Direction: isotropic (θ ∈ [0, 180]°, φ ∈ [0, 360]°).
