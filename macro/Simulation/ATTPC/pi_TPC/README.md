# Pions in the AT-TPC

Single-particle pion sandbox for UKF development. No reaction generator —
charged pions are spawned at a fixed vertex inside the drift volume and
shot isotropically. Used to characterise UKF performance for HELIOS-style
symmetric polar acceptance.

## Layout

- `pi_TPC_sim.C`, `run_digi_attpc.C`, `run_ukf_only.C` — main pipeline.
- `run_ukf_truthpra.C` — bypasses `AtPRAtask`, seeds UKF with truth-derived
  `(R, θ, φ, GeoCenter)` evaluated at the seed-cluster MC point. Use to
  measure intrinsic UKF resolution (no PRA noise).
- `run_ukf_display.C` — ROOT Eve display.
- `analysis/` — all post-processing macros (resolution tables, residual
  plots, attrition diagnostics, PRA seed-quality vs MC truth).
- `analysis/make_performance.C` — **the** one-shot 6-panel summary plot.
  Reads `data/output_ukf_only.root` (PRA+UKF) and, if present,
  `data/output_ukf_truthpra.root` (intrinsic) — overlays them.
- `analysis/make_kinematics.C` — (θ_lab, KE) scatter overlay: truth (gray),
  PRA-seeded UKF (blue), truth-seeded UKF ceiling (red).
- `data/` — fresh outputs at top level; old plots in `data/archive/`.

## Run order

```bash
# 1. Geant4 transport (default: pi+, KE 5–50 MeV, 2000 events)
root -b -q 'pi_TPC_sim.C(2000, +1)'
#  pi-:  root -b -q 'pi_TPC_sim.C(2000, -1)'

# 2. Digi + PSA + PRA (square-pad map, 200×200 mm², 2 mm pitch)
root -b -q 'run_digi_attpc.C(2000)'

# 3. UKF fit
root -b -q 'run_ukf_only.C(2000, +1)'

# 4. (optional) intrinsic-UKF baseline with truth seeds
root -b -q 'run_ukf_truthpra.C(2000, +1)'

# 5. Performance summary
root -b -q analysis/make_performance.C
```

## Generator

- Vertex: fixed at (0, 0, 500 mm) — i.e. axis, mid-drift.
- Direction: isotropic in cos θ ∈ [cos 175°, cos 5°], φ ∈ [0, 360°].
- KE: uniform [5, 50] MeV.
- HELIOS workaround: `AtVertexPropagator::SetRndELoss(1e30)` disables the
  beam-reaction trigger that otherwise stops the primary at first dE
  deposit (single-particle, no reaction generator).

## Configuration

- B-field: 0.5 T along +z (`AtConstField`).
- Gas: P10 @ 1 bar, ρ = 1.654 × 10⁻³ g/cm³ (`geometry/ATTPC_P10_1bar.root`).
- Pad plane: programmatic `AtTpcSquareMap(2 mm, 100, 100)` — uniform 2 mm
  square pads on a 200 × 200 mm² active area, 10000 channels.
- UKF: `MeasurementSigma=1`, `MomentumSigmaFrac=0.1`, energy straggling on,
  `MinClusters=5`, `NIterations=3`, `BackExtrapMaxPath=250 mm`,
  `UseClusterDirSeed=true` (needed for symmetric forward/backward — see
  the project note `project_attpcroot_pra_backward.md`).

## Particle / energy-loss notes

- PDG ±211; rest mass 139.57039 MeV/c²; cτ ≈ 7.8 m → decay-in-flight in
  the TPC volume is rare but non-zero (Geant4 handles it).
- `AtELossCATIMA` is nuclear-projectile-oriented. We pass A=Z=1 (unit
  charge, correct Bethe-Bloch Z² scaling) and `mass_amu = m_pi/u` so the
  internal T/m argument yields the correct β. Shell/density corrections
  are not pion-tuned — adequate for first-pass UKF kinematics, revisit if
  few-percent dE/dx accuracy matters.

## Performance (2026-05-09, 2000 evts, isotropic, post-fix)

| θ_MC (deg) | yield | σ/E (PRA) | σ/E (truth) | σ_θ |
|------------|-------|-----------|-------------|-----|
| 25–155     | 99–100% | **9.9–15%** | **2–7%**     | 0.4–0.9° |
| 5–25, 155–175 | 88% | 17–25% | 10–15% | 0.5° |

Bias |<1| MeV across all bins. Polar profile is symmetric (backward is
slightly *better* than forward thanks to longer track length).

The PRA-vs-truth gap (~3× on σ/E) is dominated by `AtPRA::SetTrackInitial
Parameters` radius-fit precision (σ(R)/R ≈ 10%; intrinsic geometry floor
~5% for our chord/sagitta SNR). See `analysis/pra_seed_quality.C`.
