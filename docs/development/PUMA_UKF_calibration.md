# PUMA UKF Calibration Notes (2026-05-03 / 04)

This document records the chain of fixes that brought the PUMA single-pion
sample from σ_xy ≈ 2 mm / σ_z ≈ 14 mm with a +8 mm Δz median bias to
σ_xy ≈ 2 mm / σ_z ≈ 2 mm with essentially zero Δz bias and +83% fit yield.
All changes are opt-in via `AtFitterUKF` setters so other experiments
(16C+p, e20009, …) keep their existing calibrations.

## Single-pion vertex resolution: before → after

| metric | original | after 2026-05-03/04 fixes | after late 2026-05-04 retune |
| --- | --- | --- | --- |
| Δx σ_IQR | 2.7 mm | 2.2 mm | 2.49 mm |
| Δy σ_IQR | 2.6 mm | 2.4 mm | 2.30 mm |
| Δz σ_IQR | 13.7 mm | 6.9 mm | **2.13 mm** |
| Δz median | +7.4 mm | −0.4 mm | +0.12 mm |
| median \|Δr\| | 12.4 mm | 5.5 mm | **4.05 mm** |
| fit yield (1000 evts, mult=1) | 1191 | 1211 | **2197** |

The final column adds the cluster-density retune (target=8 over [4, 14]
mm) and the measurement-sigma calibration (0.5 mm) on top of the prior
fixes. See "Cluster density + measurement sigma retune" below.

## How to reproduce

From `macro/Simulation/ATTPC/PUMA/`:

```bash
# 1. Single-pion gun (random charge, |p|=100–300 MeV/c, vertex at 0,0,75)
root -b -q 'PiGun_sim.C(1000)'

# 2. Digi+PSA+PRA+UKF in one job, π+/π− hypotheses only
root -b -q 'run_digi_attpc.C(1000, 8.0, true, true)'

# 3. Per-event inspection
root -l pi_inspect.C

# 4. Resolution / PID / θ-bias / digi-z-offset diagnostics
root -b -q 'vertex_recon.C("./data/output_digi.root", "./data/attpcsim.root", "./data/vertex_recon.png")'
root -b -q 'pid_confusion.C("./data/output_digi.root", "./data/attpcsim.root")'
root -b -q theta_bias.C
root -b -q digi_z_offset.C
```

## Five `AtFitterUKF` setters added (all default-off)

The PUMA opt-in cluster is at the bottom of `run_digi_attpc.C` and `run_ukf_multi.C`.

### `SetUseClusterDirSeed(true)`
Replaces PRA's `(GeoTheta, GeoPhi)` direction seed with the cluster-sequence
derivative AND re-sorts clusters by xy-distance from the beam axis (closest
first). Required for PUMA where `acos(sign·|dir.Y|)` collapses
upward- and downward-going tracks into the same θ half-sphere
(GeoTheta = 115° instead of 65° for upward π−'s).

### Tangent-from-PRA-circle direction
Inside the cluster-direction-seed branch: instead of the raw chord
`clusters[1] − clusters[0]`, the seed direction is the perpendicular to
`clusters[0] − GeoCenter` in (x,y), with the sign disambiguated by chord
alignment. Removes chord-vs-tangent bias on weakly-curved tracks. No
flag — applied whenever `fUseClusterDirSeed = true`.

### `SetUseArcWalk(true)`
Switches the in-fitter re-clustering from `ClusterizeSmooth3D` (Z-seeded,
gap-sensitive) to `ClusterizeArcWalk` (kNN spanning-tree + double-DFS arc
ordering, geometry-driven). Halved Δz σ_IQR (13.7 → 7.5 mm) by
eliminating the wrong-end-first-cluster bias on upward tracks. Tunable
via `SetArcWalkParams(minHits, kNN)`.

### `SetForceVertexOnBeamAxis(true)`
Forces the back-extrapolated vertex xy = (0, 0) instead of POCA on the
PRA circle. Tested on PUMA — the +8 mm bias was NOT from this offset, so
the flag is **left off by default**, but kept available for setups where
the vertex is known to be on the beam axis.

### `SetVertexZBias(8.6)`
Subtracts a constant from the final `initialPositionXtr.Z()`. Compensates
the digi z-reconstruction systematic discovered in the `digi_z_offset.C`
diagnostic (see "Root cause of the +8 mm Δz bias" below).

## Back-extrap also uses raw first hit, not cluster centroid

In `AtFitterUKF::GetFittedTrack`'s helix back-extrap branch, the starting
point is now the raw hit closest to the beam axis (`rxFirst, ryFirst,
rzFirst`) rather than the smoothed first-cluster state. The cluster
centroid sits ~half-cluster-spacing inside the track from the first
physical hit; for upward-going tracks (vz_first_hit < vz_centroid) this
introduces a positive Δz bias of size ~half_cluster · cot(θ).

The change is unconditional (no flag); it only affects the helix-back-extrap
path which is itself gated by `SetUseHelixBackExtrap(true)`.

## Root cause of the +8 mm Δz bias

`AtPulse` adds a pulse shaping delay of `peakingTime / TBTime` time-buckets
when generating the pad waveform (the response peaks `peakingTime` after
the electron arrives). `AtPSA::CalculateZGeo` then uses the *peak* time-bucket
as `peakIdx` without subtracting the shaping delay:

```cpp
// AtPSA::CalculateZGeo (current)
return fZk - (fEntTB - peakIdx) * fTBTime * fDriftVelocity / 100.;
```

This produces `z_recon = z_phys + peakingTime · vDrift`. For PUMA:
500 ns × 1.5 cm/μs = **7.5 mm**, plus ~1 mm from the pulse shape's
asymmetric peak → measured **8.62 mm** constant offset across the
entire drift volume (verified by `digi_z_offset.C`).

The `TBEntrance` parameter would in principle absorb this, but
`AtPulse` and `AtPSA` use the same `TBEntrance` value and the term cancels
in `(fEntTB_pulse + drift + peak) − fEntTB_psa`. So tweaking `TBEntrance`
in the par file has zero effect on `z_recon`.

The proper fix is inside `AtPSA::CalculateZGeo` (subtract `peakingTime/TBTime`
from `peakIdx` before computing z), but that touches a code path used by
every experiment with potentially calibrated empirical compensation. To
avoid breaking those, we apply the correction post-fit via
`AtFitterUKF::SetVertexZBias(8.6)` for PUMA.

## Cluster density + measurement sigma retune (2026-05-04 late)

Two layered wins on top of the 2026-05-03/04 fixes, both opt-in via
existing `AtFitterUKF` setters.

### `SetTargetClusters(8)` + `SetAdaptiveDistBounds(4.0, 14.0)`

Replaces the prior PUMA setting `target=5` over `[6, 18]` mm. Tightening
the floor to 4 mm (PUMA inner-ring pad pitch is ~2 mm) and reducing the
cap to 14 mm gives short K+ tracks ~8 clusters at 4–6 mm spacing instead
of the 1–2 they collapsed to before; long π− tracks cap out at the
14 mm distance bound to still hit the 8-cluster budget. Net effect on
the PiGun calibration sample: fit yield 1211 → 2219 (+83%), σ_z 7.5 →
2.96 mm. target=12 over [3, 10] gives identical numbers — the digi
resolution floor is reached.

### `SetMeasurementSigma(0.5)` (was 1.0)

The default 1.0 mm was overconservative for PUMA. The actual digi
resolution is pad pitch ÷ √12 ≈ 0.6 mm transverse and 0.75 mm in z
(time-bucket quantization). Tightening σ to 0.5 mm increases the
Kalman gain on every cluster update and drops σ_z(vertex) from
2.96 → 2.13 mm with no fit-yield or wall-time cost. σ = 0.3 plateaus
(only −3 % more on σ_z, per-event regresses), and `SetNIterations(5)`
also gives only −3 % σ_z at +30 % wall time — both rejected.

### Tested and not adopted

| knob | result |
| --- | --- |
| `SetUsePerClusterCov(true)` | σ_z went 2.96 → 4.16 mm — ArcWalk's diffusion-based per-cluster covariance is too pessimistic in z. |
| `SetNIterations(5)` | −3 % σ_z at +30 % wall time. Bad trade. |
| `SetMeasurementSigma(0.3)` | −3 % σ_z per-track but per-event slightly worse. Saturated. |
| `SetForceVertexOnBeamAxis(true)` on PiGun | PiGun vertex is uniform in `fTrapRadius = 10 mm`, so snapping to (0,0) biases the test. Worth re-testing on real PUMA reactions where the vertex IS on axis. |

## `KinematicsXtr` — at-first-cluster KE alongside at-vertex KE

`AtFitterUKF::GetFittedTrack` now also calls
`fittedTrack->SetKinematicsXtr(KE_first_cluster, theta_s, phi_s)` in
addition to the existing `SetKinematics(KE, theta_s, phi_s)`. Convention
in this codebase (slightly counter to the AtFittedTrack header doc) is:

- `Kinematics` ⇒ at-vertex (post back-extrapolation) state. KE here has
  the energy-loss correction along the back-extrap arc subtracted off,
  so `kin.kineticEnergy` is the KE the particle had at the vertex.
  Backward-compatible with every analysis macro that already reads it.
- `KinematicsXtr` ⇒ at-first-cluster (pre back-extrapolation) state.

θ and φ are the same in both — the smoothed-state angles are at the
first cluster and we don't update them during back-extrap. `pi_inspect.C`
displays both side by side per fit:

```
pdg=-211   KE 64.3 -> 82.4 MeV   p 145.7 -> 167.5 MeV/c   th=141.8  ph=-87.6
```

## Active PUMA UKF setting block

Mirror copy in `run_digi_attpc.C` and `run_ukf_multi.C`:

```cpp
ukf->SetMeasurementSigma(0.5);                  // digi-floor matched
ukf->SetMomentumSigmaFrac(0.1);
ukf->SetEnableEnergyStraggling(true);
ukf->SetMinClusters(2);
ukf->SetNIterations(3);
ukf->SetZPadPlane(0.0);
ukf->SetBackExtrapMaxPath(250.0);
ukf->SetUseHelixBackExtrap(true);
ukf->SetMinClusterSpacing(1.0);
ukf->SetTargetClusters(8);                      // short tracks ~8 clusters
ukf->SetAdaptiveDistBounds(4.0, 14.0);          // 2-mm pad-pitch aware
ukf->SetUseClusterDirSeed(true);
ukf->SetUseArcWalk(true);
ukf->SetVertexZBias(8.6);                       // AtPSA shaping-delay fix
```

## New macros

| file | purpose |
| --- | --- |
| `PiGun_sim.C` | random π+/π− gun, configurable multiplicity / momentum, same Geant4 transport as PUMA_sim |
| `pi_inspect.C` | TGMainFrame-based event-by-event viewer: 3 projection pads (XY, XZ, YZ) + info pad. PRA tracks colored, MC primary points colored by parent PDG, per-track truth-match labels |
| `inspect_fit.C` | one-shot event projection plot (saves PNG/PDF) |
| `theta_bias.C` | smoothed θ vs truth θ binned by truth θ; ruled out θ-shift as the +8 mm Δz bias source |
| `digi_z_offset.C` | matches each digi hit to nearest MC point in xy and reports Δz_digi-MC; identified the +8.6 mm constant additive bias |
| `track_loss.C` | tallies PRA → UKF acceptance and DROP categories |
| `pid_confusion.C` | confusion matrix vs MC truth |
| `vertex_recon.C` | Δx/Δy/Δz histograms colored by truth species |
| `brho_vs_dedx.C` | truth-labelled signed Bρ vs dE/dx PID plot (PRA-cross-product sign) |
| `brho_vs_dedx_ukf.C` | UKF-hypothesis-driven version |

## New generator

`AtGenerators/AtPiGunGenerator.{h,cxx}` — random π+/π− gun:
- Multiplicity uniform in `[fMinMult, fMaxMult]`
- Charge: random per pion (or forced via `SetCharge(±1)`)
- Direction: isotropic
- |p|: uniform in `[fMinP_GeV, fMaxP_GeV]`
- Vertex: same Gaussian-clipped-cylinder API as `AtPUMAGenerator`

## Outstanding issues

- **K vs π mass discrimination is poor** (~55% correct on pure pion sample
  with both K and π hypotheses enabled). This is a dE/dx model issue, not
  a fit architecture issue — confirmed by the pure-pion + π-only-hypothesis
  test, where mass ID is 100% within sign-correct fits. Worth investigating
  the `AtELossCATIMA` calibration for K+ at 150–250 MeV/c.
- **φ has a residual −17° median bias** (rotation between vertex and first
  cluster — known structural issue from `Kinematics.{theta,phi}` being
  stored at the first cluster, not at the vertex).
- **`AtPSA::CalculateZGeo`** still has the +7.5 mm systematic baked in
  globally. The post-fit `SetVertexZBias` is a workaround that only
  corrects the back-extrapolated vertex; clusters and smoothed positions
  retain the offset internally. Cleanest long-term fix: subtract
  `peakingTime/TBTime` from `peakIdx` in the formula and re-calibrate
  every experiment's reference resolution.
- **σ_xy floor at ~2.5 mm**: dominated by φ uncertainty at the first
  cluster acting over the ~63 mm back-extrap arc. Cluster-direction seed
  + tangent-from-circle already in use; further xy gains likely require
  PRA-side improvements (more accurate circle on short / kinked tracks).
- **`run_eve` events 71/72 look identical** despite different MC truth
  (e71: π+ θ=161°, e72: π- θ=96°). Pre-existing in `ATTPCROOTv2-OpenKF`
  baseline 72c24e44 — not introduced by these calibration commits.
  `pi_inspect.C` 2D viewer correctly shows both, so the data is fine —
  it's the Eve display path (`AtViewerManager` / `AtSidebarBranchControl`
  branch caching, or `FairRunAna::RunSingleEntry` not invalidating the
  `AtTabInfoFairRoot` cache).
