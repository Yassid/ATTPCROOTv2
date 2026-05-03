# PUMA UKF Calibration Notes (2026-05-03 / 04)

This document records the chain of fixes that brought the PUMA single-pion
sample from σ_xy ≈ 2 mm / σ_z ≈ 14 mm with a +8 mm Δz median bias to
σ_xy ≈ 2 mm / σ_z ≈ 7 mm with a Δz bias of −0.4 mm (essentially zero).
All changes are opt-in via `AtFitterUKF` setters so other experiments
(16C+p, e20009, …) keep their existing calibrations.

## Single-pion vertex resolution: before → after

| metric | original | final |
| --- | --- | --- |
| Δx σ_IQR | 2.7 mm | 2.2 mm |
| Δy σ_IQR | 2.6 mm | 2.4 mm |
| Δz σ_IQR | 13.7 mm | **6.9 mm** |
| Δz median | +7.4 mm | **−0.4 mm** |
| median \|Δr\| | 12.4 mm | **5.5 mm** |

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
