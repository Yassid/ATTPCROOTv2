# Unscented Kalman Filter (UKF) Development Guide

## Overview

The UKF implementation reconstructs charged particle tracks in the AT-TPC detector. It propagates a 6D state vector through magnetic/electric fields and corrects against measured hit positions using the Unscented Kalman Filter algorithm.

## Architecture

The implementation is layered in three tiers:

```
kf::UnscentedKalmanFilter<DIM_X, DIM_Z, DIM_V, DIM_N>   (generic UKF)
        |
kf::TrackFitterUKF                                        (physics-specific)
        |
EventFit::AtFitterUKF                                     (FairRoot integration)
```

### Layer 1: Generic UKF (`AtReconstruction/AtFitter/OpenKF/`)

Header-only, templated Kalman filter library with no physics dependencies.

| File | Description |
|------|-------------|
| `types.h` | Eigen type aliases (`Vector<N>`, `Matrix<R,C>`) |
| `kf_util.h` | Matrix utilities, Cholesky updates, linear solvers |
| `kalman_filter/kalman_filter.h` | Abstract base class for all KF variants |
| `kalman_filter/unscented_kalman_filter.h` | Full UKF with augmented state for process/measurement noise |
| `kalman_filter/unscented_transform.h` | Standalone sigma-point transform |
| `kalman_filter/square_root_ukf.h` | Numerically stable variant using Cholesky factors |

### Layer 2: Track Fitter (`AtReconstruction/AtFitter/OpenKF/kalman_filter/`)

| File | Description |
|------|-------------|
| `TrackFitterUKF.h` | Template base + concrete `TrackFitterUKF` class |
| `TrackFitterUKF.cxx` | Physics process/measurement models, RTS smoother |

**State vector (6D):** `[x, y, z, p_mag, theta, phi]`
**Measurement vector (3D):** `[x, y, z]` (hit position)
**Process noise (1D):** energy straggling
**Measurement noise (3D):** position uncertainty

Key methods:
- `funcF()` - Process model: propagates state through B/E fields via `AtPropagator`
- `funcH()` - Measurement model: projects state to measurement space
- `predictUKF()` - Prediction step to next measurement plane
- `correctUKF()` - Correction step with observed hit
- `smoothUKF()` - Rauch-Tung-Striebel backward smoother

### Layer 3: FairRoot Wrapper (`AtReconstruction/AtFitter/`)

| File | Description |
|------|-------------|
| `AtFitterUKF.h/cxx` | Bridges `kf::TrackFitterUKF` to the FairRoot task pipeline |

Configurable parameters:
- Magnetic/electric field vectors
- UKF sigma-point weights (alpha, beta, kappa)
- Position measurement sigma (default: 1.0 mm)
- Momentum uncertainty fraction (default: 0.1)
- Energy straggling toggle

## Dependencies

- **Eigen3** - Required. The entire OpenKF module is conditionally compiled when Eigen3 is found.
- **AtPropagator** - RK4 adaptive stepper for track propagation through fields.
- **AtELossModel** - Energy loss calculation (CATIMA backend).

## Tests

| Test file | What it covers |
|-----------|---------------|
| `OpenKF/kalman_filter/kalman_filter_test.cxx` | Base KF linear/extended |
| `OpenKF/kalman_filter/unscented_kalman_filter_test.cxx` | Generic UKF |
| `OpenKF/kalman_filter/unscented_trasform_test.cpp` | Sigma-point transform |
| `OpenKF/kalman_filter/square_root_ukf_test.cpp` | Square root variant |
| `OpenKF/kalman_filter/TrackFitterUKFTest.cxx` | Track-specific UKF |
| `AtFitterUKFTest.cxx` | FairRoot integration |

Run tests:
```bash
cd build && ctest -V -R UKF
# Or directly:
./build/tests/AtReconstructionTests
```

## Macros

Integration/demo macros in `macro/tests/UKF/`:
- `UKFTask.C` - Full task pipeline demo
- `UKFSingleTrack.C` - Single track fitting
- `SimulateManyTracks.C` / `TestManyTracks.C` - Batch validation

## Development Notes

### Resolved Issues (identified via `macro/tests/UKF/UKFDiagnostics.C`)

**Bug 1: Model noise Qmod applied uniformly across dimensions** — FIXED
- `SetInitialState` now sets `m_matQmod` per-dimension using `fPosModelNoise` (positions),
  `fMomModelNoise` (momentum), and `fAngModelNoise` (angles).
- Defaults: `fPosModelNoise=1e-4`, `fMomModelNoise=1e-2`, `fAngModelNoise=1e-4`.

**Bug 2: Sigma-point weight explosion with alpha=1e-3** — FIXED
- Default `fAlpha` in `AtFitterUKF` changed from `1e-3` to `0.5`.
- Users can still override via `SetUKFParameters()`.

**Bug 3: kReplaceFirstCov breaks RTS smoother consistency** — FIXED
- Default changed from `true` to `false`. Can still be toggled: `ukf.kReplaceFirstCov = true`.

**Bug 4: GetAugStateVector hardcodes 7 elements** — FIXED
- Now loops over `DIM_A` instead of hardcoding indices.

**Bug 5: Chi2/ndf doesn't subtract fitted parameters** — FIXED
- ndf changed to `3*(nClusters-1) - 6`, accounting for 3 measurements per cluster
  and 6 fitted state parameters.

**Bug 6 (ROOT CAUSE OF SMOOTHER DEGRADATION): Stopped sigma points set momentum to zero** — FIXED
- `AtPropagator.cxx` now preserves momentum direction with small magnitude (`fStopTol`)
  when a particle stops, instead of setting momentum to `(0,0,0)`.
- Falls back to `fLastMom` direction if current momentum is already zero.
