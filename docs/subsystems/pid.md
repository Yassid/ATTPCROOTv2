# Particle Identification (PID)

`AtPIDTask` (in `AtReconstruction/`) turns every pattern-recognized track into PID observables: a magnetic rigidity
(Bρ) and a specific energy loss (√(dE/dx)). The estimators live in `AtTools/AtPID/`. Each track gets one entry per
estimator in the `AtPIDEvent` branch; the vectors are parallel, but match on `trackID`, never on array position.

## Estimators

| Estimator | Result type | How it measures | Strengths | Known failure modes |
|-----------|-------------|-----------------|-----------|---------------------|
| `AtPIDEstimator` (legacy) | `AtPIDResult` | total charge / arc length | — | inflated dE/dx; kept only for comparison (`GetClassic()`) |
| `AtSpyralPID` | `AtSpyralResult` | Spyral method: hits sorted in z, spline along the track, first-arc circle; dE/dx over the first arc on the small pads | sharpest dE/dx on well-measured tracks; the default (`GetDefault()`) | tracks nearly flat in z (θ ≈ 90°): the z-ordered spline is ill-posed and the radius can come out much too large (see below); singular spline on z ties (`SetZTieTolerance`); `SetMinPoints` default 30 rejects many tracks |
| `AtHelix3D` | `AtHelix3DResult` | hits ordered in 3D (minimum-spanning-tree backbone), helix fit at each end, vertex end by a 2-of-3 vote (rotation sense, Spyral direction, Bragg rise); dE/dx = charge / helix path over the first arc on the small pads | does not depend on z ordering; robust direction; supplies range, containment and charge-profile flags | its circle averages a radius that shrinks along a slowing track, so Bρ of stopping tracks is biased low (more for heavier species) |

`AtPIDTask` runs all three; AtHelix3D can be switched off with `SetRunHelix3D(false)`. AtHelix3D's lab-frame outputs
use z_lab = zSign·z_digi + zOffset (`SetZPadPlane`); its internal centroids are in the same transformed frame.

### AtHelix3D options

| Setter | Default | Purpose |
|--------|---------|---------|
| `SetHitRefit(circle, alsoZ)` | off | refit the vertex-end circle (and z(φ)) on the raw window hits instead of the centroids; better precision |
| `SetArcAngle(on, frac)` | off, 0.5 | polar angle from z_lab vs arc length round the vertex-end circle over the first `frac` of the window hits → `thetaArc`, `brhoArc = B·R/sin θ`. Unbiased on simulation (0.2–0.3 % in sin θ); use with `SetHitRefit(true, true)` |
| `SetContainment(rMax, zMin, zMax)` | 145, 10, 955 mm | `contained` = far-end centroid inside the gas. The a1975 values: tracks leaving the drift window end at z_lab ≈ 964 mm otherwise look contained. **Re-measure per setup** |
| `SetMidPeak(fpkMax, minAfter, tailMax)` | 0.75, 30 mm, 0.35 | `midPeak` = the charge profile (10 mm bins from the vertex end) peaks before 75 % of the track, > 30 mm before its end, and the mean after the peak is < 35 % of it: the particle stopped and the candidate carries other hits |
| `SetEndWindowCap`, `SetAxisConstraint`, `SetRadiusGuard`, `SetBlobRemoval`, `SetBraggVote`, `SetExpectedSense` | see header | vertex-end window and robustness; `SetExpectedSense` must be measured for each field/charge/geometry |

Always computed when the backbone exists: `length` (range of a stopping track), `rStop`, `zLabStop`, `peakFrac`,
`tailRatio`, `contained`, `midPeak`.

## Calibration

1. **Per-run gain** (`SetDedxGain(gSpyral, g3D)`, factors from `AtTools::AtPIDGainTable`, one line per run). The
   AT-TPC gain drifts over a campaign (a1975: −20 %); without it every band is smeared and a fixed gate is off-centre
   for some runs. Produce the table from the data (a1975: `pid_gain_calib.C`, `pid_gain_norm.C`). A run missing from
   the table returns −1 so a production fails loudly.
2. **Pad-gain map** (`SetPadGainMap(file, column)`, read by `AtPIDGainTable::LoadPadGains`, `pad value...` per line).
   Every hit charge is divided by its pad's gain before any estimator runs (on a copy of the track; the pattern event
   is untouched). a1975: from each track's hit amplitudes over its median (`pad_gain_map_dp.C`); narrows the bands
   by 4–7 %.
3. **Saturation.** GET pulses saturate; with AtPSAMultiFit a clipped pulse is extrapolated (a1975 hit amplitudes
   pile up between ~3800 and ~4850 ADC). Slow, heavily ionizing tracks lose charge in their PID segment (a1975: ~20 %
   for deuterons at Bρ 0.15–0.20). The simulation does not clip unless `AtPulse::SetADCCeiling` is set, and the right
   value must come from the data (per-pad clip level, PSA response); until then simulated PID resolution is
   optimistic for high dE/dx.

## Combined PID (`SetCombine`)

`AtPIDTask::SetCombine(true, opt)` fills `AtPIDEvent::GetCombined()` (one `AtPIDCombinedResult` per track) with
`AtTools::CombinePID`: the Spyral result, unless its Bρ fails a cross-check against the AtHelix3D circle, in which
case the 3D result is used (`source` 0 = Spyral, 1 = 3D, −1 = none; `vetoed` marks the cross-check failure).

The cross-check (`AtPIDCombineOpt`): veto Spyral when ln(Bρ_Spyral / Bρ_3D) ≥ `lnRatioMax` and
(Bρ_Spyral < `maxBrho` or Spyral √(dE/dx) > `highSq`). The 3D Bρ is `brhoArc` when computed (`useArc`), and the 3D
√(dE/dx) is multiplied by `addBackScale` to put it on the Spyral scale. `fallback3D` uses the 3D result when Spyral
has none.

Why: on a1975 the slow elastic deuterons near 90° are nearly flat in z (≈ 10 mm); Spyral's radius comes out ≈ 1.6×
too large and they form a fake "triton" blob and streak (14 % of tracks). The circle measures them correctly, so
they are re-added rather than lost.

Recommended configuration (the a1975 D2 recipe; the last two values are gain/detector specific):

```cpp
auto pid = new AtPIDTask();
pid->SetMinPoints(15);                          // the default 30 rejects many tracks; knee at 15-20
pid->SetDedxGain(gSpyral, g3D);                 // per run, from AtPIDGainTable
pid->SetPadGainMap("padgain.txt");              // optional
pid->GetHelix3D().SetHitRefit(true, true);
pid->GetHelix3D().SetArcAngle(true);
AtTools::AtPIDCombineOpt opt;                   // lnRatioMax 0.25, maxBrho 0.45
opt.highSq = 30;                                // a1975 D2, Spyral run+pad-gain units
opt.addBackScale = 0.965;                       // a1975 D2: median Spyral/3D sqrt(dE/dx) where the Brho agree
pid->SetCombine(true, opt);
```

Validated (2026-11-02, a1975 run 0026, 2670 tracks): `brhoArc`, the 3D √(dE/dx), the veto, `midPeak` and
`contained` reproduce the analysis macros track by track (`validate_pid_framework_dp.C`). The a1975 quality cuts were
derived on raw charges; with the pad map `midPeak` changes on 4.6 % of tracks and `contained` on 0.7 %.

## Track-quality cuts (analysis level)

These depend on the species expected and stay in the analysis, built on the stored flags:

- **Too short:** a stopping track (`contained`) whose `length` is shorter than the lightest-range species allowed at
  its Bρ (a1975: e^−0.4 × the simulated triton range) cannot be any Z = 1 or 2 ion — at fixed Bρ the range scales as
  z²·Bρ⁴/M³, so tritons are the shortest. These are track-finder artefacts.
- **midPeak in the gate of interest:** removes stopped particles with other hits attached. Apply it only in the region
  being purified (a1975: the triton region); applied everywhere it also removes real slow deuterons (17 % at
  Bρ 0.10–0.15).

a1975 reference: `plot_pid_final_cuts_dp.C` (2.4 % of tracks removed; deuteron peak kept 99–100 %).

## Choosing for a new experiment

1. **Simulate every species you must separate** (not only the ejectile of interest) and run the estimators on them:
   separation and misidentification per Bρ slice, against truth (a1975: `plot_pid_dt_sim_dp.C`). A comparison
   against a computed locus misleads when an estimator biases both species alike.
2. **Calibrate the gain per run and per pad on the data** before comparing estimators.
3. **Look at the PID plane of each estimator and the tracks one removes or moves** (event displays). Spyral is the
   default; switch on the combination if Spyral's Bρ disagrees with the 3D circle for a population (typically slow,
   flat, contained tracks).
4. **Use range** for stopping tracks: Bρ vs `length` separates species independently of the gain (a1975:
   `plot_pid_range_data_dp.C`); it is the best check that a dE/dx band is real.
5. The PID limit is the dE/dx resolution: on a1975, genfit Bρ (with either hypothesis), a truncated mean, averaging
   Spyral and 3D dE/dx, and outer pads gave no PID gain. Genfit is for the kinematics once the species is known.

## Pitfalls

- AtHelix3D internal coordinates (centroids) are in z' = zSign·z + zOffset; compare hits in that frame.
- A ROOT `TTree::Draw` formula carrying per-slice lookups can exceed TTreeFormula's operator limit ("Too many
  operators") and silently return nothing; compute such flags in C++.
- Genfit output trees can repeat an entry; match fits to PID results by `trackID` and skip copies with an empty
  `AtPIDEvent`.
- `AtPIDEvent` v3 adds `fCombined`; older files read back with it empty.
