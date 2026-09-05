# Reversed detector (`ReverseDrift`)

The beam enters through the **pad plane** and leaves through the **cathode**, instead of the other
way round. The motivation is the exit aperture: with the beam leaving through the cathode, the
opening for a downstream telescope can be made large without removing pads, because a cathode is
just a plane.

Turn it on with one line in the parameter file:

```
ReverseDrift:Int_t   1     # beam enters through the PAD PLANE
```

It is **optional and defaults to 0**, so every existing par file keeps behaving exactly as before.

---

## What actually changes, and what does not

Only one physical quantity changes: **the drift length of a given point.**

| | drift length | hit `z_digi` |
|---|---|---|
| normal | `ZPadPlane - z_beam` | `ZPadPlane - z_beam` |
| reversed | `z_beam` | `ZPadPlane - z_beam` |

`z_beam` is the distance along the beam from where it enters the active volume; `z_digi` is the z
stored on `AtHit`.

**The hit z is the same expression in both modes.** That is the entire design. The drift changes,
the digitisation maps it to a time bucket, and `AtPSA` maps it back with the complementary sign, so
the hit cloud comes out in the same frame it always did. Consequently:

* pattern recognition sees the same clouds,
* the fitters' `z_lab = ZPadPlane - z_digi` still recovers distance along the beam,
* "the vertex is the highest-`z_digi` end" still picks the beam-entry end,
* the B-field handedness, `bFieldSign`, and every angle convention are untouched,
* analysis macros that do their own `ZPadPlane - z` keep working.

What *does* differ, and is the point of the exercise:

* **diffusion** — a track near beam entry now drifts ~0 mm instead of the full length, and one at
  the far end drifts the full length instead of ~0. Transverse and longitudinal spread swap ends.
* **which time bucket** a hit lands in.

The total set of drift lengths is unchanged (still 0 … `ZPadPlane`), so the readout window
requirement is identical and nothing new can fall outside it.

## Where it is implemented — two lines, both reached from the par

1. `AtClusterize::getCurrentPointLocation` — `drift = reversed ? z_beam : ZPadPlane - z_beam`.
2. `AtPSA::CalculateZGeo` — returns `ZPadPlane - drift` instead of `drift` when reversed.

**It lives in `AtDigiPar`, not in a task setter, deliberately.** Digitisation and reconstruction
must agree about which end the pad plane is on. Two independent setters would let a job digitise
reversed and reconstruct normal, and the result of that is a plausible-looking mirrored z rather
than an error. One shared parameter makes the disagreement impossible to express.

Both `AtClusterize` and `AtPSA` print the mode on every run, and `AtDigiPar` emits a warning when
it is on, so it can never be in effect silently.

---

## Traps

### `TBEntrance` stops meaning "beam entrance"

`CalculateZGeo` uses `TBEntrance` only to locate the pad plane:
`TB(pad plane) = TBEntrance - ZPadPlane / mm-per-TB`. What the parameter really fixes is **the time
bucket of the far end of the drift volume**.

In normal running the far end *is* the beam entrance window, so the name is honest. **In reversed
running the far end is the cathode, i.e. where the beam LEAVES.** The number keeps its arithmetic
meaning and the name becomes misleading.

For data this changes the calibration procedure: the beam's own ionisation now appears at *low*
time bucket where it enters at the pad plane and at *high* time bucket where it exits, the
opposite time-order from normal running. A TB calibration that finds "the beam entrance" by
looking for the beam at high TB will find the beam **exit** instead and be wrong by the full drift
length. Calibrate on the pad plane (drift = 0) and the cathode (drift = full length) explicitly.

### `CalibrateZ`'s Spyral branch uses the opposite z convention — and always has

`AtPSA::CalibrateZ` has two paths and they do not agree, independently of this feature:

* `CalculateZGeo` returns **0 at the pad plane** (the drift length),
* the Spyral two-point branch (`SetSpyralZ`) returns **0 at the window** (distance along the beam).

Since the fitters apply `z_lab = ZPadPlane - z_digi`, which is correct for `CalculateZGeo`, any
chain enabling `SetSpyralZ` would get a **mirrored** `z_lab`. No production macro does — the only
caller is `a1975/D2_UKF/_archive/unpackPSA_compare.C`, whose own comment says it enables the branch
"so z is comparable to Spyral (no flip)".

This was left exactly as it is rather than quietly corrected, because changing it would move every
number anyone has ever taken through that branch. It is documented here so the next person to
enable it knows. `ReverseDrift` does not touch that branch: in reversed running the two frames
coincide, since beam entry *is* the pad plane.

### A z offset calibrated in the digi frame flips sign

The validation shows the per-hit residual `z_digi - (ZPadPlane - z_truth)` sitting at
**+9.4 mm normal and -9.4 mm reversed**, flat against truth z in both, with an IQR of 1.8 mm.

That offset is the PSA peaking-time bias. It is applied in time buckets, i.e. in the *drift*
domain, before the z mapping — so it is handled correctly and needs no change. But because the
reversed mapping negates the drift, the same bias appears in the digi frame with the **opposite
sign**.

The practical consequence: `AtPSAMax::SetPeakingShiftTBs` and anything else expressed in time
buckets carries over untouched, but **any correction someone has tuned as a millimetre offset on
`z_digi` must have its sign flipped** when the detector is reversed. A constant that was right in
one orientation is wrong by twice itself in the other.

### Not covered

* **`AtSpaceChargeTask`** reads `ZPadPlane` and models charge along the drift; it has not been made
  reversal-aware. No macro in the repository uses it, so nothing is broken today, but it must be
  revisited before space charge is used with a reversed detector.
* **`AtPSA::CalculateZ`** (already `[[deprecated]]`, no callers) is not reversal-aware.
* **The pad-plane material** is not in the geometry. The beam now passes through the micromegas and
  its support instead of the thin entrance window, which is a real difference in material budget
  and straggling that this flag does not model. `geometry/ATTPC_*.C` has a `tpc_window` volume at
  z = 0 and no pad-plane volume at all.

---

## Validating a change here

```
# reconstruct ONE sim twice, with pars differing only by the ReverseDrift line
root -b -q 'run_reco_Ar46_TC.C("<sim>","nrm/gs_s3001_reco.root","ATTPC.46Ar_3Hed_sim_B285.par",...)'
root -b -q 'run_reco_Ar46_TC.C("<sim>","rev/gs_s3001_reco.root","ATTPC.46Ar_3Hed_sim_B285_rev.par",...)'
root -b -q 'validate_reverse_3Hed.C("<sim>","nrm/gs_s3001_reco.root","rev/gs_s3001_reco.root")'
```

`validate_reverse_3Hed.C` is anchored on MC TRUTH, not on a comparison between the two outputs.
Each hit carries the index of the `AtMCPoint` that made it and the sim file still holds that
point, so the macro checks the invariant directly, hit by hit:

    residual = z_digi - (ZPadPlane - z_truth)     must be ~0 in BOTH modes

**Do not test this by comparing the mean hit z between the two modes.** That was the first
attempt and it produces a false failure. Diffusion spreads charge over more pads at long drift, so
more hits clear threshold at the long-drift end — and that is the opposite end of the chamber in
each mode. The mean hit z therefore legitimately moves by tens of mm while the mapping is
perfectly correct. It moves for the same reason the feature exists.

Measured on 400 events of the 46Ar g.s. sample:

| | residual median | IQR | vertex-z correlation vs truth |
|---|---|---|---|
| normal | +9.37 mm | 1.81 mm | **+0.967** (185/193 fits) |
| reversed | −9.75 mm | 1.81 mm | **+0.967** (182/193 fits) |

The residual is **flat against truth z** in both modes; a flipped mapping would slope at −2 and put
the median hundreds of mm out. The ±9.4 mm is the peaking-time bias discussed above. The
correlation column is the independent end-to-end check via `fitGenfitter_Ar46.C` +
`ex_genfit_3Hed.C`: it must stay positive and the same size in both modes.

Positive control, printed by the same macro: the charge-weighted mean hit z moves
526.9 → 477.2 mm as diffusion swaps ends. Zero movement would mean the flag never reached
`AtClusterize`.

## Normal running is unchanged by construction

Both edits are guarded so the normal path is the original expression:
`AtPSA::CalculateZGeo` early-returns before the reversal, and `AtClusterize` selects the original
formula with a ternary. Nothing in normal mode goes through a new code path, so no existing result
can move.
