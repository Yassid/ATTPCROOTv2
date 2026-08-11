# Dec 2014 alpha calibration

Determines the drift velocity and the detector tilt from the alpha+alpha elastic
opening angle, which is **90 deg** for equal masses (non-relativistic: the beam is
1.95 MeV/u, beta ~ 0.06, so the correction is well under a degree).

## Result

| parameter | value | from |
|---|---|---|
| drift velocity | **2.251 +- 0.011 cm/us** | run 100, cross-checked by run 128 (2.253) |
| detector tilt | **6.47 +- 0.04 deg** | run 100 |

Both are written into `parameters/ATTPC.alpha_150torr.par` and
`parameters/ATTPC.alpha_300torr_noB.par`. The tilting was set once for the whole
experiment, so 6.47 deg applies to every Dec 2014 run.

## Why run 100

Runs 96-113 were taken with the **magnet off**. With B=0 the Langevin drift vector
collapses to `v_D z-hat`, which removes both of the things that otherwise limit this
measurement:

* no ExB lateral drift, so the beam direction measures the mechanical tilt directly,
  with no degeneracy against the electric field;
* no track curvature, so a straight-line fit over a whole arm *is* the tangent at the
  vertex -- with B on, a 300 mm arc bends 28-49 deg and the chord/tangent difference
  swamps the effect being measured.

The opening-angle peak is RMS 6.2 deg at B=0 versus 11.2 deg with B on.

The reduced field is the same for both settings (30 kV/299.5 torr vs 15 kV/150 torr =
100.2 vs 100.0 V/m/torr), so the drift velocity should agree between them -- which is
what makes run 128 a real cross-check rather than an assumption.

## Pipeline

    unpack_noB.sh 100                  # unpack a magnet-off run (hits only)
    # dump hits to CSV (event,tb,x,y,q) with a short ROOT macro, then:
    calibrate_run.py all_run0100.csv --B 0.0 --E 30000

`batch_scan.py` + `fit_params.py` + `plot_opening.py` are the equivalent chain for a
magnet-on run, where the Lorentz correction has to be modelled.

## Clustering

`cluster_tracks.py`. Plain HDBSCAN on (x,y,z) merges the two outgoing arms exactly at
the vertex, which is the one place the direction matters. What separates them there is
not position but *direction*, so each hit gets a local direction from a PCA over its k
nearest neighbours and the clustering runs in [position, orientation] space. A line
direction is sign-ambiguous, so orientation enters as the sign-invariant tensor d (x) d
rather than as d. The vertex is then recovered as the intersection of the fitted arms,
not as a cluster of its own.

Two failure modes show up as distinct populations in the opening-angle spectrum and are
useful rather than harmful: **~180 deg** is one track split longitudinally into two
halves, **~0-10 deg** is two parallel tracks, i.e. beam pile-up.

## Resolved: the Lorentz correction over-rotated (fixed in 4bfcdb5f)

This section previously recorded an unexplained 1.6 deg over-rotation and named the
drift-time zero as the prime suspect. **Both of those were wrong**, and the record is kept
here because the wrong reasoning is instructive.

The drift-time zero can never explain an *angular* discrepancy: it enters every hit
identically, so it is a rigid translation. It moves the vertex, not a direction. Verified
by scanning it -- the reconstructed beam angle is identical to four decimal places for any
origin between -100 and +200 tb.

The actual cause was a missing frame rotation. `CalcLorentzVector` solves the Langevin
equation in the **field** frame, but the result was applied straight to **pad**
coordinates, which are rotated from it by `ThetaPad` (~111 deg). The shear had roughly the
right magnitude pointing ~110 deg away from where it belonged, which is why adjusting E or
the tilt never helped -- those scale it without redirecting it. Fixed in 4bfcdb5f; the
drift-time zero was separately wrong (280 vs 320) and fixed in d20536aa, but only ever
affected positions.

Note also that the paper's specialisation of the Langevin solution assumes B lies in the
y-z plane. That is an assumption about the tilt's *azimuth* and is false here -- the beam,
hence B, sits at about -162 deg in the pad plane. `ThetaRot` now carries that azimuth.

**Caveat that remains:** `ThetaRot` and `ThetaPad` are not individually identified. Every
combination of azimuth, handedness and pad rotation giving the same net shear direction
fits equally well, so changing `ThetaPad` alone will silently break the correction.

## Requirements

The Python needs `scikit-learn` (for `sklearn.cluster.HDBSCAN`), numpy and scipy; the
system python has none of these, so run it with the Spyral environment:
`~/attpc_spyral_1.1.1/bin/python`. Paths in the scripts are absolute for this machine.

## Production: runs 128-139

The first runs reconstructed with the full chain (Lorentz shear in the pad frame,
TBEntrance = 320, v_D = 2.251, tilt 6.47 deg). All 12 runs `rc=0`, 9216 pads/event,
entry counts exact, ROOT files closed cleanly. See `plots/`.

| check | production (12 runs) | reference |
|---|---|---|
| beam polar angle | **6.42 +- 0.08 deg** | 6.47 +- 0.04 (run 100, B=0) |
| alpha+alpha opening peak | **89.66 +- 0.31 deg** | 90 (elastic, equal masses) |
| vertices inside the drift volume | 89-92% | -- |

Ten of those runs were never used to tune anything -- the correction was calibrated on
128 and validated on 130 -- so the run-to-run consistency is the substantive result.

### Two estimator traps, both of which bit during this verification

**Do not fit hit positions with a 3D total-least-squares line** when you want an angle.
TLS weights all three axes equally, but the pads are 8 x 12 mm while the time coordinate
is fine-grained, so transverse pad scatter is attributed to the track direction and the
angle is biased *up*: 6.89 vs 6.40 deg on the same run-128 beams. Regress x and y **on**
drift time instead.

**Aggregate slope components, not angles.** `hypot(sx, sy)` is positive definite, so
per-event noise inflates every individual angle and a median over them is biased up
(6.96 vs 6.42 deg). Take the median of `sx` and `sy` separately and form the angle once.

### One check that is not valid with the magnet on

The beam's closest approach to the detector axis is *not* a cross-check for the B-on runs.
With B = 0 the beam flies straight and passes through the axis at the entrance window
(2.2 mm miss -- that is how TBEntrance was calibrated). With the solenoid on the beam is
focused and steered, so where it approaches the axis is a property of the beam and not of
the detector. The production value (1157 +- 6 mm) is stable run to run but must not be
compared against the 1000 mm expected from the geometry.
