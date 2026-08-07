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

## Known issue

With tilt = 6.47 deg and E = 15 kV/m, the Lorentz-corrected beam in run 128 reads
**8.04 deg** instead of 6.47 -- the correction over-rotates by ~1.6 deg. Ruled out: the
drift geometry (matching would need a 569 mm gap, but the hit z-profile reaches ~1200 mm)
and pile-up in the beam cluster (run 128 beam clusters are statistically identical to run
100's). Prime suspect is the **drift-time zero**: `TB0 = 2.2` here was derived from
`AtPSA::CalculateZGeo`, which is the same formula whose `TBEntrance = 280` contradicts the
hardcoded `-271.0` in the deprecated `AtPSA::RotateDetector`, and the shear scales directly
with it.

v_D and the tilt are unaffected -- both come from B=0 data where there is no shear at all.

Note also that `AtPSASimple2.cxx` still has the Lorentz block commented out, so none of
this reaches the reconstructed positions yet.

## Requirements

The Python needs `scikit-learn` (for `sklearn.cluster.HDBSCAN`), numpy and scipy; the
system python has none of these, so run it with the Spyral environment:
`~/attpc_spyral_1.1.1/bin/python`. Paths in the scripts are absolute for this machine.
