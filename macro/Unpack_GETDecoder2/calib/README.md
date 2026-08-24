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

## alpha+alpha excitation function -- a negative result

Thick-target inverse kinematics: the beam enters at 7.8 MeV and slows continuously, so
the vertex position is an energy measurement and one run scans E_cm from 3.9 MeV to 0.
Beam energy at the vertex comes from CATIMA (the port in 926468f8), He:CO2 90/10 at
150 torr, rho = 7.27e-5 g/cm3. 1033 elastic events from runs 128-139, covering
E_cm = 0.94 to 3.90 MeV.

**What is validated.** The upper endpoint lands at **E_cm = 3.899 MeV** against the
3.900 expected from the beam energy. That is a real end-to-end check of the chain --
vertex reconstruction near the window, the entrance position, and the energy loss all
have to be right for the endpoint to come out there.

**What is NOT validated: there is no resonance signal here.** The yield peaks sharply
at E_cm = 2.98 MeV, invitingly close to the 8Be 2+ at 3.03. It is an artifact:

* the reconstructed vertex distribution peaks at z = 625 mm, and pushing that through
  the energy-loss curve gives E_cm = 2.96 MeV -- the yield peak is just the vertex
  distribution's image, not sigma(E);
* the observed structure is ~0.3 MeV FWHM while the 2+ has Gamma = 1.5 MeV, so it is
  five times too narrow to be that state.

The z = 625 mm feature is an acceptance/trigger edge, not physics. **A raw thick-target
yield is not a cross section**: it needs the detection efficiency as a function of
vertex position and scattering angle, which means a simulation. Until that exists,
no statement about resonances can be made from this data.

`excitation.py` builds E_cm from `prod_results.pkl`; `eloss_alpha_heco2.txt` is the
CATIMA table (regenerate with a short program against AtELossCATIMA).

## Energy deposition: use GetQHit(), not GetCharge()

`AtHit` carries two charge observables and they behave oppositely along a track:

* `GetCharge()` is the **peak** ADC amplitude minus baseline. Longitudinal diffusion spreads
  the pulse in time over a long drift, so the height falls for the same deposited energy --
  it is drift-distance dependent.
* `GetQHit()` is the **integral** over time buckets, which is not.

Checked on the run-128 beam tracks against the CATIMA energy-loss curve:

| observable | correlation with dE/dx |
|---|---|
| `GetCharge()` (peak) | **-0.36** |
| `GetQHit()` (integral) | **+0.78** |

The peak amplitude *anti*-correlates: it falls where dE/dx rises, because drift attenuation
runs opposite to the energy-loss trend along the beam. This initially looked like a genuine
inconsistency in the reconstruction and it is not -- it is the wrong observable. With the
integral, the beam's energy deposition follows the CATIMA curve. See
`plots/charge_vs_dedx.png`.

Use `GetQHit()` for dE/dx, particle ID, Bragg curves or energy from deposition. `GetCharge()`
remains fine for thresholds, hit weighting in fits, and picking the beam cluster by relative
brightness -- anywhere a monotonic proxy within a similar drift range is enough.
