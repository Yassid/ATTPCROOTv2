# ATTPCROOT vs Spyral on the same a2091 15C+d data — 2026-11-03

Both codes reconstruct the SAME raw runs. This is what they agree and disagree on. All numbers
from `dd/plots/dd_kin_C15d.root` (ATTPCROOT, 2026-09-03 production) against
`~/Data/spyral/C15d_a2091/spyral_kin.root` (Spyral 1.1.1 + patches), IC-gated, chi2/ndf < 5,
ke > 6, vz in [0,400] unless stated.

## They agree on the physics

| | ATTPCROOT | Spyral |
|---|---|---|
| resonance centroid | **7.04 MeV** (22 sigma over sidebands) | **6.79 MeV** (13 sigma) |
| g.s. resolution, all theta | 1.42 FWHM | 0.59 FWHM |
| g.s. resolution, theta > 68 | 1.36 FWHM | 0.96 FWHM |

The ~6.9 MeV resonance is seen independently by both codes. That is the one solid piece of
physics here.

## They disagree on the beam energy by 12 MeV

Measured identically in both — the elastic peak position, where Ex = 0 by definition:

| theta cut | ATTPCROOT Ebeam | Spyral Ebeam |
|---|---|---|
| > 40 | 205 | 193 |
| > 60 | 203 | 190 |
| > 68 | 205 | 190 |

Stable within each code across the angular range, so the gap is systematic, not an artefact of
acceptance. **ATTPCROOT is internally inconsistent about it**: its own `ebeam_dd_C15d.C` ridge
method reports 186-191, and its g.s. WIDTH is minimised near 185, while the g.s. CENTROID says
205. Spyral has no such split.

For reference, the a2091 PROTON-target runs (138-182, same beam) were calibrated at 198 MeV, which
sits between the two.

## Where the difference lives: a constant 0.5 deg in theta, NOT in z

42,637 tracks matched on (run, event), comparing the most energetic track per event:

```
KE ratio (ATTPCROOT/Spyral)   mean 1.006, peak 0.988   -> the energies agree to ~1 %
dTheta   (ATTPCROOT-Spyral)   mean +0.45 deg, rms 0.79 -> a systematic angular offset
```

Binned in vertex z:

```
vz [mm]      n      dTheta     dVz     KE ratio
   0- 125  9440    +0.561     -0.6      1.007
 125- 250  9221    +0.622     +1.9      0.999
 250- 375  8668    +0.588     +1.8      1.001
 375- 500  7446    +0.514     +1.7      1.015
 500- 625  5378    +0.492     +1.3      1.032
 625- 750  2292    +0.414     +1.1      1.134
dTheta slope: -0.27 deg across the full 1000 mm
```

**The two codes agree on the vertex z to within ~2 mm everywhere.** That contradicts the
expectation in README.md that the missing Garfield drift-field correction would show up as a
z non-linearity between the codes — whatever that correction does, it is not moving z here.

What is left is a NEAR-CONSTANT theta offset of ~0.5 deg. Propagating it: at theta 58 deg,
dT/dtheta = -1.26 MeV/deg and dT/dEbeam = 0.114, so 0.5 deg of theta is worth about **5 MeV of
apparent beam energy** — a large fraction of the 12 MeV gap, but not all of it.

A constant offset is characteristic of a convention or axis definition rather than a calibration:
worth checking how each code defines the beam axis (ATTPCROOT applies a pad-plane rotation,
`ThetaPad`), before assuming either number is right.

**The KE ratio is flat at 1.00 until vz ~ 500 mm and then climbs to 1.13.** That is a separate
effect, confined to the far end of the chamber, and it is the same region where the Spyral g.s.
centroid drifts and its resolution degrades (1.63 -> 3.11 FWHM). Cutting vz < 400 avoids both.

## Open

* Which beam energy is right. 198 from the proton-target runs is the only external anchor.
* The origin of the 0.5 deg offset.
* Whether the HDBSCAN re-reco (see `rereco_hdb_C15d.sh`) changes any of this — it puts 3.4x more
  hits on each track, which is exactly what should sharpen the angle.

---

## Where the theta offset is NOT: four hypotheses tested and killed (2026-11-03)

The +0.45 deg polar-angle offset is the largest single contributor to the Ebeam disagreement, so it
was worth attacking directly. The merged `dd/plots/dd_kin_C15d.root` drops `phi`, but the per-run
`~/C15d_fit/run_*_kin_d.root` trees keep it, so the azimuth is available on both sides
(Spyral: `azimuthal` in the InterpSolver parquets). 42,619 tracks match on (run, event), and
reproduce the earlier numbers exactly: `<dtheta> = +0.447`, rms 0.669.

### The azimuth offset is real, rigid, and EXPLAINED

```
phi_A - phi_S      mean  +9.55 deg   rms 4.53   circular concentration R = 0.997
phi_A + phi_S      mean  +0.27 deg   rms 104.6  R = 0.013
```

So the codes share the azimuth sense and zero, with ATTPCROOT rotated +9.55 deg. It is uniform
across the pad plane (+9.71, +9.51, +9.58, +9.45, +9.33, +9.73 in the six 60 deg sectors), so it is
a rigid rotation, not a distortion.

**It is not the pad map.** Dumping `AtTpcMap::CalcPadCenter` for all 10,240 pads and comparing
pad-for-pad against Spyral's `data/padxy.csv`:

```
best-fit rigid rotation (Spyral -> ATTPCROOT)   -0.0000 deg
residual after rotation  mean 1.003 mm   median 0.675   max 1.586     <- sub-pad
per-pad azimuth difference  mean 0.0000 deg  rms 0.350  R = 0.99998
radius ratio 1.00019
```

The two pad planes are the same object. ATTPCROOT's `ThetaPad: 110.9` in the par file is therefore
NOT applied to hit positions (and with `TiltAng: 0.0` the Lorentz correction it rotates is off
anyway).

**It is an arc-length effect.** In a solenoid the azimuth of a track's momentum rotates
continuously along the trajectory while the polar angle is conserved, so two codes quoting "the
direction at the vertex" at slightly different points on the arc differ in phi and NOT in theta.
That predicts the offset scales with curvature, i.e. with 1/p_T -- and it does:

```
KE [MeV]      n      <dphi>      <Brho>
  6-12     6,662    +11.433      0.611
 12-18     5,412    +10.170      0.790
 18-24     6,463     +9.279      0.941
 24-30     6,681     +8.853      1.064
 30-36     5,515     +8.789      1.176
 45-60     4,918     +9.083      1.470
      correlation( <dphi>, 1/p_T ) = +0.977
```

**Conclusion: the 9.55 deg is benign** -- a reference-point convention, not a misalignment. It does
not touch theta, ke or Ex. It WOULD matter to any analysis that uses phi (azimuthal asymmetries,
pad-level efficiency maps), so do not compare phi between the codes without removing it.

### The theta offset is not a tilt, not pad geometry, and not a drift scale

A tilt of the assumed beam axis by B toward azimuth phi0 gives `dtheta = B cos(phi - phi0)`.
Fitting that to dtheta(phi) in 12 bins of 30 deg:

```
flat offset    A = +0.436 +- 0.003 deg
sinusoid       B =  0.147 +- 0.004 deg  toward phi0 = 219 deg
chi2/ndf  flat only  783.9
chi2/ndf  A + Bcos   811.8      <- adding the tilt term makes the fit WORSE
```

Harmonic decomposition in 36 bins of 10 deg (amplitudes in deg, all +-0.004):

```
n=1 360deg 0.161 | n=2 180deg 0.083 | n=3 120deg 0.261 | n=4 90deg 0.109
n=5  72deg 0.219 | n=6  60deg 0.019 | n=7 51deg 0.062 | n=8 45deg 0.070
```

No harmonic dominates and **n=6 is the smallest of all** (0.019 +- 0.004) -- the hexagonal pad
plane leaves no signature, so this is not pad geometry. A drift-velocity scale error is excluded
independently: it would scale the vertex z linearly, and the codes agree on vz to 2 mm over the
whole 1000 mm chamber, bounding any scale error below 0.25 %, which is worth only 0.06 deg at
theta = 58.

### What the theta offset DOES depend on

Restricting to the clean core (azimuth obeying the rigid rotation to +-3 deg, 28,593 tracks,
`<dtheta>` +0.474, KE ratio 1.0045):

```
theta_S     n     dtheta          KE [MeV]     n     dtheta
40-45    2,949   +0.541            6-12     6,662   +0.497
45-50    2,487   +0.581           12-18     5,412   +0.505
50-55    2,598   +0.621           18-24     6,463   +0.489
55-60    2,238   +0.606           24-30     6,681   +0.468
60-65    1,073   +0.632           30-36     5,515   +0.432
65-70    1,594   +0.651           36-45     5,932   +0.383
70-75      830   +0.616           45-60     4,918   +0.363
```

It grows toward high theta and low KE -- the same trend seen twice, since the two are correlated
along the elastic locus. It is flat in fit quality (+0.474 at chi2/ndf < 0.2, +0.486 at 0.2-0.5).

**Leading remaining hypothesis: the energy-loss model in the vertex extrapolation.** At low energy
dE/dx is largest and the extrapolation from the first measured hit back to the beam axis is most
model-dependent; ATTPCROOT uses GENFIT with CATIMA, Spyral integrates its own eloss table. That is
consistent with every dependence above but is NOT yet measured -- comparing the two stopping-power
tables is the next step.

**Propagated budget.** At theta = 58, dT/dtheta = -1.26 MeV/deg and dT/dEbeam = 0.114, so
+0.45 deg is worth +5.0 MeV of apparent beam energy and the +0.45 % KE ratio another +1.2 MeV:
**~6 MeV of the 12-14 MeV gap**, about half. The rest is still unaccounted for.

⚠ One unit trap found while doing this: **Spyral's parquet stores the vertex in METRES**
(`vertex_z` 0-1), while `spyral_kin.root` holds mm. A vz binning in mm applied to the parquet
silently puts every track in the first bin.
