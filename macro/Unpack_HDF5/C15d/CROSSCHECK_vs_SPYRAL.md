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
