# UKF vs GENFIT on a1954 12Be — is the difference the covariance?

Test of the hypothesis that the UKF/GENFIT disagreement comes from the covariance
treatment. Data: the (p,p') proton-gated tracks, **all 73 clean runs** (~32 000 protons) for
the production and swapped-sigma configurations, refit from the *same* gated input so every
configuration fits identical tracks. The material-effects check is 8 runs (~3 300 protons).

## What the two fitters were actually assuming

| | measurement sigma | where it enters |
|---|---|---|
| UKF (`fitUKF_Be12.C`) | **0.5 mm** | `AtFitterUKF`: `chi2 += dist^2 / fMeasSigma_mm^2`, and the filter gain |
| GENFIT (`fitGenfit_Be12.C`) | **4.0 mm** | `AtGenfitter`: `s2 = fMeasSigmaMM^2` in the measurement covariance |

An 8x difference in assumed per-cluster resolution — 64x in chi2.

## Result 1 — the track-yield difference is ENTIRELY covariance. Hypothesis confirmed.

chi2/ndf, 73 runs, no cut:

| fitter | N | median chi2/ndf | q25 | q75 |
|---|---|---|---|---|
| UKF | 32 120 | **0.094** | 0.038 | 1.03 |
| GENFIT | 31 578 | **6.64** | 1.90 | 16.4 |

With **no** chi2 cut the two fitters return the same number of tracks (32 120 vs 31 578,
1.7 % apart) — GENFIT is *not* failing to converge. The shared `chi2/ndf < 5` cut keeps
83 % of UKF but only 40 % of GENFIT, and that single cut is the whole 26 604 vs 12 532 gap.

Swapping the sigmas moves the chi2 scale exactly as 1/sigma^2 — **all 73 runs**:

| configuration | N | median chi2/ndf | N(chi2<5) | kept |
|---|---|---|---|---|
| UKF sigma=0.5 (production) | 32 120 | 0.0945 | 26 604 | 82.8 % |
| UKF sigma=4.0 (swapped) | 32 173 | 0.0024 | 32 002 | 99.5 % |
| GENFIT sigma=4.0 (production) | 31 578 | 6.639 | 12 533 | 39.7 % |
| GENFIT sigma=0.5 (swapped) | 31 573 | 424.1 | 1 315 | 4.2 % |

GENFIT's chi2 ratio is **63.9** against the predicted (4.0/0.5)^2 = 64 — the measurement
covariance rescales chi2 and nothing else. UKF's is 40.1, not 64, because its measSigma
also enters the filter gain, so the estimate itself shifts slightly.

> **A `chi2/ndf` cut is not portable between fitters, or even between sigma settings of one
> fitter.** Set it per fitter from a percentile of that fitter's own chi2 distribution.
> At sigma=4.0 the cut at 5 throws away 59 % of perfectly good GENFIT tracks.

## Result 2 — the energy-scale bias is NOT covariance. Hypothesis refuted.

g.s. centroid, no chi2 cut. All 73 runs, except the matEffects row (8 runs):

| configuration | N | mu [MeV] | sigma [MeV] |
|---|---|---|---|
| UKF sigma=0.5 | 32 120 | +0.037 | 0.255 |
| UKF sigma=4.0 | 32 173 | +0.061 | 0.252 |
| GENFIT sigma=4.0 | 31 578 | **+0.451** | 0.262 |
| GENFIT sigma=0.5 | 31 573 | **+0.451** | 0.262 |
| GENFIT sigma=4.0, matEffects=ON (8 runs) | 3 385 | +0.637 | 0.757 |

Changing GENFIT's assumed covariance by 8x moves its g.s. centroid by **+0.000 MeV** — the
two agree to every digit fitted, on 31 500 tracks each. The +0.45 MeV offset therefore has
a different cause entirely. (UKF moves by +0.025 MeV over the same 8x change, again because
its measSigma feeds the gain as well as the chi2.)

`matEffects=kTRUE` does not fix it — it makes the centroid worse and the fits degenerate
(chi2 median 1e9 = ndf <= 0), so material effects are not usable in this configuration.

## Result 3 — the bias grows with lab angle (the lead for what it IS)

g.s. centroid per theta_lab bin, 73 runs, no chi2 cut:

| theta_lab | UKF mu | GENFIT mu | GENFIT - UKF |
|---|---|---|---|
| 60-70 deg | -0.020 | +0.374 | **+0.394** |
| 70-80 deg | +0.058 | +0.499 | **+0.441** |
| 80-90 deg | -0.205 | +2.550 | **+2.755** |

(The 20-60 deg bins hold no elastic peak — the recoil protons are at 60-90 deg — so no fit.)

The offset scales with angle, i.e. with path length in the gas. That points at energy-loss
handling / the reference point at which KE is reported (vertex vs first measured point),
not at the covariance. This matches the earlier UKF/GENFIT KE ratio observation
(~0.9 forward rising to 1.5-1.6 at 70-85 deg).

## Reproduce

```bash
# swapped sigmas (UKF at 4.0, GENFIT at 0.5) over the already-gated input
./sigma_test_Be12.sh "run_0143 run_0144 ..."          # -> ~/a1954_Be12_sigtest/
# caches with the chi2 cut removed, per fitter
root -b -q 'pp/ex_Be12.C("<runs>","<fitdir>",155,1e9,"_tag",1.007825,12.026921,"12Be(p,p&apos;)","<fitter>")'
```

---

# Part 2 — what the Ex bias actually IS (energy loss, not the reference point)

Follow-up to Result 2/3 above. Two candidates were tested: the point at which KE is
reported, and the energy-loss treatment inside the fit.

## How each fitter reports KE

| | KE is taken at | energy loss in the fit |
|---|---|---|
| UKF (`AtFitterUKF`) | smoothed state at the first cluster, then **back-extrapolated to the beam axis with the lost energy added back** (`KE_at_vertex = KE_at_cluster + dEdx*pathLength`, CATIMA) | **yes** — the propagator carries an `AtELossModel` |
| GENFIT (`AtGenfitter`) | `gfTrack->getFittedState()` = the state at the **first track point**, no vertex extrapolation | **no** — production runs `matEffects=kFALSE` |

UKF conveniently stores both: `GetKinematics()` = vertex, `GetKinematicsXtr()` = first
cluster. That separates the two effects with no refitting.

## The reference point is a SMALL effect (13 % of the bias)

g.s. centroid, 73 runs, no chi2 cut:

| | mu [MeV] |
|---|---|
| UKF, KE at vertex (production) | +0.037 |
| UKF, KE at first cluster (vertex correction removed) | +0.092 |
| GENFIT | +0.451 |

Removing UKF's vertex correction moves the centroid by only **+0.055 MeV** — the energy UKF
adds back is a median of **0.082 MeV** per track. GENFIT still sits **+0.359 MeV** away from
UKF *at the same reference point*, so ~87 % of the bias is not the reference point.

## It is the missing energy loss in the track model

Per-track matching (same run+event, nearest theta within 3 deg): **25 397 matched pairs**.
Both sides taken at the first cluster, so the reference point is factored out:

| UKF KE bin [MeV] | n | median KE(UKF) − KE(GENFIT) | as a fraction of KE | median track length [mm] |
|---|---|---|---|---|
| 0–2 | 6 131 | +0.491 | **31 %** | 322 |
| 2–4 | 9 701 | +0.970 | **35 %** | 774 |
| 4–6 | 3 519 | +1.416 | **30 %** | 1 609 |
| 6–10 | 2 037 | +0.678 | 9.5 % | 1 403 |
| 10–15 | 807 | +0.159 | 1.3 % | 669 |
| 15–25 | 1 986 | −0.333 | −1.7 % | 403 |
| 25–40 | 1 179 | −0.479 | −1.7 % | 449 |

GENFIT is ~**30 % low below 6 MeV**, converges to UKF by ~10–15 MeV, and is marginally high
above that. That is the signature of fitting a **constant-momentum helix to a decelerating
track**: with no energy loss in the model the curvature it returns is a track average, which
sits well below the initial momentum when the particle loses a large fraction of its energy
over the measured length. It explains everything seen so far:

* the Ex bias grows toward 80–90 deg — those are the lowest-energy recoil protons;
* (p,d) deuterons (forward, ~10–30 MeV) show UKF and GENFIT agreeing;
* the bias is completely independent of the assumed covariance.

## Why turning material effects on does not rescue it

`matEffects=kTRUE` on 8 runs: only **468 / 3 385 tracks (13.8 %)** come back with ndf > 0,
and even that good subset keeps the bias (mu = +0.469 vs +0.438 with effects off).
Two reasons it cannot work as configured:

1. the fit is unstable with material effects on (86 % of tracks end with ndf <= 0, and
   `AtGenfitter` silently retries them with effects off);
2. the geometry genfit would use is `geometry/ATTPC_H1bar_geomanager.root` — **H2 at 1 bar,
   while a1954 ran at 600 torr**, i.e. ~27 % too dense. Even a working correction would be
   miscalibrated.

## Recommendation

**Use UKF for this data set** — it is the only one of the two that models dE/dx during the
fit and reports KE at the vertex. GENFIT is usable above ~10–15 MeV (where the two agree to
~1 %), which is why it is fine for the forward (p,d) deuterons and not for the backward
(p,p') recoils. Making GENFIT correct at low energy needs a 600 torr H2 geometry *and*
whatever fixes the ndf <= 0 instability — until then its low-energy KE is not trustworthy.

---

# Part 3 — THE FIX. GENFIT now agrees with UKF.

Both problems identified in Part 2 turned out to have a single root cause plus a wrong gas.

## Root cause of the "material effects are unusable" instability

`AtGenfitter::Init()` called `genfit::MaterialEffects::useEnergyLossParam()`
**unconditionally**. That switches genfit to parameterization-ONLY energy loss, whose
validity window is `maxKinEnergy_` — a member that starts at **0** and is set *only* by
`setEnergyLossFile()`. Our macros never pass an eloss table, so `maxKinEnergy_` stayed 0 and

```
if (E_kin > maxKinEnergy_) { Exception exc("MaterialEffects::dEdx ==> kinetic energy out of
                             parameterization boundaries!"); exc.setFatal(); throw exc; }
```

fired on **every single step**. genfit's `KalmanFitterRefTrack` swallows that while building
the reference track and simply truncates it, so the fit came back "converged" with ndf = −3
(two measurements minus five parameters) instead of throwing — which is why `AtGenfitter`'s
exception fallback never triggered and 87 % of tracks silently produced garbage.

That also explains why the failure count was **exactly 81/603 surviving tracks regardless of
gas density, drift-volume position, or seed direction**: it never depended on the geometry.

**Fix** (`AtGenfitter.cxx`): only request the parameterization when a table is actually
loaded. Otherwise genfit stays in Param+Bethe-Bloch mode, valid down to beta*gamma = 0.05
(~1.2 MeV protons).

```cpp
if (!fELossFile.empty()) {
   mat->setEnergyLossFile(fELossFile, fPDG);
   mat->useEnergyLossParam();
}
```

## Plus the right gas

New `geometry/ATTPC_H600torr.C` + an `H_600torr` entry in `media.geo` at
**6.616e-5 g/cm^3** (H2, 600 torr, 293 K — the H_1bar value 8.27e-5 reproduces PM/RT at
1 bar exactly, so this is the same number scaled by 79993/100000 Pa). `fitGenfit_Be12.C`
takes a `geoName` argument, now defaulting to `ATTPC_H600torr`. Build it once with
`root -b -q geometry/ATTPC_H600torr.C`.

## Result — all 73 runs, no chi2 cut

| configuration | N | ndf > 0 | mu [MeV] | sigma [MeV] |
|---|---|---|---|---|
| UKF (production) | 32 120 | 100 % | +0.037 | 0.255 |
| GENFIT old (matEff OFF, 1 bar) | 31 578 | 99.5 % | **+0.451** | 0.262 |
| **GENFIT fixed (matEff ON, 600 torr)** | 32 179 | **99.9 %** | **+0.050** | 0.250 |

The +0.45 MeV bias is gone — GENFIT and UKF now agree to **13 keV** on the g.s. centroid,
with GENFIT marginally the narrower of the two.

Per-track, matched to UKF (28 734 pairs), the energy deficit that Part 2 measured:

| UKF KE bin | old deficit | fixed deficit |
|---|---|---|
| 0–2 MeV | +30.9 % | **+15.9 %** |
| 2–4 | +34.8 % | **−2.7 %** |
| 4–6 | +29.6 % | **−1.0 %** |
| 6–10 | +9.6 % | **+0.3 %** |
| 10–15 | +1.3 % | −0.3 % |
| 15–25 | −1.7 % | −2.1 % |
| 25–40 | −1.7 % | −1.9 % |

Everything above 2 MeV now agrees with UKF to ~1–3 %. The residual **+16 % below 2 MeV is
expected and is the known limit of the fallback**: under beta*gamma = 0.05 (~1.2 MeV for
protons) genfit leaves Bethe-Bloch and calls `dEdxParam`, which returns 0 when no table is
loaded — i.e. no energy loss is modelled there. Loading a proton dE/dx table via
`SetELossFile()` would close that last gap.

## Conclusion

The two fitters were never really disagreeing about physics. The gap was
(1) a chi2 scale set by two different assumed sigmas, and (2) genfit running with its
energy loss effectively disabled by a misconfigured parameterization mode. With both fixed,
UKF and GENFIT agree on the g.s. centroid to 13 keV and on per-track KE to a few percent
above 2 MeV.
