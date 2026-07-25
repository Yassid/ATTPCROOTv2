# UKF vs GENFIT on a1954 12Be — is the difference the covariance?

Test of the hypothesis that the UKF/GENFIT disagreement comes from the covariance
treatment. Data: the (p,p') proton-gated tracks, 73 clean runs; the swapped-sigma and
material-effects refits use the first 8 runs (~3 300 protons) from the *same* gated input,
so every configuration fits identical tracks.

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

Swapping the sigmas moves the chi2 scale exactly as 1/sigma^2:

| configuration | N | median chi2/ndf | N(chi2<5) |
|---|---|---|---|
| GENFIT sigma=4.0 | 3 297 | 6.48 | 1 340 (41 %) |
| GENFIT sigma=0.5 | 3 302 | 415.15 | 159 (5 %) |
| UKF sigma=0.5 | 3 367 | 0.078 | 2 813 (84 %) |
| UKF sigma=4.0 | 3 372 | 0.0020 | 3 349 (99 %) |

GENFIT's ratio is **64.0**, i.e. exactly (4.0/0.5)^2. UKF's is 38, not 64, because its
measSigma also enters the filter gain, so the estimate itself shifts slightly.

> **A `chi2/ndf` cut is not portable between fitters, or even between sigma settings of one
> fitter.** Set it per fitter from a percentile of that fitter's own chi2 distribution.
> At sigma=4.0 the cut at 5 throws away 59 % of perfectly good GENFIT tracks.

## Result 2 — the energy-scale bias is NOT covariance. Hypothesis refuted.

g.s. centroid, same 8 runs, no chi2 cut:

| configuration | N | mu [MeV] | sigma [MeV] |
|---|---|---|---|
| UKF sigma=0.5 | 3 367 | +0.025 | 0.248 |
| UKF sigma=4.0 | 3 372 | +0.043 | 0.234 |
| GENFIT sigma=4.0 | 3 297 | **+0.438** | 0.257 |
| GENFIT sigma=0.5 | 3 302 | **+0.437** | 0.258 |
| GENFIT sigma=4.0, matEffects=ON | 3 385 | +0.637 | 0.757 |

GENFIT's centroid is **identical to 1 keV** at sigma = 0.5 and 4.0. An 8x change in assumed
covariance does not move it at all, so the +0.4 MeV offset has a different cause.

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
