# 15C(d,d'): fitting the states seen in 15C(p,p')

Done 2026-11-09 on the user's saved explorer parameters
(`~/Downloads/explorer_params_15C_d_d__2026-11-09-14-30-49_{atttpcroot,spyral}.json`).
Figures: **`figures/`** (`dd_state_fit.png` is the main one; see `figures/README.md`).
Scripts: `~/Data/C15d_statefit/` (`exlib.h`, `ddbuild.C`, `statefit2.C`, `syst.C`).

## The reference: what 15C(p,p') shows

From `C15p/pp/plots_pgate/pp_kin_C15p.root` (proton-GATED, 131,532 tracks after the a2091 proton
IC window, chi2 < 5) at Ebeam 198:

```
g.s.  mu +0.119   FWHM 0.511 MeV      <- three times better resolution than (d,d')
inelastic structures, >5 sigma over a fitted background:
   4.75 (19 sigma)   7.05 (53)   7.75 (14)   10.35 (46)   11.05 (49)
```

⚠ The 10.35 and 11.05 structures sit on the falling edge where the spectrum ends abruptly near
12 MeV, so an acceptance origin is not excluded for them. 4.75 and 7.05 are well clear of it.

## The user's parameters reproduce offline

`exlib.h` re-implements the page's own transforms exactly (page lines 408-575):

```
theta_corr = theta - (360/kcDenom)*(KE + keOff - kcPivot)        [deg]
Ebeam(vz)  = ebeam - (ebzSlope/1000)*(vz - ebzPivot)             [MeV, vz in mm]
Ex         from the same kine2b expression, U = 931.49401
```

With the saved parameters this gives

```
              tracks    g.s. centroid   FWHM
ATTPCROOT    129,450       +0.053       1.012      Ebeam 191.25, chi2<0.2, IC [1074,1201],
                                                   kcDenom 3913, kcPivot -1.5, ebz 11 MeV/m
Spyral        75,249       -0.253       0.548      Ebeam 195, chi2<5, no IC, vz<400,
                                                   kcDenom 2535.2, kcPivot 14.5, ebz 11 MeV/m
```

**That settles the beam-energy question left open earlier.** A naive scan without the kc and ebz
corrections suggested the ATTPCROOT g.s. could only be zeroed at 212 MeV and only by paying
0.38 MeV of width. With the corrections, **191.25 puts it at +0.053 with FWHM 1.012** -- better
than anything that scan reached. The earlier ambiguity was an artefact of fitting without them.

## The fit

Smooth continuum (`exp(a+bx) + c + dx`) plus two Gaussians, over 2.3-13.5 MeV.

```
                background only   +2 peaks   delta-chi2 (6 dof)
ATTPCROOT     chi2/ndf  14.18       1.78         1350
Spyral        chi2/ndf   8.45       1.80          729
```

With the positions left FREE, both datasets land on the same two structures:

```
              state A                    state B                 yields
ATTPCROOT   4.657 +- 0.038 +- 0.031   6.879 +- 0.022 +- 0.005   3260 +- 128 / 2394 +- 97
Spyral      4.500 +- 0.042 +- 0.019   6.643 +- 0.026 +- 0.005   1390 +-  85 / 1532 +- 86
(p,p')           4.75                      7.05
```

(second error = spread over 12 fit ranges, 2.0-3.0 low edge x 12.0-15.0 high edge.)

**Only these two of the five (p,p') structures are present.** Fixing all five and freeing the
amplitudes returns yields consistent with zero for 7.75, 10.35 and 11.05 in both datasets --
consistent with those being acceptance features in (p,p') rather than states, or simply not
populated in (d,d').

## ★ The cross-check that matters

**ATTPCROOT's state A, 4.657 +- 0.038, lands exactly on the 4.657 MeV state the a1975
16C(d,t)15C analysis extracts** (`macro/Simulation/ATTPC/16C_dt/run_states_dt.sh`: gs, 0.740,
3.103, **4.657**, 6.358). Different reaction, different dataset, different experiment.

**And Spyral's state B, 6.643 +- 0.026, matches the 6.648 that the same analysis records as
"the feature the data puts at 6.648"** against a literature 6.358.

## The ATTPCROOT - Spyral offset is the known fitter systematic

```
state A  4.657 vs 4.500   ->  +0.157
state B  6.879 vs 6.643   ->  +0.236
```

Same sign and magnitude as the +0.45 deg polar-angle offset documented in
`CROSSCHECK_vs_SPYRAL.md`, which was localised to the FITTERS (the hit clouds and geometry agree
to +-0.25 deg). So the two codes disagree on these level energies by ~0.2 MeV for the same reason
they disagree on the beam energy -- **do not average them, and do not quote one with the other's
Ebeam.**

## Widths: observed, not deconvolved

```
             g.s. FWHM    state A FWHM    state B FWHM
ATTPCROOT      1.012         1.898           1.117
Spyral         0.548         1.488           1.169
```

State A is clearly broad in both. State B is close to the ATTPCROOT resolution but not to Spyral's.
⚠ **Deconvolving a natural width is NOT justified here**: the resolution is measured only at
Ex = 0, and in both codes it degrades with theta (and therefore with Ex). Subtracting the g.s.
width in quadrature gives 1.61/0.47 (ATTPCROOT) against 1.38/1.03 (Spyral) -- the state-B numbers
disagree by 2x, which is the signature of the assumption failing, not of a measurement.
A natural width needs the resolution measured AT that excitation, e.g. from simulation.

---

## The ground state and the 5/2+ at 0.740

The 5/2+ IS present in both datasets. It appears as a SHOULDER on the high-Ex side of the elastic
peak, not a resolved second peak (FWHM 0.92 / 0.52 against a 0.740 spacing).

```
                 separation FREE     asym 1-peak -> +companion   dchi2 (1 dof)
ATTPCROOT all     0.762 +- 0.019         5.02 -> 3.30               283
ATTPCROOT 65-85   0.732 +- 0.016         4.67 -> 3.04               267
Spyral all        0.727 +- 0.013         6.43 -> 1.42               807
Spyral 65-85      0.689 +- 0.012         5.93 -> 1.32               738
```

Four independent measurements give **0.73 +- 0.04** against the known **0.740**.

⚠ **The yield ratio is NOT trustworthy**: 0.106 +- 0.005 (ATTPCROOT) against 0.041 +- 0.001
(Spyral), 2.6x apart for the same physical quantity. The shoulder's amplitude depends on how the
main peak's tail is modelled, and the broader ATTPCROOT peak assigns more of it to the companion.
A real resolution function from simulation is needed before quoting a ratio.

### ⚠ Method: a second Gaussian always improves a non-Gaussian peak

The first pass used a symmetric main peak and gave dchi2 = 750-1100 for 1 dof — overwhelming and
nearly meaningless. Two checks are what make a shoulder credible, and both should be routine:
**float the separation** (a real state lands at its known energy; an artefact lands wherever the
tail is) and **let the main peak be asymmetric** (if the companion only absorbed a tail it loses
significance). A third warning sign prompted them: the amplitude ratio jumped 5x between angular
selections, which a real state cannot do.

---

# The continuum under the states: phase space (2026-11-09 overnight)

Figure `~/C15d_continuum_phasespace.png`. Scripts `scripts/psmodel.h`, `psmc.C`, `psalt.C`,
`tmplfit.C`, `effcheck.C`.

The states sit on a large continuum -- **12.1 % of the gated sample between Ex 8 and 22 MeV for
ATTPCROOT, 17.0 % for Spyral** -- so what shape that continuum has is a systematic on every state
energy. Two three-body hypotheses were simulated with `TGenPhaseSpace`, each reconstructed exactly
as the analysis does (two-body missing mass, deuteron mass hypothesis, same KE and theta cuts):

```
(A)  15C + d -> d + 14C + n     target excitation above Sn, the DEUTERON detected
(B)  15C + d -> p + n + 15C     deuteron BREAKUP, the PROTON detected and mis-assigned a
                                deuteron mass -- what a breakup proton leaking the gate looks like
```

Sanity check on the model: Sn(15C) computed from the masses is **1.2180 MeV** against the
literature 1.218, and the kinematic limit at Ebeam 191.25 is Ex = 22.61 MeV, where the data ends.

## (A) is excluded by its shape

Uniform d + 14C + n phase space RISES from its 2 MeV threshold to a maximum near 14-16 MeV. The
data FALLS from 6 MeV onward. The ratio data/(A) drops 6x across the range.

**That cannot be reconstruction efficiency**: within the acceptance the mean deuteron energy
*increases* with Ex (7.1 MeV at Ex = 2 to 13.5 at Ex = 20), so the tracks are getting FASTER where
the data is losing counts. A short-track efficiency loss would do the opposite.

## (B) matches the shape, but NOT the yield

The breakup template is nearly flat and follows the data from 8 to 22 MeV including the fall-off
at the kinematic limit. As a fit template it beats (A) in both datasets:

```
                      chi2/ndf     state A           state B
ATTPCROOT  tmpl B       2.23    4.559 +- 0.040    6.933 +- 0.023
           tmpl A       3.39    4.705 +- 0.044    8.000 +- 0.098   <- railed at the bound
           empirical    1.66    4.679 +- 0.037    6.874 +- 0.022
Spyral     tmpl B       1.87    4.459 +- 0.043    6.654 +- 0.026
           tmpl A       2.10    4.796 +- 0.073    6.690 +- 0.028
           empirical    1.59    4.501 +- 0.042    6.644 +- 0.026
```

⚠ **But the yield rules it out as stated.** The continuum is 12-17 % of the gated sample while the
gate contains **1.0 % protons** (+2.2 % unassigned), measured against the species rails. And a
breakup proton cannot hide on the deuteron rail: the rails are a dE/dx separation ~1.78x apart,
set by the particle, not by the reaction. A flat template matching a flat continuum is weak
evidence on its own.

## What the continuum most likely is

Genuine DEUTERONS exciting the 15C continuum above Sn non-resonantly. That keeps the particle
identity consistent with the gate, and explains why neither uniform phase space describes it: the
shape of direct inelastic excitation into the continuum is set by the reaction dynamics, not by a
uniform three-body final state. **Not established here** -- it would need a reaction model.

## The useful outcome: a background systematic on the states

The three background models give

```
             state A span            state B span
ATTPCROOT    4.559 - 4.705  (0.15)   6.874 - 6.933  (0.06, excluding the railed A fit)
Spyral       4.459 - 4.796  (0.34)   6.644 - 6.690  (0.05)
```

so **state B is robust against the background model and state A is not.** Quote state A with a
background systematic of ~0.1 (ATTPCROOT) to 0.2 (Spyral) MeV, on top of the fit-range systematic
already in the table above.

---

# The complete spectrum in one fit (2026-11-10)

`scripts/globalfit.C` + `scripts/plotglobal.C` -> `figures/global_ex_fit.png`

Everything above fits pieces of the spectrum separately. This fits **all of it at once**: the
15C ground state, the 0.740 MeV 5/2+, state A, state B and the background, for both codes.

![global fit](figures/global_ex_fit.png)

```
                 g.s. mu          FWHM(L/R)     5/2+ amp/g.s.         state A              state B
ATTPCROOT   +0.005 +- 0.005   0.906 / 0.906   0.100 +- 0.006     4.681 +- 0.037 +- 0.065   6.871 +- 0.022 +- 0.021
                                                  +- 0.012(bkg)     FWHM 1.78                 FWHM 1.10
Spyral      -0.287 +- 0.003   0.479 / 0.567   0.033 +- 0.001     4.510 +- 0.043 +- 0.009   6.643 +- 0.027 +- 0.004
                                                  +- 0.009(bkg)     FWHM 1.47                 FWHM 1.15
```
(first error statistical, second the background systematic defined below; MeV throughout.
Background adopted: the **empirical straight line for BOTH codes** -- see the warning below.)

## Two windows, not one minimisation -- and why

The g.s. holds ~9000 counts in its peak bin, state A ~300 in total. In a **single** Poisson
likelihood over the whole range the g.s. lineshape -- which is not Gaussian: sharp core, longer
right tail -- dominates the gradient and drags the weak states. All three of these were tried and
all three failed, on this data:

| attempt | what happened |
|---|---|
| one fit, free widths | A collapsed to the width floor, B blew up to FWHM 3.06 (ATTPCROOT); both pinned at 0.283 (Spyral) |
| one fit, widths bounded <= 1.8 FWHM | the background exponential slope railed at its bound for **both** codes |
| `erf x (exp + const)` background | Spyral's A railed at 5.800, B at 8.000, background scale at its floor -- background and states swapped roles |

The peak region and the state region are physically disjoint; the only thing linking them is the
5/2+ on the g.s. right shoulder. So they are fitted in **two windows sharing one background**
(split at 2.2 MeV), iterated to self-consistency -- the state window's background frozen into the
peak window, the peak window's g.s. tail frozen into the state window, four passes -- and drawn as
one curve. The per-region chi2/bin shows where the model is and is not adequate:

```
                g.s. core   1.5-3.5   state A   state B   bkg 8-20
ATTPCROOT          9.33       5.25      1.27      1.13      1.39
Spyral             2.91       2.51      1.15      1.08      1.58
```

**The states sit at chi2/bin ~1.1-1.3.** The large numbers are confined to the g.s. core, where a
bifurcated Gaussian is simply not the true response, and to its right shoulder. That mismatch is a
lineshape limitation, not a fit failure, and it does not propagate into A and B (next section).

## The background is a systematic, not a choice

Three forms were fitted independently, none of which can grow a Gaussian bump:

* `LIN` &nbsp; `c0 + c1*Ex` -- empirical, monotonic, 2 parameters
* `PS` &nbsp;&nbsp; `c + N*dens(Ex)*exp(-k*Ex)` -- flat floor + the three-body `15C+d -> d'+(14C+n)`
  phase space of `psmodel.h`, threshold fixed at Sn, one acceptance factor, 3 parameters
* `THR` &nbsp;`c + A*Phi((Ex-Sn)/sigma_res)*exp(-k*Ex)` -- flat floor + a continuum turning on at Sn
  smeared by that code's **own measured** g.s. sigma_R (not refitted), 3 parameters
* `FLAT` &nbsp;`c` -- a constant, 1 parameter. No threshold, no phase space.
* `FLATEXP` &nbsp;`c + A*exp(-k*Ex)` -- 3 parameters. No threshold, no phase space.

```
bkg        -2lnL   chi2/ndf   g.s. mu    5/2+ ratio      state A          state B
                        -- ATTPCROOT --
LIN        281.2     3.05     +0.0052   0.0997(58)   4.681(37) FWHM 1.78   6.871(22) FWHM 1.10
PSATT      204.6     4.35     +0.0340   0.1228(45)   4.626(39) FWHM 1.97   6.878(21) FWHM 1.00
THR        305.0     7.83     +0.0806   0.1221(39)   4.685(37) FWHM 1.76   6.871(22) FWHM 1.11
FLAT       602.4     5.14     +0.0156   0.1118(49)   4.555(41) FWHM 2.19   6.913(26) FWHM 1.40
FLATEXP    306.4     3.22     +0.0038   0.0984(59)   4.684(37) FWHM 1.76   6.871(22) FWHM 1.12
                        --  Spyral   --
LIN        272.4     2.68     -0.2872   0.0326(14)   4.510(43) FWHM 1.47   6.643(27) FWHM 1.15
PSATT      255.7     1.86     -0.2854   0.0356(14)   4.500(42) FWHM 1.52   6.640(25) FWHM 1.08
THR        272.6     3.15     -0.2662   0.0513(14)   4.510(43) FWHM 1.47   6.643(27) FWHM 1.15
FLAT       279.8     2.31     -0.2865   0.0338(14)   4.492(42) FWHM 1.55   6.649(27) FWHM 1.20
FLATEXP    272.6     2.68     -0.2872   0.0326(14)   4.510(43) FWHM 1.47   6.643(27) FWHM 1.15
```

`-2lnL` is the state-window likelihood, which is what the fit actually minimises; lower is better.

**Spread across all five background forms:**

```
             state A    state B    5/2+ ratio
ATTPCROOT    0.130      0.043       22 %
Spyral       0.018      0.009       45 %
```

So **A and B are robust and the 5/2+ yield is not.** This still improves on the earlier background
systematic for state A (0.15-0.34 MeV, from fits to the state region alone) -- constraining the
background over the full 20 MeV instead of a narrow window is what buys that -- though less
dramatically than the three-form version suggested, because `FLAT` is an outlier for ATTPCROOT.

**Adopted background: `LIN` for BOTH codes.** It is empirical and makes no physical claim. `PSATT`
fits better in both codes but is *not* adopted, for the reason in the next section.

⚠ **Two statistics traps, both of which produced a wrong answer here first:**

1. **Do not compare these models by a Pearson chi2 restricted to bins with n >= 10.** The fits
   minimise a Poisson likelihood over *all* bins, so the two are inconsistent, and `FLATEXP` came
   out WORSE than the `FLAT` model nested inside it (chi2/ndf 2.68 vs 2.31) -- impossible. Compare
   by `-2lnL`.
2. **Do not filter the systematic by a likelihood-ratio test.** With 129k events every background
   form differs from every other at far beyond 5 sigma, so an LR cut keeps exactly one model and
   reports a systematic of **zero**. None of these shapes is the truth; the spread between
   reasonable ones *is* the uncertainty.

A third: a single starting point lands in a bad local minimum often enough to invert a nested
comparison, so `globalfit.C` now runs 12 starts per model and keeps the best likelihood.

A real consistency check on the floor: the fitted flat floor is compared against each code's own
measured counts/bin **below** threshold, which nothing in the fit is told about. The two codes
differ by an order of magnitude there (ATTPCROOT ~70/bin at Ex = -1.9, Spyral ~7/bin, because
Spyral cuts vz < 400 mm and vertexes better), so it is not a tautology. `PS` on ATTPCROOT fits
66.7 +- 4.5 against a measured 69.5 -- it passes. `THR` drives the floor to zero in both codes and
has the worst chi2 for ATTPCROOT, which is why it is not adopted.

## ⚠⚠ Why `PSATT` wins the fit and still means nothing

`PSATT` has the best likelihood in both codes -- and at *equal parameter count* it beats a plain
exponential decisively:

```
                FLAT     LIN     FLATEXP   PSATT    THR      (state-window -2lnL, lower is better)
ATTPCROOT      602.4    281.2     306.4    204.6    305.0
Spyral         279.8    272.4     272.6    255.7    272.6
```

`PSATT` vs `FLATEXP` is 3 parameters against 3, so the phase-space *shape* is worth
Delta(-2lnL) = 102 (ATTPCROOT) and 17 (Spyral). It is tempting to read that as evidence that the
continuum **is** three-body 14C+n phase space. **It is not, and it directly contradicts the
2026-11-09 result above only if you ignore where the win comes from.**

The model is `c + N*dens(Ex)*exp(-k*Ex)`. The fitted attenuation is

```
ATTPCROOT  k = 0.0909 /MeV  ->  Ex = 20 suppressed 5.1x relative to Ex = 2
Spyral     k = 0.0662 /MeV  ->                     3.3x
```

That **is** the discrepancy the earlier study measured: *"Uniform d + 14C + n phase space RISES
from its 2 MeV threshold to a maximum near 14-16 MeV. The data FALLS from 6 MeV onward. The ratio
data/(A) drops 6x across the range."* And that study showed the 6x **cannot be reconstruction
efficiency**, because within the acceptance the mean deuteron energy *rises* with Ex (7.1 MeV at
Ex = 2 to 13.5 at Ex = 20) -- the tracks get faster exactly where the counts are lost.

So `PSATT` does not describe the continuum with phase space. It takes a shape that is known to be
wrong by ~6x and multiplies it by a free parameter that absorbs exactly that factor. What it
actually buys over `FLATEXP` is a *threshold at Sn with an extended sqrt-like rise* rather than no
threshold at all -- a real shape statement about the data, but not a phase-space interpretation.

**Conclusion: the 2026-11-09 exclusion of uniform d + 14C + n phase space STANDS.** `PSATT` is
retained only as a shape systematic, and the adopted background is the empirical straight line.
Had I adopted `PSATT` on its likelihood -- which is what I did first, labelling it "floor + 3-body
phase space" -- the figure would have asserted a continuum mechanism this analysis has excluded.

For the record, how much of the background the phase-space term even carries when it is used:

```
                floor   phase-space term   PS share of the background (Sn to 20 MeV)
ATTPCROOT        66.7      N = 0.0077                50 %
Spyral           78.0      N = 0.0015                17 %
```

For Spyral the "phase space" background is 83 % flat floor.

## ★ Referencing to each code's own g.s. removes most of the inter-code difference

```
                    state A    state B
ATTPCROOT absolute   4.681      6.871
Spyral    absolute   4.500      6.640
difference          +0.181     +0.231

ATTPCROOT rel. g.s.  4.676      6.866
Spyral    rel. g.s.  4.786      6.926
difference          -0.110     -0.060
```

Spyral's g.s. sits at -0.286 and ATTPCROOT's at +0.005, and that single offset accounts for most of
the disagreement: quoted **relative to each code's own ground state**, A agrees to 0.11 MeV and B to
0.06 MeV. This is the expected behaviour if the difference is the known fitter/beam-energy scale
(see `c15d-attpcroot-vs-spyral`): a scale error moves the whole spectrum, so it cancels in a
difference. **Quote level energies relative to the g.s. centroid, not on the absolute Ex axis.**

## The 5/2+: visible, but its yield is not measurable here

The bottom row of the figure plots **data minus (fitted g.s. + background)**. The 0.740 MeV state
then appears as the bump it is, in both codes, with the 5/2+ Gaussian from the full fit drawn over
it (not refitted to the residual).

⚠ **The residual is not independent evidence.** The g.s. and the 5/2+ were fitted together, so
subtracting the g.s. necessarily leaves something at 0.740. What the panel shows is that the
leftover is a clean, centred, roughly Gaussian bump rather than a lineshape artefact -- and that it
survives in Spyral, where the g.s. FWHM is 0.48 and the two states are 1.5 FWHM apart.

The amplitude ratio differs by 3x between the codes (0.100 vs 0.035) and carries a 21-45 %
background systematic, so **it is not a yield measurement.** What limits it is quantified in the
next section -- and it is *not* an internal degeneracy of the fit, which was my first guess and is
wrong.

# What limits the 5/2+, and what does NOT (2026-11-10)

Two hypotheses were tested directly. The first failed, the second is the answer.

## ✗ REFUTED: the non-Gaussian g.s. is NOT a superposition over theta slices

`scripts/thetaslice2.C`. Resolution is strongly theta-dependent in both codes, so the obvious guess
is that the integrated g.s. is a sum of Gaussians of different width and is non-Gaussian only for
that reason. If so, fitting *within* a narrow theta slice should give chi2/bin ~ 1. It does not:

```
theta slice   64-68  68-72  72-76  76-80  80-84  84-88   |  ALL
ATTPCROOT      0.79   2.67   3.48   3.59  18.97   0.96   |  8.06
Spyral         0.31   1.34   2.45   3.74   4.29   4.14   |  1.47
```

The slices carrying the statistics are **no better than the integral and often worse** (ATTPCROOT's
dominant 80-84 slice is 18.97 against 8.06 for all angles). **The non-Gaussian response is
intrinsic within a 4 degree slice**, not an artefact of summing. So a per-slice response model
would not fix it, and that avenue is closed.

⚠ Method note: a first pass fitted each slice with a plain `gaus+pol1` and got chi2/bin 10-49,
which was meaningless -- the 5/2+ sits *inside* the fit window, so omitting it inflates chi2
precisely where the statistics are. The table above uses the same bifurcated-Gaussian + 5/2+ +
linear model as the global fit, so the comparison is like for like.

Incidentally the forward slices (40-64 deg) peak at Ex = 4.6 to 11.6 MeV, not at the g.s. -- the
elastic peak is simply not populated there. **69 % of all g.s. counts sit at theta = 76-88 deg.**

## ★ The 5/2+ ratio is limited by the g.s. LINESHAPE, not by statistics

`scripts/degen.C`, `scripts/degenfig.C` -> `figures/fivehalf_lineshape.png`

![5/2+ lineshape sensitivity](figures/fivehalf_lineshape.png)

Fix the g.s. right-hand width to a grid of values, refit everything else, and watch the 5/2+ ratio
and the likelihood:

```
                    best FWHM_R   Delta(-2lnL)=1 band   ratio there   slope       10 % width error
ATTPCROOT              0.906      0.894 - 0.917 (+-11 keV)  0.0997    -0.327/MeV  -> 30 % on the ratio
Spyral                 0.565      0.570 - 0.576 (+- 3 keV)  0.0353    -0.099/MeV  -> 16 % on the ratio
```

**The likelihood is steep, so the width is pinned to a few keV and the ratio is NOT internally
degenerate** -- my first explanation (that the two peaks overlap too much to separate) is refuted by
this scan. Within its own model each code determines the ratio sharply.

The problem is the *model*. The ratio slides fast with the assumed width, and the two codes are
bridged by a change of only **+21 %** in the right-hand width: put ATTPCROOT at FWHM_R = 1.09
instead of 0.906 and it reproduces Spyral's 0.0353 exactly. Given that the bifurcated Gaussian
leaves chi2/bin 9.3 in the ATTPCROOT g.s. core, a 20 % error in the right-hand width is entirely
plausible.

⇒ **The quoted statistical error on the ratio (0.0058 / 0.0014) understates the truth by an order
of magnitude.** The 5/2+ yield is a lineshape-limited quantity and cannot be improved by more
statistics -- only by a correct response function. Trust Spyral's number over ATTPCROOT's (its
peaks are 1.54 FWHM apart against 0.82, and its core chi2/bin is 2.9 against 9.3), and neither
better than a factor two.

## What is still missing

* A **simulated resolution function** at Ex ~ 0, 0.74, 4.6, 6.9. Without it the quoted FWHMs are
  observed widths, not natural widths, and the 5/2+ ratio stays a factor-two quantity. The scan
  above makes this concrete: the response has to be right to ~5 % in the right-hand width before
  the 5/2+ ratio is good to 15 %.
* The g.s. **lineshape** itself: a bifurcated Gaussian leaves chi2/bin 9.3 in the ATTPCROOT core
  and 2.9 in Spyral's. This is the single largest modelling deficiency left. It is NOT a theta
  superposition (refuted above), so it must come from the fitter's own KE/theta response -- which
  means a simulation is the only way to get at it.
* State A and state B are unaffected by all of this: chi2/bin 1.1-1.3 in their regions, and
  <= 0.06 MeV of movement across three background forms.

---

# Why we do not see the phase space (2026-11-11)

`scripts/whynops.C`, `scripts/vzcut.C`. Model-free -- no fit in the way. Counts per 100 keV bin,
both curves normalised to the data at Ex 8-10:

```
Ex window     ATTPCROOT   Spyral  | uniform PS | PS x halo form factor
  2.0-3.4        158        92    |     82     |       969
  8.0-10.0       160       108    |    160     |       160
 10.0-12.0       138        95    |    166     |       110
 14.0-16.0       113        87    |    159     |        57
 18.0-20.0        93        92    |    125     |        28
```

**1. The data is flat; uniform phase space rises ~2x from threshold to 11 MeV.** That is the
original exclusion, now shown bin by bin rather than as a single ratio.

**2. Uniform phase space was never the right expectation.** It counts available final states and
ignores the transition amplitude. In a direct reaction the amplitude is the overlap of the 15C
ground state with the final 14C+n continuum, and its momentum scale is

```
kappa = sqrt(2*mu*Sn) = 46.2 MeV/c      (mu = 876.4 MeV, Sn = 1.2180 MeV)
```

a relative-energy scale of only 1.22 MeV. The honest expectation for direct breakup therefore
**peaks just above threshold and falls** -- it does not rise to 14-16 MeV. A statistical
(phase-space) distribution would be right for compound decay, not for this.

**3. But direct halo breakup fails too, in the opposite direction.** `PS x q/(q^2+kappa^2)^2`,
with the shape fixed entirely by Sn and no free shape parameter, overshoots by **6x** at
Ex = 2-3 and undershoots by **3x** at Ex = 18-20. The data sits between the two extremes.

## ★ The real answer: the continuum shape is not a measured quantity yet

Ex < Sn is kinematically forbidden for any 15C decay, so counts there measure mis-reconstruction
directly:

```
                 floor (Ex -2.5..-1.5)   8-10    18-20   floor-subtracted fall 8-10 -> 18-20
ATTPCROOT vz<1100        69.5            159.8    93.3            3.79x
ATTPCROOT vz<400         26.9            121.2    71.5            2.11x
Spyral    vz<400          6.9            107.9    91.6            1.19x
Spyral    vz<1100        20.8            146.7   117.6            1.30x
```

* **A large part of the continuum is instrumental**, and it lives at large vertex z -- the floor
  drops 2.6-3x when vz < 400 mm is imposed, in **both** codes.
* **The two codes disagree on the shape even at the same vz cut**: at vz < 400 ATTPCROOT's
  floor-subtracted continuum falls 2.11x from 8-10 to 18-20 while Spyral's falls 1.19x. The
  residual is the other cuts (chi2 < 0.2 vs < 5, IC gate on vs off) and the fitter.

**A factor ~2 of reconstruction-dependent shape means no continuum model can be tested yet.**
Comparing a raw yield against a cross-section shape without an acceptance x efficiency calculation
is the error underneath this whole question. That simulation is the blocker -- and it is the *same*
one needed for the resolution function. Until it exists, use an empirical background and quote the
spread as a systematic, which is what `globalfit.C` does.
