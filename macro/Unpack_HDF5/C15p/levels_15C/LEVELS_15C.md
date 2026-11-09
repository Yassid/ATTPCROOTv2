# 15C(d,d'): fitting the states seen in 15C(p,p')

Done 2026-11-09 on the user's saved explorer parameters
(`~/Downloads/explorer_params_15C_d_d__2026-11-09-14-30-49_{atttpcroot,spyral}.json`).
Figures: `~/C15d_dd_statefit.png`, `~/C15p_pp_inelastic.png`, `~/C15d_dd_ex_both.png`.
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
