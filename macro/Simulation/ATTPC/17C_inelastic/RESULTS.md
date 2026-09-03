# 17C(p,p') and 17C(d,d') — results for the M_n/M_p proposal

Campaign of 2026-09-02: 2 channels × 3 levels × 2 fields, 16000 events per sample, in
`/media/yassid/Seagate Hub/ATTPC/C17_inel`. Reproduce with

```bash
root -b -q 'inel_summary_C17.C'
root -b -q 'decompose_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b285")'
```

χ²/ndf < 5, truth-matched fits, median and IQR throughout; σ ≡ IQR/1.349.

## The headline

**The measurement works, and the thing that makes it work is the vertex beam-energy correction,
which the proposal does not currently mention.** Under the constant-beam-energy analysis σ(E_x) is
0.22 MeV and the three low-lying 17C states are one bump. Taking the beam energy at the
reconstructed vertex — pure software, no hardware — gives **σ(E_x) = 0.098 MeV** for (p,p'), and the
two inelastic states become a resolved shoulder on the elastic peak.

At that resolution a fixed-position three-component fit recovers the two yields, in one day of
beam, to:

| configuration | δ(1/2⁺)/N | δ(5/2⁺)/N |
|---|---|---|
| **17C(p,p') 2.85 T, all angles** | **5.6 %** | **9.8 %** |
| 17C(p,p') 4.00 T, all angles | 4.4 % | 7.8 % |
| 17C(d,d') 2.85 T, θ_lab > 60° | 10.3 % | 15.6 % |
| **17C(d,d') 4.00 T, θ_lab > 60°** | **7.8 %** | **12.1 %** |

at N_elastic/N_217 = 10, against a pure-statistics floor of ~3 % and ~5 %. The proposal's target is
~10 % per state, so **the proton day makes it at the nominal field and the deuteron day needs
either 4 T or a second day.**

## Corrections to the proposal text

1. **The 1/2⁺–5/2⁺ separation is 115 keV, not 130.** The proposal quotes the states at 217 and
   332 keV and the gap at 130 keV; 332 − 217 = 115. (ENSDF: 217 and 331 → 114 keV.)
2. **It is a triplet, not a doublet.** The ground state is 217 keV below the 1/2⁺ and the elastic
   channel is far stronger than either inelastic, so the two peaks to be extracted sit on the flank
   of a much larger one. The Ref. [24] precedent does not cover this — there the states of interest
   are at 6–7 MeV, some twenty resolution widths from elastic.
3. **"300 keV resolution achievable in the AT-TPC" is conservative by 3–7×** once the vertex
   correction is applied: 98 keV integrated, and 44–57 keV in the θ_lab 50–70° window where the
   FRESCO yield sits. Worth saying, because it changes the claim from "we hope to fit a doublet" to
   "we resolve a shoulder and fit it".
4. **The 5/2⁺ yield is about half the quoted "~1000 counts".** From the supplied FRESCO
   cross sections, 1379 (1/2⁺) and 579 (5/2⁺) at 4π in one day; with acceptance folded, **1163 and
   487 detected**. The 5/2⁺ sets the statistical error.
5. **FRESCO ran at 136 MeV, the proposal says 8.37 MeV/u = 142.29.** Not a discrepancy — 142.29
   enters and ~14.6 MeV is lost over the metre, so 136 is the mid-chamber mean — but as written it
   reads as two beam energies. This simulation measures the loss from its own truth:
   **E_beam(z) = 141.51 − 0.01462·z [mm]**, 14.6 MeV over the metre, within 1 % of the (d,p) arm's
   independent CATIMA figure.
6. **The field is never stated.** It matters: see below.

## σ(E_x), the campaign in one table

Vertex-corrected, all angles:

| configuration | σ(g.s.) | σ(217) | σ(332) | sep 217–332, const E_beam | sep, vertex-corr |
|---|---|---|---|---|---|
| 17C(p,p') 2.85 T | 0.098 | 0.096 | 0.093 | 0.26 | 0.61 |
| 17C(p,p') 4.00 T | 0.060 | 0.060 | 0.059 | 0.26 | 0.97 |
| 17C(d,d') 2.85 T | 0.494 | 0.471 | 0.444 | 0.10 | 0.13 |
| 17C(d,d') 4.00 T | 0.352 | 0.338 | 0.344 | 0.12 | 0.17 |

Separation is Δ/(σ₁+σ₂); above ~2 a pair is resolvable, below ~1 it is one bump. **Nothing here
reaches 2**, which is why the deliverable is the constrained fit and not a resolved doublet.

The constant-E_beam column is flat at 0.22 MeV for (p,p') at *both* fields, against a method floor
of 0.206 MeV — the detector is nearly invisible and the assumption of one beam energy for every
vertex is essentially the whole width. This reproduces the 14C(p,p') finding on an independent
channel.

## Resolution is strongly angle-dependent, and that is the whole story

σ(E_x) [MeV] per θ_lab slice, ground state, vertex-corrected:

| θ_lab | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 |
|---|---|---|---|---|---|---|
| (p,p') 2.85 T | 0.252 | 0.163 | 0.118 | **0.057** | **0.044** | 0.062 |
| (d,d') 2.85 T | 1.531 | 1.130 | 0.747 | 0.394 | **0.142** | **0.046** |

because σ(KE) itself varies by two orders of magnitude across the range:

| θ_lab | 20–30 | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 |
|---|---|---|---|---|---|---|
| σ(KE) proton | 0.486 | 0.345 | 0.242 | 0.085 | 0.023 | 0.026 |
| σ(KE) deuteron | 2.981 | 2.311 | 1.480 | 0.802 | 0.267 | 0.046 |

A slow recoil makes a tight helix and is measured well; a fast one is stiff and is not. Since
|dE_x/dKE| is *flat* at 0.533 (p) / 0.563 (d) — the kinematics gate's one solid output — σ(E_x)
simply follows σ(KE).

**This is why the deuteron channel needs an angular cut and the proton channel does not.** The
deuteron is twice as energetic as the proton at the same θ_cm, so its forward-angle tracking is
6× worse; integrated over everything that gives σ(E_x) = 0.47 MeV and a 27 % tail. Above
θ_lab 60° it is 0.05–0.14 MeV, as good as the proton channel.

Cutting at θ_lab > 60° (δ(1/2⁺), δ(5/2⁺) at R = 10):

| | all angles | θ_lab > 60° |
|---|---|---|
| (p,p') 2.85 T | 5.6 / 9.8 | 6.1 / 10.9 |
| (p,p') 4.00 T | 4.4 / 7.8 | 6.8 / 11.9 |
| (d,d') 2.85 T | 13.6 / 25.5 | **10.3 / 15.6** |
| (d,d') 4.00 T | 11.3 / 20.2 | **7.8 / 12.1** |

**Protons: use everything. Deuterons: cut at 60°.** The cut costs statistics and buys resolution,
and the trade only pays where the forward angles were unusable anyway.

## The field

Not the lever it first appeared to be. On the constant-E_beam analysis 4 T changes σ(E_x) by
nothing (0.222 vs 0.222) and costs acceptance. Vertex-corrected it is worth 1.6× on (p,p')
(0.098 → 0.060) and it matters most where the measurement is weakest — the deuteron day, where it
takes the 5/2⁺ from 15.6 % to 12.1 % with the angle cut.

**Recommendation: 4 T if SOLARIS will give it, since it is a magnet setting and not a development.
But 2.85 T is not disqualifying** — the proton day makes the 10 % target there, and the deuteron
day can instead be bought with beam time (a second D₂ day is worth √2, i.e. 15.6 % → 11.0 %).

Acceptance, elastic: 0.827 (p, 2.85 T), 0.798 (p, 4 T), 0.653 (d, 2.85 T), 0.665 (d, 4 T). The
proton value reproduces the 14C(p,p') campaign's 0.826 at the same field and pad plane, which is
the cheapest available check that nothing drifted between the two campaigns.

## The kinematic lines, window by window

`kine_lines_C17.C` shows the same result in the variable an experimenter looks at first — recoil KE
against θ_lab (`plots/kinemap_*.png` is that map on its own, full size). **The three loci are not
separable raw**: each event carries the beam energy at its own vertex, and the 14.6 MeV lost across
the chamber spreads the band by several MeV, an order of magnitude more than the level spacing.
Rebuilding every event at a common reference beam energy collapses that spread.

The projections are in **excitation energy**, so the peaks land on 0, 0.217 and 0.332 MeV and can be
read off directly. That also makes the two channels comparable — the same 115 keV gap is 221–269 keV
of proton recoil but only 207–223 keV of deuteron recoil, so a KE residual silently rescales between
them — and it makes these widths the same quantity as the tables above rather than a proxy.

σ(E_x) [MeV] and, beneath it, the separation Δ/(σ₁+σ₂) of the 1/2⁺ and 5/2⁺, per θ_lab window:

| θ_lab | 30–40 | 40–50 | 50–60 | 60–70 | 70–80 |
|---|---|---|---|---|---|
| (p,p') 2.85 T | 0.152 | 0.114 | **0.052** | **0.047** | 0.067 |
| | 0.38 | 0.52 | **1.15** | **1.20** | 0.88 |
| (p,p') 4.00 T | 0.065 | **0.050** | **0.047** | **0.044** | 0.073 |
| | 0.84 | **1.18** | **1.28** | **1.31** | 0.75 |
| (d,d') 2.85 T | 1.136 | 0.754 | 0.372 | 0.135 | **0.047** |
| | 0.05 | 0.08 | 0.16 | 0.44 | **1.23** |
| (d,d') 4.00 T | 0.792 | 0.534 | 0.270 | **0.074** | **0.047** |
| | 0.07 | 0.11 | 0.23 | **0.76** | **1.21** |

**This refines the angular recommendation made from the integrated numbers.** The usable window is
narrower and further back than "θ_lab > 60°":

- **(p,p') 2.85 T: 50–70°**, with 70–80° usable but degraded.
- **(p,p') 4.00 T: 40–70°** — the field buys the 40–50° bin outright, which is real added coverage.
- **(d,d') 2.85 T: 70–80° essentially alone.** The 60–70° bin gives 0.44, not the ~1 that a
  θ_lab > 60° cut implied. That is a much tighter window than the proposal's angular coverage
  assumes, and it is the strongest argument for 4 T.
- **(d,d') 4.00 T: 60–80°**, roughly doubling the deuteron coverage.

Note that at 70–80° the deuteron reaches σ(E_x) = 0.047 MeV — identical to the proton's best. The
deuteron channel is not intrinsically worse; it is worse everywhere *except* the backward window,
because its recoil is twice as energetic and only becomes slow enough to track well at the very
back.

**Correction to an earlier statement in this document:** the 70–80° (p,p') bin was previously quoted
as separation 0.30, a collapse. That number came from a KE-residual width that was not a proper
residual for the excited levels; in excitation energy it is 0.88. Resolution does still stop
improving at the last bin — the recoil falls below ~2 MeV and the track gets short — but the bin
remains usable rather than lost.

The figures are `plots/kinelines_*.png`: the locus map, the same in E_x (flat bands at the level
energies), and the five angular windows on a log axis with the components drawn separately. The
elastic there is drawn at an **assumed** N_elastic/N_217 = 10 with the 1/2⁺ angular shape — the
1/2⁺ : 5/2⁺ ratio of 2.38 is real, the elastic normalisation is not, and it is labelled on the
figure.

## What is still missing

1. **The elastic FRESCO curve at 136 MeV.** R = N_elastic/N_217 is the single number that decides
   the result and no calculation was supplied, so it is scanned. The dependence is steep: at
   R = 10 the (p,p') 5/2⁺ is 9.8 %, at R = 30 it is 12.7 %, at R = 100 it is 21.7 %. **Anything
   above ~30 breaks the 5/2⁺ and the proposal should not be submitted without knowing which side
   of that the real number falls on.**
2. **The (d,d') angular distributions.** The proposal asserts "a similar amount for deuterons";
   the deuteron numbers above use the *proton* cross sections for want of anything else.

## Caveats on these numbers

- The yields are **angle-integrated**. That is the right figure of merit for the deformation length
  and hence M_n/M_p, which is a normalisation of the whole angular distribution; the per-bin points
  needed to confirm the shape will be ~√N_bins worse.
- The toy fit is a **linear least squares with per-bin variance max(n,1)**, not a Poisson
  likelihood. It leaves a −1 to −5 % bias on the 5/2⁺ (worse, −10 to −14 %, once the angle cut
  thins the bins), which a proper Poisson likelihood would mostly remove. The quoted *widths* are
  reliable; the biases are an artefact of the estimator and should not be read as a physics
  systematic.
- Angular distributions are **flat in cos θ_cm** at generation and reweighted by FRESCO in the
  analysis, so a revised DWBA calculation reweights these results rather than requiring a new
  campaign.
- `MEASSIGMA = 4.0` mm, matching the 14C(p,p') reference so the comparison holds. It is not the
  measured hit residual (0.59–0.64 mm on this pad plane), at which the χ²/ndf < 5 cut would
  actually bite instead of passing nearly everything. One cell with `MEASSIGMA=0.6` is the control
  to run before quoting a tail fraction.

## Traps hit while building this, worth not repeating

1. **`acceptance_C14.C` and `ex_res_C14_hf.C` truth-match the ejectile by hard-coded PDG 2212.** On
   (d,d') that matches nothing, and both report `generated reactions 0 ... acceptance 0.000` on
   six samples whose genfit output was perfect (5965 fitted tracks). Both now take an `ejPdg`
   argument defaulting to proton, so every (p,p') and (d,p) caller is byte-identical.
   `check_beam_C17.C` had the same defect and the same fix.
2. **`decompose_C17.C` labelled its output "vertex-corrected" while building templates from the raw
   `exReco`,** which is the constant-E_beam reconstruction. The label was right about intent and
   wrong about fact, and it overstated the difficulty by the full 2.3× the correction is worth.
   The correction has to be *redone* from the tracked (θ, KE) at E_beam(z_reco); there is no
   corrected column in the tree.
3. **Normalising to the 4π cross section without folding acceptance** treats all 1379 counts as
   detected and understates the errors by ~10 %. The generator is flat in cos θ_cm, so the
   acceptance per angular bin is recoverable from the accepted θ_cm distribution and the
   generated-reaction count in the acceptance log.
4. **`gSystem->GetFromPipe("ls -1 \"" + pattern + "\"")` does not glob.** The products live under
   `.../Seagate Hub/...`, so the path needs quoting for the space — but quoting the whole pattern
   stops the shell expanding the wildcard and every sample silently reads as MISSING. Quote the
   directory only.
5. **Extrapolating one angle-averaged σ(KE) from another campaign predicts the wrong answer.** The
   kinematics gate fed 14C(p,p')'s single elastic σ(KE) = 0.343 MeV through this channel's leverage
   and concluded σ(E_x) was flat in angle at 0.18 MeV and that the field was the only lever. Both
   were wrong: σ(KE) varies 20× across the angular range, and the constant-E_beam term dominates
   until it is corrected. The leverage half of the gate (dE_x/dKE, dE_x/dθ, dE_x/dE_beam) is exact
   and was worth computing; the σ(E_x) prediction was not, and the macro now says so in its header.

## Realistic counts: the ledger (2026-09-03)

`yields_C17.C` writes out every factor between a FRESCO cross section and a fitted differential
point, each one either **measured** from the campaign or **named** as a parameter. Reproduce with

```bash
root -b -q 'yields_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b285")'
```

The proposal's "∼1000 counts of scattered protons" is not derived from a cross section — it is an
analogy to Ref. [24] scaled by the beam rate. Replacing it with arithmetic:

| step | 1/2⁺ | 5/2⁺ | source |
|---|---|---|---|
| σ(4π) = 8.582 / 3.601 mb × 940 pps × 86400 s × 1.978e21 cm⁻² | 1379 | 579 | FRESCO + gas |
| × duty 0.70 × purity 1.00 | 965 | 405 | **[par]**, no source supplied |
| × pile-up-free 0.930 (77 µs drift at 940 pps) | 898 | 377 | [par], computed |
| × usable angular window × acceptance × χ² × z-fiducial | **352** | **145** | **[sim]** |

So **the counts that can actually enter the fit are ~350 and ~145, i.e. 0.35× and 0.15× the quoted
~1000** for (p,p') at 2.85 T. Two factors do all of that work and neither is in the proposal:

1. **A requested day is not 86400 s of beam on target.** Every yield in this study and in the
   proposal multiplies 940 pps by a full day. The simulation cannot know the duty factor; it is a
   parameter with a *guessed* 0.70 default and it is the cheapest thing in the whole study to fix —
   it needs a number from ReA6 operations.
2. **Counts outside the usable angular window are not counts.** The FRESCO distributions put only
   **42.7 %** of the (p,p') yield in the 50–70° window that the resolution study leaves usable
   (53.7 % in the 40–70° window at 4 T, 13.0 % in the (d,d') 70–80° window at 2.85 T, 37.3 % at
   4 T). The earlier "1163 and 487 detected" integrates over angles whose E_x spectrum is a
   featureless bump.

| config | window | N(1/2⁺) | N(5/2⁺) | δ(1/2⁺) | δ(5/2⁺) | best per-point δ(5/2⁺) |
|---|---|---|---|---|---|---|
| (p,p') 2.85 T | 50–70° | 352 | 145 | 7.3 % | 11.5 % | 21.0 % (60–65°) |
| (p,p') 4.00 T | 40–70° | 429 | 179 | 6.0 % | 9.5 % | 18.9 % (60–65°) |
| (d,d') 2.85 T | 70–80° | 91 | 36 | 12.7 % | 20.5 % | 25.7 % (70–75°) |
| (d,d') 4.00 T | 60–80° | 241 | 106 | 10.5 % | 16.7 % | 23.9 % (65–70°) |

at duty 0.70, purity 1, one day, R = 10. The window-integrated errors are close to the earlier
all-angle numbers (5.6 / 9.8 % for (p,p') 2.85 T) — the forward angles the window cut removes were
contributing little — but the **per-angular-point** column is new and it is the one that matters:
M_n/M_p is the normalisation of a dσ/dΩ *curve*.

**The per-point conclusion is the harshest result in this study.** Reaching 10 % on the 5/2⁺ at
even the best single angular point needs 3.6–6.6× the integrated luminosity of one day at duty
0.70. That is more than duty = 1 can give back, so with the requested 2 days the deliverable is a
**window-integrated yield ratio**, not a measured angular distribution. Either the proposal should
say that — the deformation length is a normalisation, so an integrated ratio is a defensible
deliverable — or it should ask for more days, or the window has to widen, which is the strongest
remaining argument for 4 T (0.427 → 0.537 of the yield on the proton day, 0.130 → 0.373 on the
deuteron day, i.e. **2.9× the deuteron statistics**, on top of the 1.6× resolution gain).

Days-equivalent value of each factor, (p,p') 2.85 T: window ×2.34, duty ×1.43, pile-up ×1.08,
purity ×1.00. **The window is the largest and the only one a detector setting can move.**

`plots/yields_*.png` is δN/N per angular point against θ_lab, with the 10 % target drawn.

### What the ledger still cannot include

- **The elastic angular distribution.** R is still scanned and still assumed *constant in angle*,
  which is wrong in a known direction: elastic is forward-peaked, so the true R rises towards the
  forward end of every usable window — exactly where the per-point errors are already worst.
- **The trigger.** The AT-TPC trigger counts *pads* above threshold, and a backward-angle recoil in
  a 2.85–4 T field is a tightly curled helix that revisits the same pads. The digitised events are
  already on disk in the campaign's sim files and nothing has counted their pad multiplicity against
  a threshold. This would bite hardest in precisely the backward window the resolution study
  selected, so it is the cheapest unrun check in the study.
- **Competing channels and the 17O contaminant**, neither of which exists in the simulation at any
  rate. On the deuterium day (d,p), (d,t), (d,³He) and breakup all compete.
- **The assumed deformation length behind the FRESCO cross sections.** Every count above scales with
  it. If the DWBA β is 20 % smaller the yields fall by ~36 %, which moves the ledger further than
  any detector factor in it.

## What the two angular distributions would actually look like (2026-09-03)

`dsdo_C17.C` draws the supplied DWBA curves with the data points and error bars this experiment
would put on them. The elastic is left arbitrary (R = 10, flat) — it enters only through the
overlap penalty on the bars.

```bash
root -b -q 'dsdo_C17.C("/media/yassid/Seagate Hub/ATTPC/C17_inel","pp_b400",40,70,10)'
```

**The first thing the curves say, before any detector: the two shapes are nearly proportional.**
dσ(217)/dσ(332) runs only **1.80 → 2.62 across θ_cm 2–150°** — flat to ±15 % about 2.2. Both have
their first maximum at **θ_cm 45° (θ_lab 68°)**, the same diffraction minimum at **θ_cm 93°
(θ_lab 44°)** and the same secondary maximum near θ_cm 130°. The 5/2⁺ is the 1/2⁺ divided by ~2.2
with identical structure.

So **the angular shape carries essentially no power to separate the two states.** What the
measurement extracts is two *normalisations*; the shape's only job is to confirm the assumed L
transfer. That reframes the deliverable, and it is good news for the count problem: a window that
covers one maximum well beats a wider window that smears both.

**And the first maximum lands inside the good-resolution window.** θ_cm 45° is θ_lab 68°, which sits
in the 50–70° window on the proton day — a piece of luck the proposal does not currently claim.

Expected points, one day, duty 0.70, R = 10, **10° θ_lab bins**:

| config | θ_lab | θ_cm | N(1/2⁺) | N(5/2⁺) | dσ/dΩ(217) | ± | dσ/dΩ(332) | ± |
|---|---|---|---|---|---|---|---|---|
| (p,p') 2.85 T | 50–60 | 60–80 | 161 | 66 | 0.769 | 10.7 % | 0.316 | 18.6 % |
| | 60–70 | 40–60 | 191 | 79 | 1.251 | 8.4 % | 0.537 | 15.5 % |
| (p,p') 4.00 T | 40–50 | 80–100 | 92 | 36 | 0.431 | 15.9 % | 0.179 | 27.9 % |
| | 50–60 | 60–80 | 148 | 66 | 0.769 | 9.7 % | 0.316 | 15.8 % |
| | 60–70 | 40–60 | 188 | 78 | 1.251 | 8.6 % | 0.537 | 15.0 % |
| (d,d') 4.00 T | 60–70 | 40–60 | 170 | 73 | 1.251 | 13.6 % | 0.537 | 22.7 % |
| | 70–80 | 20–40 | 71 | 34 | 1.019 | 14.1 % | 0.457 | 23.4 % |
| (d,d') 2.85 T | 65–75 | 30–50 | 135 | 57 | 1.234 | 14.9 % | 0.538 | 22.4 % |
| | 75–85 | 10–30 | 36 | 14 | 0.741 | 21.7 % | 0.351 | 35.5 % |

At 5° binning the same windows give 12–18 % (1/2⁺) and 21–31 % (5/2⁺) per point — see
`plots/dsdo_*_d5.png`. **10° bins are the right choice**: nothing in a curve this smooth needs 5°
sampling, and the shape is fixed by DWBA anyway.

**The coverage, config by config:**

- **(p,p') 2.85 T** (θ_lab 50–70 = θ_cm 40–80): covers the first maximum and its fall-off. **Two
  points.** Misses the minimum.
- **(p,p') 4.00 T** (40–70 = θ_cm 40–100): covers the maximum *and reaches the diffraction
  minimum at θ_cm 93°*. **Three points spanning peak → minimum**, and that is the only real shape
  information the experiment can get. A third argument for 4 T, independent of resolution and of
  statistics.
- **(d,d') 4.00 T** (60–80 = θ_cm 20–60): the maximum and the backward side of it, no minimum.
- **(d,d') 2.85 T** (70–80 = θ_cm 20–40): one, at best two usable points, all on the far side of
  the maximum. Sliding the window back to 65–85° gives two points at 14.9 / 22.4 % and
  21.7 / 35.5 %, so even there the deuteron day at 2.85 T is a one-good-point measurement.

**What the proposal should therefore claim.** Not "differential cross section curves to both the
1/2⁺ and the 5/2⁺", which two or three points on a shape that is degenerate between the states
does not support. The defensible claim is **two normalisations, measured at 2–3 angular points that
bracket the DWBA maximum, to 8–11 % (1/2⁺) and 15–19 % (5/2⁺) in one day** — which is exactly what
M_n/M_p needs, since the deformation lengths are normalisations. Stated that way the proposal's
~10 % per state is met on the proton day and nearly met on the deuteron day at 4 T, and the
per-point pessimism of the previous section stops being the binding constraint.

`plots/dsdo_<cfg>_<window>_d<bin>.png`: left panel the two DWBA curves with the expected points
(inner bar statistics only, outer bar with the three-component overlap penalty) and the usable
window shaded; right panel the ratio of the two curves, which is the figure that shows the shapes
are degenerate.
