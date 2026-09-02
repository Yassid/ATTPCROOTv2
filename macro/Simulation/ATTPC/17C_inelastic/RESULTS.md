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
