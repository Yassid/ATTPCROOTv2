# Is the Spyral PID plane better? — measured 2026-11-04

Prompted by the impression that Spyral's plane separates the **triton** band better and that its
**deuteron** band is thinner. Short answer: **the deuteron band is the same width in both, the
triton sits below its local background in both, and most of the visual difference is statistics
and axis scaling.** Images:

| file | what it shows |
|---|---|
| `~/C15d_plane_compare.png`   | both planes, matched axes, each with its own statistics |
| `~/C15d_plane_equalstats.png`| both planes downsampled to the SAME 476,704 tracks, identical colour scale |
| `~/C15d_triton_region.png`   | zoom on the triton region with the Bethe p/d/t loci overlaid |

## The bands are the same bands

The two planes are on different dE/dx scales — **Spyral / ATTPCROOT = 1.726** (median sqrt_dEdx in
Brho 0.6-0.8: 18.76 vs 10.87). Nothing can be compared by eye until that is divided out; the two
GUIs also default to different x ranges (90 vs 45), which alone changes how tight the bands look.

Identified by anchoring on the proton and predicting the rest from Bethe (same momentum at a given
Brho, so the ratios follow from beta alone):

```
ATTPCROOT Brho 0.60-0.75   peaks 5.4  10.9  15.2 ...   predicted  d 9.7   t 13.6
Spyral    Brho 0.60-0.75   peaks 5.7  11.2  16.1 ...   predicted  d 10.2  t 14.2
observed t/d = 1.39        predicted t/d = 1.40
```

So the third band **is** the triton, in both.

## Deuteron width: no difference

Single-Gaussian-plus-linear fits at each band's own peak, five Brho slices:

```
sigma/mu    0.40-0.50  0.50-0.60  0.60-0.70  0.70-0.85  0.85-1.05   mean
ATTPCROOT     0.105      0.085      0.090      0.099      0.111     0.0979
Spyral        0.077      0.097      0.087      0.089      0.124     0.0948
```

3 % apart, with slice-to-slice scatter larger than the difference. **Not a real difference.**

## Triton: below background in both

The deuteron-triton separation `R_dt = (mu_t - mu_d)/(sig_d + sig_t)` could only be measured in 1
slice for ATTPCROOT (2.71) and 3 for Spyral (2.39) — the fits do not converge elsewhere, because
the triton peak-to-background ratio is **0.0-0.3 in both planes**. A band that sits below its own
local background cannot be said to be "better separated" in one code than the other on this data.

## Band purity: ATTPCROOT is not worse

Fraction of tracks in a Brho slice lying within +-2 sigma of a fitted p/d/t peak:

```
overall in-band fraction      ATTPCROOT 65.7 %   Spyral 56.8 %
overall gaussian-area frac.   ATTPCROOT 60.1 %   Spyral 49.5 %
```

⚠ One Spyral slice returned 102 % and 145 %, i.e. overlapping fits double-counting, so treat the
per-slice numbers as indicative and the direction (ATTPCROOT not worse) as the result.

## What the impression probably was

ATTPCROOT's plane has **1.9x more tracks** (2,268,873 vs 1,732,192; 823,742 vs 476,704 after each
plane's own IC cut, from 102 runs vs 85). More tracks means a proportionally denser background
filling the gaps between bands, at identical band width. Downsampled to equal statistics with the
same colour scale (`C15d_plane_equalstats.png`) the two look very similar.

**Conclusion: no evidence that switching the gate to the Spyral plane would improve the PID.**
The ATTPCROOT plane has nearly twice the statistics at equal band quality.

## What the gain matching does and does not explain

Asked separately, and worth recording: **the per-run gain matching is equivalent.**

```
Spyral    83 runs, factors 0.638-1.175, spread 1.84x, median 1.021
ATTPCROOT 102 runs, factors 0.679-1.252, spread 1.84x, median 1.060
correlation of the two per-run shapes on 80 common runs: +0.983, residual scatter 3.4 %
```

Both track the same physical run-to-run drift. The real differences are elsewhere, and two of them
were tested and do NOT change band width:

- **charge per hit**: Spyral integrates the pulse, ATTPCROOT uses the peak (`GetCharge()`;
  `AtSpyralPID::fUseTraceIntegral` defaults false and `AtPIDTask` does not expose it). Switching it
  on is a **pure 1.72x rescale** at track level — every scale-free width is identical to three
  decimals, because averaging ~130 hits washes out the per-hit difference.
- **dE/dx integration region**: ATTPCROOT truncates at `fSmallPadRadius` = 152 mm so the two pad
  sizes never mix; Spyral uses the whole first arc with a per-pad size gain (1.3 vs 1.0). Widening
  ATTPCROOT's radius makes the bands **worse** (IQR/median 0.896 -> 0.982 at 250 mm), so 152 mm is
  already the better choice.

Still different and untested: the IC window (**[1045,1225] vs [1000,1250]**), the clustering, and
Spyral's Garfield field correction.

## Open: there is a FOURTH band

Both planes show a strong band beyond the triton. It is present in BOTH codes, so it is in the
data rather than in a reconstruction, and it is bright enough to matter for a yield.

**It is not yet identified, and the obvious method does not work.** Measuring its sqrt_dEdx ratio
to the proton per Brho slice gives inconsistent answers, because the per-slice peak finder is not
reliable past the first three bands:

```
4th/proton ratio   Brho 0.45-0.55  0.55-0.65  0.65-0.80  0.80-0.95
  ATTPCROOT             3.57         2.56        3.62       5.27
  Spyral                2.58         2.84        2.88       3.99
  predicted:  t 2.49    3He 2.81     4He 3.57    6Li 5.36
```

The values bracket 3He and 4He but the scatter is far too large to claim either, and in the Spyral
0.45-0.55 slice the finder split the PROTON band into two peaks (7.5 and 8.3), which shows the
method is at its limit. Identifying this band needs continuous ridge-following in 2D rather than
per-slice peak finding. **Left open** -- do not quote a yield that includes this region until it
is settled.

## Reproduce

```bash
S=xcheck
root -b -q $S/pid_identify_bands.C    # which band is which, from Bethe
root -b -q $S/pid_band_metrics.C      # per-band width and separation
root -b -q $S/pid_purity.C            # in-band fraction
root -b -q $S/pid_plot_equalstats.C   # the equal-statistics image
root -b -q $S/pid_plot_triton.C       # the triton zoom with loci
```

⚠ Two traps hit while doing this, both of which produced plausible wrong answers:
- a max-bin FWHM on ~300 entries in 200 bins measures fluctuations, not width (it claimed the
  trace integral narrowed the bands 3-7x; the robust IQR showed it changes nothing);
- clipping a histogram range at 3 or 5 when sqrt_dEdx runs to 60 piles everything into overflow
  and makes the distribution look ten times smaller than it is.

---

# Overnight follow-up, 2026-11-04 night

Built a **rail-guided ridge follower** (`~/Data/C15d_pidstudy/ridge.C`) because every earlier
measurement was limited by per-slice peak finding. Free-running continuity DRIFTS -- seeded on the
triton at 16.08 it slid onto the deuteron at 9.5-11.2 and happily reported 24 bins of the wrong
band. The fix is to search near the band's own **Bethe locus**, whose normalisation is fixed once
from the proton: the locus is a rail the search cannot leave, so a lost bin is dropped instead of
dragging the chain onto a neighbour.

## The planes are equivalent — now with an error bar

```
deuteron sigma/mu   ATTPCROOT 0.0683   Spyral 0.0801   difference -0.0118 +- 0.0092  (1.3 sigma)
peak/background     p 4.18 / 3.00      d 3.64 / 3.84      t 1.08 / 1.09
Brho distribution   medians 0.461 vs 0.403, quartiles within a few percent
IC window           harmonising [1045,1225] vs [1000,1250] moves nothing material
```

The deuteron width difference is **1.3 sigma, i.e. not significant**, the intensity-independent
S/B is equal band for band, and the two codes have the same Brho acceptance. Combined with the
earlier results (gain matching equivalent, trace integral a pure rescale, 152 mm already the best
integration radius), **there is no measurable PID-quality difference between the two planes.**
The visual impression is 1.9x statistics plus a 1.726x scale and a 90-vs-45 axis range.

## ★ THE FOURTH BAND IS 4He — and it is the strongest band in the plane

Identified by asking which (A,z) rail each band follows, with the dE/dx scale calibrated once from
the proton so every candidate locus is fully predicted:

```
ATTPCROOT   species  mu@Brho 0.60   bins   3-sigma   <S/B>
            p             6.15        27      18      2.38
            d            10.99        27      20      1.19
            t            15.33        14       6      1.55
            3He          17.28        18       5      0.64
            4He          21.98        59      56      5.55    <- strongest of all
```

Confirmed by SHAPE, not just position -- t, 3He and 4He rails differ in Brho dependence, so
following the band over 40 bins discriminates:

```
mean |residual| over 40 Brho bins:   t 35.2 %   3He 20.7 %   4He 5.7 %   6He 32.4 %
```

15C(d,alpha)13B is a real channel and alphas have high dE/dx, so this is physical. See
`~/C15d_bands_identified.png`.

## ★★★ THE SAVED DEUTERON GATE IS MORE TRITON THAN DEUTERON

The consequence of having the bands identified. Counting tracks inside `pid/deuteron_C15d.json`
within +-2 sigma of each rail (the d and t rails sit ~39 % apart in sqrt(dE/dx) at every Brho, so
they cannot be confused -- only 0.1 % of tracks are ambiguous):

```
inside the saved gate (IC 1045-1225, single pulse): 450,931 tracks
  deuteron rail only   169,948   37.7 %
  TRITON   rail only   193,058   42.8 %   <-
  alpha    rail only     5,315    1.2 %
  proton   rail only       855    0.2 %
  ambiguous                636    0.1 %
  background            81,119   18.0 %
  -> deuteron purity 46.0 %
```

`~/C15d_gate_vs_bands.png` shows the gate polygon enclosing both the deuteron and triton rails.
**Every yield and spectrum drawn through this gate is a d+t mixture**, and the triton tracks are
fitted under a deuteron mass hypothesis.

A physics-based alternative following the deuteron rail at +-2.5 sigma is written to
`pid/deuteron_rail_proposal.json` (the hand-drawn gate is untouched):

```
                      tracks    deuteron   triton   background
saved gate           450,931     37.7 %    42.8 %     18.0 %
rail proposal        175,095     76.8 %     0.8 %     22.5 %
```

Twice the purity, triton contamination essentially gone, at 39 % of the statistics (~134k real
deuterons against ~170k, so most of the loss is contamination rather than signal). It is a
PROPOSAL, not a replacement -- a hand-drawn gate can follow the band's real shape better than a
constant-width rail, and the 22.5 % background inside it says the rail is not the last word.

---

# You were right: Spyral DOES separate better. Here is why. (2026-11-05)

My earlier "no measurable difference" was measured the wrong way -- per-Brho-slice fits, averaged.
Redone with the **a1975 (d,p) projection** (`pid_width_plot.C`): plot

    ln( sqrt(dE/dx) * Brho )

For a particle of charge z, dE/dx ~ z^2/beta^2 ~ z^2 m^2/p^2 and p = z*Brho*c, so
**dE/dx ~ m^2/Brho^2 and sqrt(dE/dx)*Brho ~ m with the charge cancelling exactly.** The whole
plane collapses into one spectrum, every band becomes a peak, all the statistics are used at once,
and the peak spacings are known constants. That is what my slice-averaging was washing out.

## The measurement

Local maxima, spacings quoted from the dominant (deuteron) peak:

```
                        p          t        3He       4He
expected (full Bethe) -0.578     +0.330    +0.455    +0.693
ATTPCROOT             -0.590     +0.410              +0.700
Spyral                -0.560     +0.530              +0.900
```

**The 4He spacing is a parameter-free check**: in this variable it must be exactly
ln(4/2) = +0.693 at EVERY Brho, because the charge and mass dependencies cancel. ATTPCROOT gives
+0.700 (1 % off). **Spyral gives +0.900 -- 30 % too wide.**

So the bands really are further apart in Spyral, by 20-30 % in ln. Confirmed track by track on the
433,453 tracks both codes reconstruct from the same events:

```
region        ln(ATTPCROOT)    ln(Spyral)     stretch
p                 -0.605         -0.582        x0.96
d                 -0.006         -0.010        x1.00
t                 +0.402         +0.482        x1.20
alpha             +0.698         +0.867        x1.24
```

Spyral's dE/dx is **non-linear**: flat at the deuteron, expanding above it. The deuteron peak WIDTH
is identical on those same tracks (sigma 0.0946 vs 0.0948), so this is extra SEPARATION, not better
resolution -- and it is why the plane looks better.

## ★ The cause: ATTPCROOT's "trace integral" is not an integral

`AtSpyralPID` can use `AtHit::GetTraceIntegral()` instead of the peak amplitude, and the comment
says that matches Spyral's scale. It cannot, in this configuration. `AtPSAMultiFit` sets

```cpp
double amp = A * fRespPeakVal;            // peak amplitude
double Q   = A * fRespIntegral * tauTB;   // "integral"
hit->SetTraceIntegral(Q);
```

with `tauTB` the peaking time, which is **fixed** unless `fFloatTau` or `fFitChi2Cut` is set --
and `AtPSAMultiFit{}` defaults `fFloatTau{false}`, which the C15d macro never changed. So Q is a
constant multiple of amp. Measured on 57,151 hits of run_0016:

```
integral/peak :  min 2.9618   median 2.9618   max 2.9618   rms/mean 0.0000
                 and 2.9618 in every charge quintile from <q> 87 to <q> 3225
```

**Exactly constant.** ATTPCROOT's dE/dx is therefore a pure AMPLITUDE measurement and
`SetUseTraceIntegral` is a no-op up to an overall scale -- which is why switching it on changes
nothing at all (the ln spectrum is byte-identical once normalised). Spyral integrates the real
trace, so its dE/dx carries the pulse WIDTH as well, and the width grows with the deposit; that is
the natural candidate for the non-linear expansion at high ionisation.

`unpackReco_C15d.C` now takes **`psaFloatTau`** (appended LAST so every positional caller keeps
working; default off = the production). With it on, `tauTB` floats per pulse and the trace integral
becomes a genuine independent measurement.

## Which code is right?

On the 4He parameter-free check, **ATTPCROOT is linear in the true dE/dx and Spyral is not.**
That has opposite implications depending on the use:

- **for GATING, Spyral's stretch is an advantage** -- the bands are further apart at equal width,
  so a gate is easier to draw and leaks less;
- **for anything quantitative in dE/dx**, ATTPCROOT is the one on a linear scale, and Spyral's
  high-ionisation response is compressed into a 30 % over-spread.

Neither is "better" outright, and the choice should be made deliberately rather than by which plot
looks nicer.
