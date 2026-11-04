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
