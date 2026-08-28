# C15d — 15C + d, D2 target at 300 torr

**Beam is ¹⁵C** (confirmed 2026-08-28), so the channel is **¹⁵C(d,p)¹⁶C, Q = +2.026 MeV**. Note
that `a1975/D2_UKF/README.md` describes the same run range as ¹⁶C(d,p)¹⁷C — that is a different
reading of this dataset and is **not** what this workspace analyses. Do not "correct" the masses
against it.

Self-contained workspace. It shares **no macro, no par file and no output directory** with
any other analysis in this repo; the only framework additions it relies on are the new
`AtGainMatchTask` (opt-in, see below) and the `ATTPC_D300torr_v2` geometry.

Source the environment first:

```bash
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh
```

## Data and run set

| | |
|---|---|
| Raw HDF5 | `/media/yassid/Seagate Hub/ATTPC/Data/a1975/h5/run_XXXX.h5` |
| Runs | **17–103, 75 runs** (`runs_d2.txt`) — the D2 target set |
| Format | legacy remerged (`/get` + `/frib` + `/meta`) — `AtHDFUnpacker` reads it directly; no merger-format reader needed |
| Excluded | **runs ≥106 — HYDROGEN target, see below** (`runs_h2_EXCLUDED.txt`); runs 104–105 (thin changeover); **run_0090** (truncated HDF5); **run_0047** has no `/frib`, so no IC |
| Outputs | `~/C15d_reco`, `~/C15d_ic`, `~/C15d_fit` — symlinks onto the Seagate, so no path here contains a space |

### ★ Runs ≥106 are a hydrogen target

The D2 target set ends at run 103. Measured at the boundary:

| | runs ≤103 | runs ≥106 | ratio |
|---|---|---|---|
| dE/dx median | 240–256 | 122–136 | **0.52** |
| arclength median | 162–166 mm | 182–185 mm | 1.13 |

H2 and D2 at equal pressure have the same number density and half the A, so H2 should give
**0.500** of the dE/dx — and less stopping means longer tracks. Both are target changes, not gain.
Every macro defaults to `runMin=17, runMax=103`; a gain match spanning that boundary would force
the H2 runs onto the D2 scale and merge two different datasets.

## Working point

`parameters/ATTPC.C15d_D2_300torr.par`, taken from the Spyral configuration that produced
the reference analysis (`15CdAnalysis/RunScript_d.py`) so both pipelines reconstruct z the
same way: B = 2.85 T, E = 45000 V/m, drift length 1000 mm, micromegas TB 25, window TB 300,
3.125 MHz (`SamplingRate 3` → 320 ns/TB). Drift velocity is then **not free**: 1.13636 cm/µs
is forced by those anchors. Density 6.5643e-5 g/cm³ matches `media.geo`'s `TargetD2_300`
exactly, so the UKF's CATIMA model and genfit's TGeo material describe the same gas.

**Known difference from Spyral:** Spyral additionally applies a Garfield drift-field
correction (`do_garfield_correction`, a (ρ,z) → Δρ/Δtrans/Δz grid). ATTPCROOT has no
equivalent, so a residual z non-linearity between the two is expected and is not a bug.

## Pipeline

```bash
# 1. geometry, once
root -b -q ../../../geometry/ATTPC_D300torr_v2.C

# 2. reco:  raw -> AtPatternEvent
root -b -q 'unpackReco_C15d.C("run_0017", -1, false, "/home/yassid/C15d_reco/")'

# 3a. GENFIT with the CATIMA material model (material effects ON by default)
root -b -q 'fitGenfit_C15d.C("run_0017", -1, "/home/yassid/C15d_reco/", "", "/home/yassid/C15d_fit/")'

# 3b. UKF, same reco and same PID plane -> a controlled comparison on the fitter alone
root -b -q 'fitUKF_C15d.C("run_0017", -1, "proton", "/home/yassid/C15d_reco/", "/home/yassid/C15d_fit/")'
```

Ejectile defaults to the **proton**, i.e. the (d,p) channel. Both fitters take the species as
arguments for (d,d) / (d,t).

### GENFIT + CATIMA

`fitGenfit_C15d.C` defaults to `matEffects = kTRUE` with `catimaMSC`, `catimaStraggling` and
`catimaELoss` all on. Three things that silently defeat it:

1. **The CATIMA flags are inert with `matEffects = kFALSE`** — the noise terms are never
   reached, so an A/B against the material-effects-off configuration reads as a perfect null.
2. **MSC and straggling must both be on**, or genfit's own model still handles the other term.
3. **It needs a GenFit built with `-DGENFIT_USE_CATIMA=ON`.** Check:
   `ldd $GENFIT/lib/libgenfit2.so | grep catima`. This install (`~/fair_install/GenFitInst`,
   branch `catima-scattering`) links `~/fair_install/catima-inst`.

`matFallback` defaults to **off** so a failed material-effects fit drops out instead of being
silently refitted without material effects and kept.

Backward tracks are kept (θ window 5–178°, `backwardSeedFix` on): in (d,p) inverse kinematics
the proton goes largely backward. **The B sign is unverified for this run set** — it is an
argument; pick it on the fitted vertex and a physical KE, not on χ² alone.

## Gain matching — measured from this analysis's own data

The micromegas/GET gain drifts run to run, so the same particle deposits a different measured
charge in different runs. Uncorrected, one PID gate cannot serve the run set.

**The correction** (`AtGainMatchTask`, opt-in FairTask, added *after* `AtPIDTask`):

```
dEdx *= f(run)        sqrtdEdx *= sqrt(f(run))
```

It modifies no existing framework class — it goes through `AtPIDEvent`'s public `Clear`/`Add`
API — and it is generic: table path and run number are arguments, so any experiment needing gain
matching can add it. Verified exact on run_0017: every track scaled by the table value and its
square root, none lost.

**The factors** (`measure_gain_C15d.C`) are measured here, from this plane. The anchor is a
quantity that should be constant across runs, so any drift in it is gain:

```
sel(r)  = valid && brho in [bLo, bHi]
peak(r) = location of the dE/dx distribution of sel(r)
f(r)    = peak_ref / peak(r)
```

The Bρ window makes it **species-anchored**: inside a slice where one band dominates, the true
dE/dx of those tracks is physics and does not depend on how many of them there were, so a run at
a different rate or with a different overall species mix gives the same peak unless the gain
moved. Anchoring on the median of *everything* fails exactly here — a run with more deuterons
has a higher median for reasons that are not gain, and that shift would be absorbed into `f(r)`.

**Choose the window off the plane first.** It must sit where a single band dominates; a window
straddling two bands measures their ratio, not the gain, and looks perfectly smooth while being
wrong. Profiling the cached plane in Bρ slices gives two clean candidates:

| Bρ window | N | MPV | median | shape |
|---|---|---|---|---|
| 0.20–0.25 | 2806 | 23.67 | 25.78 | 9 peaks |
| **0.25–0.30** | 2281 | 21.00 | 21.09 | **SINGLE** |
| **0.30–0.40** | 2431 | 17.75 | 17.77 | **SINGLE** ← adopted |
| 0.40–0.50 | 835 | 14.33 | 14.40 | 4 peaks |

```bash
root -b -q 'measure_gain_C15d.C(0.30, 0.40)'
```

### ★ Use a robust estimator, not the peak bin

Measured on runs 17/19/20/21 (~790 tracks each, Bρ 0.30–0.40):

| estimator | run 17 | run 19 | run 20 | run 21 | spread |
|---|---|---|---|---|---|
| max bin (MPV) | 303.6 | 308.3 | 284.8 | 341.3 | **20 %** |
| median | 305.1 | 302.9 | 301.4 | 309.4 | 2.6 % |
| interquartile truncated mean | 307.7 | 299.8 | 302.6 | 306.0 | 2.6 % |

The real drift over those adjacent runs is ~2.6 %. The 20 % is pure estimator noise: the Landau
peak is broad and flat, so on one run's statistics the highest bin wanders and a 3-point parabola
happily refines the wrong bin. Fed into `f(run)` that becomes a 20 % per-run rescaling — it would
smear the very bands the gain match exists to sharpen, while looking like a sound measurement.
The default is therefore the **interquartile truncated mean** (`estimator=0`); median, peak bin
and Landau fit are available for comparison.

Two more things the macro does deliberately: it ranges the dE/dx axis on **3× the median** rather
than a high quantile (a 97 % quantile puts the upper edge ~4.5× the peak and the bins then
straddle it, and that coarseness lands straight in the factor), and it **interpolates** thin runs
from measured neighbours with a `interp`/`held` label rather than silently assigning 1.0, which
would leave a run unmatched inside a matched set and be invisible downstream.

### Result on the 75 D2 runs

`gainmatch_C15d.csv`, 67 measured + 8 interpolated:

```
f = 0.874 (run 17)  ->  1.128 (run 100)      29 % monotonic drift
mean |Δf| between consecutive runs = 0.0129,  max 0.0433
```

**The validity test passes**: genuine gain drift is smooth in run number, and here the
step-to-step noise is ~20× smaller than the total drift. Scatter would have meant the anchor was
measuring estimator noise instead.

Applying it narrows the bands where a single band dominates — Bρ 0.25–0.30 by 11.1 %, 0.30–0.40
by 7.7 %, 0.40–0.50 by 6.7 % — and does nothing at Bρ 0.15–0.20 where bands overlap.

### The gain is applied at read time

The 72 GB of reco holds **raw** dE/dx. `gain_C15d.h` loads the table and `mkpid_C15d.C` /
`apply_gate_C15d.C` scale on the fly; `AtGainMatchTask` does the same inside a framework chain.
The correction therefore stays re-derivable rather than frozen into the reco.

⚠ **Anything that reads the caches and cuts on dE/dx must apply the same table.** A gate drawn on
a matched plane and applied to raw values silently selects the wrong tracks. Both macros warn when
a run has no table entry rather than passing it through as 1.0 unnoticed.

## Gates

### Browser explorer (no X11)

```bash
root -b -q 'make_pid_explorer_html.C()'      # -> ~/C15d_pid_explorer.html, ~6.7 MB
```

One self-contained HTML file with the plane baked in — same idea as the a2091 kinematics explorer
(`a2091/UKF/pp/make_explorer_html.C`), but showing the PID plane, because that is what exists at
this stage and drawing gates on it is what blocks the analysis. No server, no X11, no CDN.

Click **draw**, click round a band, and close the polygon (double-click, or click the first
vertex). The page reports the in-gate count live and emits **spyral_utils Cut2D JSON** with Z/A —
byte-compatible with `AtCut2D`, `draw_gate_C15d.C` and `apply_gate_C15d.C`. Save it as
`gates/<name>.json`.

Live controls: run range, `nClusters`, vertex R, axis ranges, binning, log z, and a
**matched/raw toggle**. Raw dE/dx is baked in with the gain table beside it and the factor applied
in the page, so the toggle is live and the underlying measurement is always present — the page
cannot show a matched plane it cannot reproduce, and it says so in red when no table is loaded.

It also plots the **per-run in-gate fraction** for the current gate, and flags a spread > 1.5× as a
gain-matching problem rather than physics.

Two things to know. The default axes (√dE/dx ≤ 60, Bρ ≤ 2) hold **87.5%** of tracks; Bρ carries a
1/sin(polar) divergence for near-axis tracks, so 2.6% sit above 5 T·m — raise the y-max slider to
see them. And the values are baked as scaled integers (√dE/dx ×100, Bρ ×1000), which is what keeps
400k tracks to ~6.7 MB.

### In ROOT

`draw_gate_C15d.C` (interactive — run without `-b`), applied with
`apply_gate_C15d.C`, which writes a `sel` tree of `(run, event, trackID)` so the selection can be
used by the fitting stage without recomputing `AtSpyralPID`. A gate you can only draw is half a
gate.

Do not import a gate from another pipeline. Bρ reproduces between pipelines to well under a
percent, but dE/dx does not — it is charge over arclength, and charge is not measured the same
way. A foreign gate lands somewhere on this plane and returns a plausible-looking count for the
wrong species.

`apply_gate_C15d.C` reports the **per-run** in-gate fraction: a gate that holds over part of the
run set and not the rest is a gain-matching failure, and it hides inside a healthy total.

### Band identification — checked, not assumed

At fixed Bρ every Z=1 species carries the same momentum, so the heavier ones are slower and sit
higher in dE/dx by 1/β. That predicts the band spacing with no free parameter but one overall
scale:

| Bρ slice | N | observed ratios | predicted (d, t) |
|---|---|---|---|
| 0.70–0.90 | 1982 | 1.95, 3.00 | 1.95, 2.91 |
| 0.90–1.20 | 2479 | 1.93 | 1.92 |

So the dominant lowest band is **protons**, then **deuterons**, then **tritons**, and Bρ and
dE/dx are internally consistent. (Leading-order 1/β only; the low-statistics slices are messier,
and Bρ 0.45–0.55 shows a fourth peak at ratio 4.23 where a Z=2 alpha would predict 3.95 —
suggestive but not established.)

## Where it stands

| | |
|---|---|
| IC beam window | `pid/ic_C15d.json`, **931–1413 ADC** (the 1168 peak), edges at 0.1% of peak height |
| PID gate | `pid/proton_C15d.json`, 27 vertices, Z=1 A=1 |
| Selection | 220,229 IC-gated tracks → **132,815 in gate (60.3%)** → `pid/sel_proton_C15d.root` |
| Per-run yield | 58.5%–63.7% — flat, so the gain match is holding |

`apply_gate_C15d.C` reads `pid/points_C15d.root`, which already carries the gain match and the IC
join with its length checks — deriving either again in the consumer would be two implementations
that can disagree. **Pass the same IC window the gate was drawn with**: a gate drawn on a
single-beam plane and applied to the full cocktail counts tracks from a beam it never meant to
select, and the number looks perfectly reasonable.

## Open

- Reaction masses are not yet fixed anywhere in the workspace — the fitters only need the
  ejectile. The excitation-energy step will need beam, target and residual.
- No IC gate and no PID gate yet. Build both from the persisted (gain-matched) `AtPIDEvent`.
- `run_0047` has no FRIB data, so it cannot receive an IC gate.
