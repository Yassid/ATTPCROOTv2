# 14C(p,p') at 2.85 / 4 / 7 T with AT-TPC and 2 mm pads

A design study: what would the a1954 measurement have looked like in a higher-field magnet with
a finer pad plane. Everything except the field and the pitch is held at the adopted 2026-08-25
a1954 configuration, so a difference between two cells of the matrix is the thing that was
varied and nothing else.

**Results: [RESULTS.md](RESULTS.md).** The short version is that under the adopted analysis no
configuration is detector-limited — the constant beam energy assumed for every vertex sets the Ex
width in all six cells — and once the beam energy is taken at the reconstructed vertex the matrix
becomes meaningful and 4 T + 2 mm comes out ahead of 7 T.

## The matrix

| config       | B [T] | pads   | sims       |
|--------------|-------|--------|------------|
| `b285_attpc` | 2.85  | AT-TPC | `sims_b285` |
| `b285_2mm`   | 2.85  | 2 mm   | `sims_b285` |
| `b400_attpc` | 4.00  | AT-TPC | `sims_b400` |
| `b400_2mm`   | 4.00  | 2 mm   | `sims_b400` |
| `b700_attpc` | 7.00  | AT-TPC | `sims_b700` |
| `b700_2mm`   | 7.00  | 2 mm   | `sims_b700` |

`b285_attpc` is the anchor: it is the a1954 experiment as it was actually run, so the campaign
can be checked against results that already exist rather than only against itself.

Both pitches of a field read the **same** simulation files. Transport depends on the magnetic
field but not at all on the pad plane, so a pitch comparison is free of any generation
difference — and only fifteen generations are needed for thirty reconstructions.

Levels: `gs` (elastic), `ex6094` (1⁻), `ex6728` (3⁻), `ex7012` (2⁺, the state carrying the
B(E2)), `ex8317` (2⁺, isolated above the multiplet). The question the matrix exists to answer is
whether any configuration separates 6.728 from 7.012, which are 284 keV apart.

## Running it

```bash
./run_C14_hf_campaign.sh -j 4          # 8000 events/sample, ~8 h on this box
./hf_status.sh                         # what has actually completed
root -b -q 'compare_hf_C14.C("/mnt/f/a1954_C14_hf")'
```

`run_C14_hf_campaign.sh` runs in two waves: `gs` + `ex6094` across all six configurations first,
so the headline resolution numbers exist before the campaign is half over, then the remaining
three levels. Every stage is resumable on evidence that it finished — a written file is not a
finished job — so an interrupted run costs at most the sample it was inside.

`accumulate_C14_hf.sh <state> <B> <pad_mm> <seed> [nEvents]` is one sample, and is the thing to
call by hand when re-running a single cell.

## What is held fixed

The adopted a1954 recipe, unchanged: RT gas (ρ = 3.308e-5 g/cm³, `ATTPC_H300torr_RT`), beam
14C at 161 MeV, 1 m drift, 3 cm beam hole, `AtPSAMultiFit` + `AtDirDeDxCleaner` +
`AtTrackFinderHDBSCAN(mcs 20, ms 8, ε 10 mm, mover)`, genfit with material effects and native
CATIMA dE/dx, no manual gap eloss, acceptance and resolution evaluated at χ²/ndf < 5 on
`GetKinematicsXtr` with the same truth match the production acceptance uses.

**The pattern-recognition parameters are NOT re-tuned for 2 mm pads.** ε = 10 mm spans one pad
at 8 × 12 mm and five at 2 mm, and a track deposits several times as many hits. Before reading a
pitch difference as a resolution difference, check with `cluster_eval_C14.C` that pattern
recognition is not the thing that changed.

## Diffusion is field-dependent, and that is on purpose

The digitiser turns `CoefT` into `σ_T[mm] = 10·sqrt(2·CoefT·t_drift[µs])`. At the full 1 m drift
(t = 77 µs) the legacy 9e-4 cm²/µs gives σ_T = 3.7 mm — nearly twice a 2 mm pad. Holding it
fixed across the campaign would test fine pads in the one regime where the pitch cannot matter,
and would hide what a strong field does to a TPC besides bending tracks: it magnetises the drift
electrons, and transverse diffusion falls as 1/(1+(ωτ)²).

So `make_hf_pars.sh` takes the field dependence from Magboltz (H₂, 300 torr, 293 K, 50 V/cm,
B ∥ E; `~/attpc_diff_H2`) and **anchors it on the a1954-tuned value**: `CoefT(B) = 9e-4 ·
[D_T(B)/D_T(2.85 T)]²`. The absolute Magboltz coefficient is not substituted for 9e-4, because
9e-4 was tuned so the simulation reproduces the measured cluster widths at 2.85 T and that
agreement is worth more than an ab-initio number. `CoefL` is untouched — a field parallel to E
does nothing to longitudinal diffusion, and the scan includes B = 0 so this is checked rather
than asserted.

`DIFFMODE=fixed ./accumulate_C14_hf.sh ...` selects the `_fixdiff` par set, which freezes
`CoefT` at 9e-4 for every field. That is the control that separates "the field bent the tracks
better" from "the field stopped the electrons spreading"; run it on one or two cells if the
high-field gain needs to be attributed.

## Reading the output

`compare_hf_C14.C` prints, per level and configuration, the median Ex residual, its IQR, the
gaussian-equivalent σ = IQR/1.349, the fraction within ±0.5 MeV, and the overall acceptance —
then a separation table, ΔE/(σ_a+σ_b) for the neighbouring pairs. Above ~2 a pair is resolvable
in a fit, below ~1 it is one bump.

**Median and IQR, never a walk-outward FWHM.** On the 46Ar campaign that estimator returned
0.765 vs 1.800 MeV for two histograms differing only by a constant shift of every entry. These
distributions have a narrow core on broad tails and the walk estimator is noise on them.

The centroid is quoted next to the width because both matter: a configuration that halves the
width and moves the centroid by 0.5 MeV has not resolved anything.

## Files

| file | what |
|---|---|
| `accumulate_C14_hf.sh` | one sample, all stages, resumable |
| `run_C14_hf_campaign.sh` | the whole matrix, in two waves |
| `hf_status.sh` | completed markers per configuration |
| `make_hf_pars.sh` | writes the six par files from the Magboltz scan |
| `ex_res_C14_hf.C` | Ex residual per θ_lab slice; writes a flat tree so nothing needs re-running |
| `compare_hf_C14.C` | per-level tables and the as-analysed multiplet figure |
| `exres_ebeamz_C14.C` | one sample under a vertex-dependent beam energy, and the method floor |
| `summary_hf_C14.C` | the whole matrix in one set of tables, both reconstructions |
| `hf_figures_C14.C` | the multiplet figures, as analysed and corrected |
| `meassigma_check.sh` | is genfit's assumed per-hit error hiding the fine-pitch gain? (no) |

Products live in `/mnt/f/a1954_C14_hf/` (override with `HF_ROOT`).
