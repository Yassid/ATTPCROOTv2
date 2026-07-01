# 16C(d,p) GNN training-data diagnostics (2026-07-01, autonomous run)

Goal: answer "was the training data good enough?" before retraining. Compared
**real a1975 (d,p)** clouds (run_0016 `AtEventH`, AtPSAMultiFit) against **regenerated
sim signal** (C16_dp_sim 3000 → run_digi_16C(600) → extract_labels). All scripts +
plots in this folder; data regenerable.

## 1. Labels are fine — realism is the problem
- Sim labels: **100% of hits labeled, 0% ambiguous** (`shared` flag). Not a label problem.
- Dominant domain gap = **per-hit charge**: real median 298 vs sim 36 (~8×), **KS=0.76**.
  (`dist_compare.png`) This alone plausibly explains the prior 73%→13% sim→real
  proton-recovery collapse — the model trained on sim charge, met 8×-different real charge.

## 2. Charge gap is fixable for free — per-event normalization (`chargenorm_test.png`)
| charge feature | KS (real vs sim) |
|---|---|
| raw `q` | 0.76 |
| `log10(q)` | 0.76 (no help — monotonic, can't remove an offset) |
| **per-event z-score of log q** | **0.11** |
| **per-event quantile rank** | **0.005** (≈identical) |

The 8× gap is pure per-event gain/scale, NOT relative charge structure. **Recipe:**
charge → per-event quantile rank (bounded [0,1], outlier-robust); positions x,y,z →
keep physical / global-standardize (NOT per-event — absolute geometry matters).

## 3. Composition gap — the overlay does NOT close it (`overlay_verify.png`)
Sim reaction events lack a beam track + electronic noise. The existing `build_overlay.py`
was supposed to fix this by pasting real beam-only canvases. **It fails:**
- beam-only canvases are nearly empty: **median 20 hits** (vs reaction-like real 245)
- truncation discards ~70% → only ~6 noise hits/event
- vertex translation clips ~39% of sim signal (backward protons pushed past z bound)
- → overlay ~89 hits/event, can't reach real reaction events (~245). KS(nhits) 0.40→0.51 (worse).

**Charge-consistency check passed:** without the charge fix, per-event normalization
trivially separates sim(label 0/1) from real-noise(label 2) — a domain giveaway
(KS 0.43). With quantile-matched charge (`build_overlay_cc.py`) it drops to **0.14**. ✅

## 4. Composition gap CAN be closed — overlay v2 (`overlay_v2_verify.png`)
`build_overlay_v2.py` adds: **synthetic beam track** (dense axial line entrance→vertex,
README TODO #1) + reduced clipping (39%→16%) + charge consistency + sparse real-noise haze.

| nhits KS | sim-only | overlay-cc | **overlay-v2** |
|---|---|---|---|
| vs real-rxn (245 hits) | 0.40 | 0.51 | **0.16** ✅ |
| vs real-all (120 hits) | 0.18 | 0.21 | 0.47 |

## 5. The irreducible limitation → need labeled real events
v2 matches **real-rxn** (busy events, 245 hits); cc matches **real-all** (120 hits).
**Which is the true target for a real (d,p) reaction event is unknown without
hand-labeled real events.** Marginal-matching is underdetermined. So:
- charge normalization: solved, bake into the feature builder.
- composition: a synthetic beam closes it, but tuning the target needs a small
  hand-labeled real (d,p) gold set — which is also the only honest validation set.

## Next steps (recommended)
1. Bake per-event quantile charge + global-standardized positions into the feature builder.
2. Build the manual labeling tool (matplotlib lasso on XY/ZY, HDBSCAN pre-seed) for a
   small real (d,p) gold set — defines the composition target AND validates the model.
3. Settle synthetic-beam params (spacing/charge band) against that gold set, then retrain.

## Files
`dumpRealH.C` (real cloud dump), `compare_distributions.py`, `test_chargenorm.py`,
`build_overlay_cc.py`, `build_overlay_v2.py`, `verify_overlay.py`, `verify_v2.py`,
`overlay_v2_sample.parquet` (sample output). Sim regen: C16_dp_sim.C → run_digi_16C.C →
extract_labels.C (all in parent dir). Python env: `~/Spyral/venv` (polars), `~/gnn_env`
(pandas/scipy). Backup of prior training: `~/16C_dp_gnn_training_backup_20260630.tar.gz`.
