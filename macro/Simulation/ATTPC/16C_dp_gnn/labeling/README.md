# Manual labeling tool — real 16C(d,p) GNN gold set

Hand-label real a1975 (d,p) events to (a) validate any GNN on real data and (b)
fine-tune / set the overlay composition target. Label **spirals (curved particle tracks)
as instances** — the GNN's job is track SEPARATION, not species ID (that's genfit/PID
downstream). Labels: **0 = background (beam/noise), 1,2,3… = individual spiral tracks**.
Collapse to binary track/background any time via `is_track`. See `../diagnostics/SUMMARY.md`
for why this gold set is the missing piece.

## Environment
Uses `~/gnn_env` (pandas, scipy, sklearn HDBSCAN, pyarrow) + **PySide6** (installed) for
the interactive Qt backend. GUI works under WSLg (DISPLAY :0).

## 1. Prepare events to label (already run once → `data/label_input.parquet`)
```
~/gnn_env/bin/python prepare_label_data.py           # 150 reaction candidates, 40–500 hits
```
Picks real events with a dense off-axis track (a reaction, not empty beam), size-windowed
so they're actually labelable, simplest-first. Source: `data/real_events.csv` (run_0016
`AtEventH`, from `../diagnostics/dumpRealH.C`). For more/other runs, re-dump then re-prepare.

## 2. Label (interactive)
```
~/gnn_env/bin/python label_tool.py
```
Two projections shown: **X–Y pad plane** and **Z–Y drift**. Point size ∝ charge.
Per event: **lasso a spiral** (it becomes the current track); press **T** to start the next
spiral and lasso it; repeat for crossing tracks. Everything untouched stays background.
Press **N** for the next event. Fast path is 1 lasso/event (single track), a few for crossings.

| key | action | key | action |
|-----|--------|-----|--------|
| *(lasso)* | paint hits = current track id | **G** | toggle HDBSCAN guide colors (reveal tracks) |
| **T / Enter** | start a NEW spiral (next id) | **C** | reset event to all background |
| **1..9** | select a track id (extend/fix) | **D** | toggle 'reviewed' (e.g. pure-beam events) |
| **B / 0** | paint selection = background | **N / →** | next event (autosaves) |
| **S** | save now | **M / ←** | previous event / **Q** quit |

Assigning any spiral auto-marks the event **reviewed**. Autosaves to
`data/labels.parquet` on every navigation and on quit — safe to stop and resume anytime.
Tip: press **G** first to let HDBSCAN color the clusters, spot the proton track, then **G**
off and lasso it.

## 3. Export the gold set
```
~/gnn_env/bin/python export_gold.py                  # reviewed events -> data/gold.parquet
```
Adds the domain-robust `q_qrank` charge feature (per-event quantile). Use `gold.parquet`
as the real validation/fine-tune set for the GNN.

## Files
`prepare_label_data.py`, `label_tool.py`, `export_gold.py`, `data/label_input.parquet`
(events to label), `data/label_input_srcids.csv` (label→source event id map),
`data/labels.parquet` (your work, autosaved), `data/gold.parquet` (export).
