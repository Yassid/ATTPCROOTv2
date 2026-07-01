# 16C(d,p)17C labeled simulation — GNN training data

Generates **labeled point clouds** (per-hit truth track IDs) for training a learned
clustering / track-finding model (see `project_gnn_sim_tracking`). Matches a1975
deuterium-target conditions: D2 gas, B = 2.85 T, beam ~12 MeV/u.

## Pipeline
```
run_batch.sh N            # chains the three steps below for N sim events
```
1. **`C16_dp_sim.C(N)`** — Geant4 sim of 16C(d,p)17C. Beam 16C, products 17C + p,
   ThetaCMS 2–178°. Geometry `ATTPC_D300torr_v2.root` (D2 gas), B = 28.5 kG.
   → `data/attpcsim.root` (AtMCPoint truth).
2. **`run_digi_16C.C`** — Clusterize → Pulse(`SetSaveMCInfo`) → PSA. Pad map
   `ANL2023.xml`, par `ATTPC.a1975_deuterium.par`. → `data/output_digi.root`
   (AtRawEvent w/ pad→mcPointID map, AtEventH hits, AtTpcPoint_1 truth).
3. **`extract_labels.C`** — joins hit.padNum → SimMCPointMap → mcPointID →
   AtMCPoint.trackID/A/Z. → `data/labeled.csv` (100% of hits labeled).
4. **`package_training.py`** — → `data/train.parquet` (GNN-ready, dense per-event
   cluster `label`). Keeps reaction events (proton present).

## Truth label = composite (trackID, A, Z)
The beam-propagator gives the 16C beam AND the 17C recoil the **same Geant4 trackID 0**
(species changes at the vertex), so the cluster label must combine trackID with (A,Z).
Particles: `16C_beam`, `17C_recoil`, `proton`, occasional secondaries (α, d).

## Current set (run_batch.sh 3000)
`data/train.parquet`: **1149 reaction events, 50877 hits**. 1103 two-track
(proton+17C, vertex-shared), 60 spiral protons (>150 hits). Plots: `data/train_gallery.png`.

## Backward-lab focus (FLAGGED)
Current focus = backward protons (lab 90–180°). For this inverse-kinematics reaction
**forward CM → backward lab** (CM[2,60]→lab 76–168°, 61% backward), so `C16_dp_sim.C`
defaults to **CM[2,55]**. The definitive cut is `select_lab_angle()` in
`package_training.py` (`BACKWARD_ONLY=True, LAB_MIN/MAX=90/180`). Revert to all angles:
`C16_dp_sim(N,"TGeant4",2,178)` + `BACKWARD_ONLY=False`.

## Data-driven beam+noise overlay (`build_overlay.py`)
Beam particles produce a lot of spread-out noise that **scales with how far the beam
travels before reacting**. So we don't simulate the beam — we sample it from data:
1. `dump_real.C` → real event clouds. Select **beam-only canvases** (PRA-free: reject
   events with a dense off-axis cluster). 271/600 run_0016 events qualify.
2. Translate each sim event's vertex (proton max-|p| hit) to a sampled depth on the beam
   axis, then **truncate the canvas at the vertex** (keep the beam path upstream).
   → noise ∝ beam length. Verified: **corr(vertex_z, n_noise) = −0.55**.
3. Backward protons go toward the entrance → naturally overlap the beam (the hard case).

Output `data/overlay.parquet`: labels **0=proton, 1=17C, 2=beam/noise**.
Current set: **1457 events, 299k hits** (~47% signal / 53% real beam+noise), 518 spirals.
Gallery `data/overlay_backward_gallery.png`. All cuts are flagged params at the top of
`build_overlay.py` (vertex z range, beam direction, canvas density, detector z bounds).

## Known limitations / TODO for realism
- **No beam track in reaction events.** Sim alternates beam-only (even) / reaction
  (odd) events, and the beam does NOT connect to the reaction vertex (deposits near
  the entrance z~950 vs vertex z~150-400) — so they can't be merged. Add the beam as a
  synthetic straight axis track, or fix the generator/propagator to co-simulate.
- **No detector noise** (real clouds have ~beam + electronic noise). Add before
  training a model that must run on real a1975 data.
- **Sim→data domain gap** — validate any sim-trained model on real events (cf. the
  AT-TPC point-cloud detector-response translation work).
- Many protons are short (median ~10 hits, wide CM range). Narrow ThetaCMS / tune beam
  energy if more substantial tracks are wanted.

## Event-by-event scanner (ROOT macro)
```
~/Spyral/venv/bin/python convert_overlay_root.py   # parquet -> data/overlay.root (flat TTree)
source $REPO/build/config.sh
root -l 'view_overlay.C(0)'
```
On-canvas buttons: |< First / < Prev / Next > / Last >| / "Next spiral >" (proton >150 hits).
Also from the prompt: vnext() vprev() vfirst() vlast() vnext_spiral() view_overlay(N).
Colors: blue=proton, red=17C, gray=beam/noise. XY pad-plane + ZY drift side views.
(scale_up.sh refreshes overlay.root automatically.)
