#!/usr/bin/env python
"""Select real REACTION-candidate events for hand labeling.

From a real event dump (data/real_events.csv: event,x,y,z,q) keep only events that
look like a reaction (a dense off-axis track present), so the labeler spends time on
useful events, not empty beam triggers. Writes data/label_input.parquet ordered by a
"labelability" score (off-axis hits first) so the clearest events come up first.

Run:  ~/gnn_env/bin/python prepare_label_data.py [--max N] [--in CSV]
Regenerate other runs first with dumpRealH.C (see ../diagnostics/dumpRealH.C).
"""
import argparse, os, numpy as np, pandas as pd
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("--in",  dest="inp", default="data/real_events.csv")
ap.add_argument("--out", default="data/label_input.parquet")
ap.add_argument("--max", type=int, default=150, help="max events to keep for labeling")
ap.add_argument("--r-offaxis",  type=float, default=35.0)
ap.add_argument("--dens-r",     type=float, default=12.0)
ap.add_argument("--dens-thresh", type=int,  default=7)
ap.add_argument("--nhits-min", type=int, default=40,  help="drop near-empty events")
ap.add_argument("--nhits-max", type=int, default=500, help="drop pileup/haze monsters")
ap.add_argument("--keep-existing", action="store_true",
                help="preserve events already in --srcids at their current label ids (0..K-1), "
                     "then append NEW candidates up to --max. Keeps existing labels.parquet valid.")
ap.add_argument("--srcids", default="data/label_input_srcids.csv",
                help="label_event->src_event map to preserve when --keep-existing")
a = ap.parse_args()

df = pd.read_csv(a.inp)
df["r"] = np.hypot(df.x, df.y)
scored = []
for ev, g in df.groupby("event"):
    if not (a.nhits_min <= len(g) <= a.nhits_max):
        continue                      # keep clearly-labelable events (not empty, not haze-bombed)
    off = g[g.r >= a.r_offaxis]
    if len(off) < 5:
        continue
    Q = off[["x", "y", "z"]].to_numpy()
    dens = cKDTree(Q).query_ball_point(Q, r=a.dens_r, return_length=True)
    if dens.max() < a.dens_thresh:
        continue                      # no dense off-axis track -> not a reaction candidate
    scored.append((int(ev), len(g)))

# simplest events first (fewest hits) so labeling warms up easy, then take a spread up to --max
scored.sort(key=lambda t: t[1])

if a.keep_existing and os.path.exists(a.srcids):
    # PRESERVE the events already selected (in their exact label order 0..K-1) so any hand
    # labels in data/labels.parquet stay aligned; APPEND fresh candidates as new ids K..N-1.
    ex = pd.read_csv(a.srcids).sort_values("label_event")
    existing = [int(s) for s in ex.src_event.tolist()]
    existing = [ev for ev in existing if (df.event == ev).any()]   # must still be in the dump
    existing_set = set(existing)
    new_pool = [ev for ev, _ in scored if ev not in existing_set]  # simplest-first, not-yet-used
    n_new = max(0, a.max - len(existing))
    if len(new_pool) > n_new:
        idx = np.linspace(0, len(new_pool) - 1, n_new).round().astype(int)
        new_pool = [new_pool[i] for i in idx]                      # even spread over the new ones
    keep_ids = existing + new_pool
    print(f"keep-existing: {len(existing)} preserved (labels stay valid) + {len(new_pool)} new "
          f"= {len(keep_ids)} total")
else:
    if len(scored) > a.max:
        idx = np.linspace(0, len(scored) - 1, a.max).round().astype(int)
        scored = [scored[i] for i in idx]     # even spread across the size range
    keep_ids = [ev for ev, _ in scored]
out = pd.concat([df[df.event == ev] for ev in keep_ids])[["event", "x", "y", "z", "q"]]
# re-index events 0..K-1 in labeling order (stable, tool-friendly)
remap = {ev: i for i, ev in enumerate(keep_ids)}
out = out.assign(event=out.event.map(remap)).sort_values(["event"]).reset_index(drop=True)
out.to_parquet(a.out, index=False)
print(f"reaction candidates: {len(scored)} found, kept {len(keep_ids)} -> {a.out}")
print(f"  {out.shape[0]} hits, median {out.groupby('event').size().median():.0f} hits/event")
print(f"  (original event ids preserved in data/label_input_srcids.csv)")
pd.DataFrame({"label_event": range(len(keep_ids)), "src_event": keep_ids}
             ).to_csv("data/label_input_srcids.csv", index=False)
