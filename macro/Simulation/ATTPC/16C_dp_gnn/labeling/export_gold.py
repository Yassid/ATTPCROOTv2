#!/usr/bin/env python
"""Export reviewed hand-labeled events as the GNN gold set (spiral-instance scheme).
label: 0 = background (beam/noise), 1,2,... = spiral track instances. Adds `is_track`
(binary) and the domain-robust `q_qrank` charge feature (per-event quantile; see
../diagnostics/SUMMARY.md). Run: ~/gnn_env/bin/python export_gold.py
"""
import pandas as pd
df = pd.read_parquet("data/labels.parquet")
gold = df[df.reviewed].copy()
if gold.empty:
    raise SystemExit("no reviewed events yet — label some with label_tool.py first")
gold["is_track"] = (gold.label > 0).astype("int8")
gold["q_qrank"] = gold.groupby("event")["q"].rank(pct=True)
gold = gold[["event", "x", "y", "z", "q", "q_qrank", "label", "is_track"]]
gold.to_parquet("data/gold.parquet", index=False)
nev = gold.event.nunique()
nsp = gold[gold.label > 0].groupby("event").label.nunique()
print(f"gold set: {nev} reviewed events, {len(gold)} hits -> data/gold.parquet")
print(f"  track hits: {int(gold.is_track.sum())}  background hits: {int((gold.is_track==0).sum())}")
print(f"  spirals/event: median {nsp.median():.0f}, max {nsp.max() if len(nsp) else 0}; "
      f"events with >=2 spirals (crossings): {(nsp>=2).sum()}/{nev}")
