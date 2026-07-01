#!/usr/bin/env python
"""Export the reviewed hand-labeled events as the GNN gold set.
Keeps only reviewed=True events, drops the flag. Adds the domain-robust charge feature
(per-event quantile rank, see ../diagnostics/SUMMARY.md) so the gold set is training-ready.
Run: ~/gnn_env/bin/python export_gold.py
"""
import pandas as pd, numpy as np
df = pd.read_parquet("data/labels.parquet")
gold = df[df.reviewed].copy()
if gold.empty:
    raise SystemExit("no reviewed events yet — label some with label_tool.py first")
gold["q_qrank"] = gold.groupby("event")["q"].rank(pct=True)      # domain-robust charge feature
gold = gold[["event", "x", "y", "z", "q", "q_qrank", "label"]]
gold.to_parquet("data/gold.parquet", index=False)
nev = gold.event.nunique()
byc = gold.label.value_counts().to_dict()
has_p = gold[gold.label == 0].event.nunique()
print(f"gold set: {nev} reviewed events, {len(gold)} hits -> data/gold.parquet")
print(f"  hits per class: proton={byc.get(0,0)} 17C={byc.get(1,0)} beam/noise={byc.get(2,0)}")
print(f"  events with a proton track: {has_p}/{nev}")
