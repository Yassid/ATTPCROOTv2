#!/usr/bin/env python
# Compare ATTPCROOT-GNN (embed_attpc, trained on triplclust labels) vs Spyral HDBSCAN on the
# HELD-OUT run (run_0305). Prints a summary table + saves a side-by-side gallery.
#   Run:  ~/gnn_env/bin/python run_comparison.py
import sys, json
import pandas as pd, numpy as np
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RUN = 305
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM = EMB_DIM = 4; K = 20; EPS = 1.1

def mlp(ch):
    L = []
    for i in range(len(ch) - 1):
        L += [nn.Linear(ch[i], ch[i + 1]), nn.BatchNorm1d(ch[i + 1]), nn.ReLU()]
    return nn.Sequential(*L)

class Net(nn.Module):
    def __init__(s, k=20):
        super().__init__()
        s.ec1 = DynamicEdgeConv(mlp([2 * IN_DIM, 64, 64]), k, 'max')
        s.ec2 = DynamicEdgeConv(mlp([2 * 64, 64, 64]), k, 'max')
        s.ec3 = DynamicEdgeConv(mlp([2 * 64, 128]), k, 'max')
        s.head = nn.Sequential(nn.Linear(256, 128), nn.BatchNorm1d(128), nn.ReLU(), nn.Linear(128, EMB_DIM))
    def forward(s, x, b):
        x1 = s.ec1(x, b); x2 = s.ec2(x1, b); x3 = s.ec3(x2, b)
        return s.head(torch.cat([x1, x2, x3], 1))

def clean_proton_rate(df, labelcol):
    # heuristic: an event "has a clean proton" if it has a cluster that is off-axis, elongated and
    # of physical size (25-200 hits) -- a separated ejectile, not the beam blob.
    good = 0; tot = 0
    for e, g in df.groupby('event'):
        tot += 1
        for cl in g[g[labelcol] >= 0][labelcol].unique():
            c = g[g[labelcol] == cl]
            n = len(c)
            if not (25 <= n <= 200): continue
            r = np.hypot(c.x - c.x.mean(), c.y - c.y.mean())
            zspan = c.z.max() - c.z.min()
            if r.max() > 40 and zspan > 150:  # off-axis + extends in z
                good += 1; break
    return 100.0 * good / max(tot, 1)

def main():
    # --- GNN inference on ATTPCROOT clouds (held-out run) ---
    a = pd.read_parquet("data/attpc_all.parquet"); a = a[a.run == RUN].copy()
    m = Net(K).to(DEV); m.load_state_dict(torch.load("sup/embed_attpc.pt", map_location=DEV)); m.eval()
    gn = []
    for gid, g in a.groupby('gid'):
        x = g.x.values; y = g.y.values; z = g.z.values; q = g.q.values; n = len(x)
        lab = np.full(n, -1)
        if n >= 10:
            feat = np.stack([x / 250, y / 250, (z - 600) / 450, np.log1p(np.clip(q, 0, None)) / 9.0], 1)
            with torch.no_grad():
                emb = m(torch.tensor(feat, dtype=torch.float, device=DEV),
                        torch.zeros(n, dtype=torch.long, device=DEV)).cpu().numpy()
            lab = DBSCAN(eps=EPS, min_samples=6).fit_predict(emb)
        gg = g.copy(); gg['gnn'] = lab; gn.append(gg)
    gn = pd.concat(gn)

    sp = pd.read_parquet("data/spyral_all.parquet"); sp = sp[sp.run == RUN].copy()

    def stats(df, col):
        nev = df.gid.nunique() if 'gid' in df else df.event.nunique()
        ncl = df[df[col] >= 0].groupby(['gid' if 'gid' in df else 'event', col]).ngroups
        key = 'gid' if 'gid' in df else 'event'
        cpe = df[df[col] >= 0].groupby(key)[col].nunique()
        sz = df[df[col] >= 0].groupby([key, col]).size()
        multi = (cpe >= 2).mean() * 100
        return nev, cpe.median(), sz.median(), multi

    g_nev, g_cpe, g_sz, g_multi = stats(gn, 'gnn')
    s_nev, s_cpe, s_sz, s_multi = stats(sp, 'label')
    gn_cp = clean_proton_rate(gn, 'gnn')  # gn already has per-run 'event' col (unique within run_0305)
    sp_cp = clean_proton_rate(sp, 'label')

    print(f"\n===== ATTPCROOT-GNN  vs  Spyral   (held-out run_0{RUN}) =====")
    print(f"{'metric':32s} {'ATTPCROOT-GNN':>14s} {'Spyral':>10s}")
    print(f"{'events clustered':32s} {g_nev:>14d} {s_nev:>10d}")
    print(f"{'median clusters / event':32s} {int(g_cpe):>14d} {int(s_cpe):>10d}")
    print(f"{'median hits / cluster':32s} {int(g_sz):>14d} {int(s_sz):>10d}")
    print(f"{'events with >=2 clusters (%)':32s} {g_multi:>13.0f}% {s_multi:>9.0f}%")
    print(f"{'clean-proton rate (%)':32s} {gn_cp:>13.0f}% {sp_cp:>9.0f}%")

    # --- style gallery: NOT event-matched (chunk numbering doesn't map to Spyral global #s).
    # Left column: ATTPCROOT-GNN events; right column: Spyral events. Shows clustering STYLE. ---
    PAL = plt.cm.tab10.colors; col = lambda c: (0.7, 0.7, 0.7, 0.3) if c < 0 else PAL[int(c) % 10]
    g_multi_ev = [gid for gid, g in gn.groupby('gid') if 80 < len(g) < 600][20:26]
    s_multi_ev = [e for e, g in sp.groupby('event') if 80 < len(g) < 600][20:26]
    fig, ax = plt.subplots(6, 2, figsize=(13, 19))
    for r in range(6):
        gg = gn[gn.gid == g_multi_ev[r]]; gs = sp[sp.event == s_multi_ev[r]]
        for cl in sorted(gg.gnn.unique()):
            mk = gg.gnn == cl; ax[r, 0].scatter(gg.z[mk], gg.y[mk], s=7, color=col(cl))
        for cl in sorted(gs.label.unique()):
            mk = gs.label == cl; ax[r, 1].scatter(gs.z[mk], gs.y[mk], s=7, color=col(cl))
        ax[r, 0].set_title(f"ATTPCROOT-GNN  ({gg[gg.gnn>=0].gnn.nunique()} clusters)", fontsize=9)
        ax[r, 1].set_title(f"Spyral  ({gs.label.nunique()} clusters)", fontsize=9)
        for c in range(2):
            ax[r, c].set_xlabel('z (beam axis)'); ax[r, c].set_ylabel('y')
    fig.suptitle("Clustering STYLE comparison (independent events, not matched)", y=1.001)
    plt.tight_layout(); plt.savefig("comparison_gnn_vs_spyral.png", dpi=80)
    print("\nsaved style gallery -> comparison_gnn_vs_spyral.png (independent events, not matched)")

if __name__ == "__main__":
    main()
