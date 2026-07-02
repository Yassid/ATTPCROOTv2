#!/usr/bin/env python
"""Connect-the-dots track builder: directional track-following (cellular-automaton style).

Seed from the most 'linear' unused hit, then GROW in both directions by repeatedly picking
the neighbour that best continues the current momentum (straightest, within a cone + step),
gated by dE/dx continuity. Direction is updated each step so it FOLLOWS curvature (spirals).
Greedy used-marking => at a crossing the straight continuation wins and the other track's hits
are left for a later seed. Over-segments at genuine breaks (stitchable). No training, fast.

  ~/gnn_env/bin/python trackfollow.py                    # gallery
  ~/gnn_env/bin/python trackfollow.py 2147,1198,1718,2272 # specific events
"""
import sys, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors

# ---- parameters ----
K        = 12
R_STEP   = 18.0    # mm, max step to next hit
COS_CONE = 0.60    # ~53 deg forward cone (relative to current momentum)
QRATIO   = 0.45    # dE/dx continuity gate (smoothed local charge)
ALPHA    = 0.55    # momentum smoothing (higher = stiffer; lower = turns faster)
MIN_HITS = 6

def follow(P, Q, k=K, r_step=R_STEP, cos_cone=COS_CONE, qratio=QRATIO, alpha=ALPHA, min_hits=MIN_HITS):
    n = len(P)
    lab = -np.ones(n, dtype=int)
    if n < k + 1:
        return lab
    nn = NearestNeighbors(n_neighbors=k+1).fit(P)
    _, idx = nn.kneighbors(P)
    nb = idx[:, 1:]
    d = P[nb] - P[nb].mean(1, keepdims=True)
    cov = np.einsum('nki,nkj->nij', d, d) / k
    w, v = np.linalg.eigh(cov)
    tang = v[:, :, 2]
    lin = w[:, 2] / (w.sum(1) + 1e-9)               # linearity (seed order)
    Qloc = np.median(Q[nb], axis=1)

    used = np.zeros(n, bool)
    order = np.argsort(-lin)
    nxt = 0
    for seed in order:
        if used[seed]:
            continue
        track = [seed]; used[seed] = True
        for sign in (1.0, -1.0):
            cur = seed
            dirv = tang[seed] * sign
            cq = Qloc[seed]
            while True:
                cand = nb[cur]
                vv = P[cand] - P[cur]
                L = np.linalg.norm(vv, axis=1)
                with np.errstate(invalid='ignore', divide='ignore'):
                    dots = (vv @ dirv) / (L + 1e-9)
                    qok = np.minimum(cq, Qloc[cand]) / np.maximum(cq, Qloc[cand] + 1e-9) > qratio
                valid = (~used[cand]) & (L <= r_step) & (L > 1e-6) & (dots > cos_cone) & qok
                if not valid.any():
                    break
                scores = np.where(valid, dots, -2.0)
                j = int(scores.argmax())
                best = int(cand[j])
                track.append(best); used[best] = True
                step = vv[j] / (L[j] + 1e-9)
                dirv = alpha * dirv + (1 - alpha) * step
                dirv = dirv / (np.linalg.norm(dirv) + 1e-9)
                cq = alpha * cq + (1 - alpha) * Qloc[best]
                cur = best
        if len(track) >= min_hits:
            lab[track] = nxt; nxt += 1
    return lab

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    D = pd.read_csv("labeling/data/real_events.csv")
    if len(sys.argv) > 1 and any(c.isdigit() for c in sys.argv[1]):
        evs = [int(e) for e in sys.argv[1].split(",")]; out = "diagnostics/trackfollow_events.png"
    else:
        sz = D.groupby('event').size()
        cand = np.array(sz[(sz > 100) & (sz < 1200)].index.to_numpy(), copy=True)
        np.random.seed(3); np.random.shuffle(cand); evs = list(cand[:16]); out = "diagnostics/trackfollow_gallery.png"
    ncol = 2 if len(evs) <= 6 else 4
    nrow = int(np.ceil(len(evs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4*ncol, 3.1*nrow)); axes = np.atleast_1d(axes).flatten()
    for i, ev in enumerate(evs):
        g = D[D.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        lab = follow(np.stack([x, y, z], 1), q)
        ncl = len(set(lab[lab >= 0])); ng = int((lab < 0).sum())
        axes[i].scatter(x[lab < 0], y[lab < 0], s=3, c='lightgray')
        cols = cm.tab20.colors
        for j2, cl in enumerate(sorted(set(lab[lab >= 0]))):
            m = lab == cl; axes[i].scatter(x[m], y[m], s=7, color=cols[j2 % 20])
        axes[i].set_title(f"evt {ev}: {ncl} tracks, {ng} gray, {len(x)} hits", fontsize=8)
        axes[i].tick_params(labelsize=6)
    for j in range(len(evs), len(axes)): axes[j].axis('off')
    fig.suptitle("Directional track-following (connect-the-dots) on real data", fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=95); print("wrote", out)
