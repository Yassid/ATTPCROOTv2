#!/usr/bin/env python
"""Local direction-continuity clustering (offline, no training, shape-aware).

Per hit: local tangent from PCA of its k 3D neighbours. Connect two neighbouring hits
only if BOTH tangents align with the segment between them AND with each other (the pair
is locally collinear = same trajectory), AND their charge (dE/dx) is continuous. Then
connected components. Cuts at crossings (tangents differ) and proton/17C junctions
(dE/dx jumps); follows each curve. Conservative thresholds -> over-segment (stitchable).

  ~/gnn_env/bin/python dircluster.py               # gallery on real data
  ~/gnn_env/bin/python dircluster.py 2147,1198,1718,2272   # specific real events
"""
import sys, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# ---- parameters ----
K          = 12      # neighbours for tangent PCA + candidate edges
COS_SEG    = 0.78    # ~39 deg
COS_TAN    = 0.72    # ~44 deg
RMAX       = 32.0    # mm
QRATIO     = 0.65    # smoothed dE/dx; balances merge ~2% with proton ~95% + 17C ~50% recovery
MIN_HITS   = 3       # keep small fragments (stitchable) instead of dropping

def cluster(P, Q, k=K, cos_seg=COS_SEG, cos_tan=COS_TAN, rmax=RMAX, qratio=QRATIO, min_hits=MIN_HITS,
            smooth_q=True):
    n = len(P)
    if n < k + 1:
        return -np.ones(n, dtype=int)
    nn = NearestNeighbors(n_neighbors=k+1).fit(P)
    _, idx = nn.kneighbors(P)
    nb = idx[:, 1:]                                   # (n,k) neighbour indices (drop self)
    # local tangent per hit = principal axis of neighbourhood
    neigh = P[nb]                                     # (n,k,3)
    d = neigh - neigh.mean(1, keepdims=True)
    cov = np.einsum('nki,nkj->nij', d, d) / k
    w, v = np.linalg.eigh(cov)                        # ascending; largest = col 2
    tang = v[:, :, 2]                                 # (n,3) unit tangent
    # local dE/dx: median neighbourhood charge rides over the heavy-ion Bragg profile (keeps
    # the 17C recoil whole) while still separating the ~5-8x proton/17C gap at the vertex.
    Qeff = np.median(Q[nb], axis=1) if smooth_q else Q
    # candidate edges from kNN
    src = np.repeat(np.arange(n), k)
    dst = nb.ravel()
    seg = P[dst] - P[src]
    dl = np.linalg.norm(seg, axis=1) + 1e-9
    su = seg / dl[:, None]
    ai = np.abs((su * tang[src]).sum(1))             # segment vs tangent_i
    aj = np.abs((su * tang[dst]).sum(1))             # segment vs tangent_j
    at = np.abs((tang[src] * tang[dst]).sum(1))      # tangent_i vs tangent_j
    qi, qj = Qeff[src], Qeff[dst]
    qr = np.minimum(qi, qj) / (np.maximum(qi, qj) + 1e-9)
    keep = (ai > cos_seg) & (aj > cos_seg) & (at > cos_tan) & (dl < rmax) & (qr > qratio)
    s, t = src[keep], dst[keep]
    A = coo_matrix((np.ones(len(s)), (s, t)), shape=(n, n))
    A = A + A.T
    ncomp, lab = connected_components(A, directed=False)
    # relabel: drop tiny components to noise (-1)
    out = -np.ones(n, dtype=int); nxt = 0
    for c in range(ncomp):
        m = lab == c
        if m.sum() >= min_hits:
            out[m] = nxt; nxt += 1
    return out

def draw(ax, x, y, lab, title):
    ax.scatter(x[lab < 0], y[lab < 0], s=3, c='lightgray')
    cols = cm.tab20.colors
    for i, cl in enumerate(sorted(set(lab[lab >= 0]))):
        m = lab == cl; ax.scatter(x[m], y[m], s=7, color=cols[i % 20])
    ax.set_title(title, fontsize=8); ax.tick_params(labelsize=6)

if __name__ == "__main__":
    D = pd.read_csv("labeling/data/real_events.csv")
    if len(sys.argv) > 1 and any(ch.isdigit() for ch in sys.argv[1]):
        evs = [int(e) for e in sys.argv[1].split(",")]; out = "diagnostics/dircluster_events.png"
    else:
        sizes = D.groupby('event').size()
        cand = np.array(sizes[(sizes > 100) & (sizes < 1200)].index.to_numpy(), copy=True)
        np.random.seed(3); np.random.shuffle(cand); evs = list(cand[:16]); out = "diagnostics/dircluster_gallery.png"
    ncol = 2 if len(evs) <= 6 else 4
    nrow = int(np.ceil(len(evs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4*ncol, 3.1*nrow)); axes = np.atleast_1d(axes).flatten()
    for i, ev in enumerate(evs):
        g = D[D.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        lab = cluster(np.stack([x, y, z], 1), q)
        ncl = len(set(lab[lab >= 0])); ng = int((lab < 0).sum())
        draw(axes[i], x, y, lab, f"evt {ev}: {ncl} clusters, {ng} gray, {len(x)} hits")
    for j in range(len(evs), len(axes)): axes[j].axis('off')
    fig.suptitle("Local direction-continuity clustering (real data)", fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=95); print("wrote", out)
