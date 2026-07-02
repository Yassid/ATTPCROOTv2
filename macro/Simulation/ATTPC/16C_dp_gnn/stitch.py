#!/usr/bin/env python
"""Algo 2 of the dual: STITCHER. Recombine over-segmented fragments (from dircluster) that
belong to the same track: their endpoints must continue smoothly (trajectory) AND share a
dE/dx level (same particle). This fixes fragmentation without re-merging crossings, because a
proton fragment won't stitch to a 17C fragment (dE/dx) nor across a kink (direction).

  combined(P,Q) = stitch(dircluster.cluster(...))
  ~/gnn_env/bin/python stitch.py [events]   # gallery of stitched output on real data
"""
import sys, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
from dircluster import cluster


def attach(P, Q, lab, k=12, cos=0.60, qratio=0.40, r_step=22.0, iters=4):
    """Grow confident track cores by absorbing gray hits that continue a labelled
    neighbour's local direction with compatible dE/dx. Iterated so tracks dilate outward.
    Raises the clustered fraction (recovery) without re-merging (dir+dE/dx gate)."""
    lab = lab.copy(); n = len(P)
    nn = NearestNeighbors(n_neighbors=k+1).fit(P); _, idx = nn.kneighbors(P); nb = idx[:, 1:]
    d = P[nb] - P[nb].mean(1, keepdims=True)
    cov = np.einsum('nki,nkj->nij', d, d) / k
    _, v = np.linalg.eigh(cov); tang = v[:, :, 2]
    Qloc = np.median(Q[nb], axis=1)
    for _ in range(iters):
        gray = np.where(lab < 0)[0]
        if len(gray) == 0: break
        newlab = lab.copy(); changed = 0
        for i in gray:
            neigh = nb[i]; ln = neigh[lab[neigh] >= 0]
            if len(ln) == 0: continue
            best = -1; bestscore = cos
            for j in ln:
                vv = P[j] - P[i]; L = np.linalg.norm(vv)
                if L > r_step or L < 1e-6: continue
                al = abs((vv / L) @ tang[j])
                qok = min(Qloc[i], Qloc[j]) / (max(Qloc[i], Qloc[j]) + 1e-9) > qratio
                if al > bestscore and qok:
                    bestscore = al; best = lab[j]
            if best >= 0:
                newlab[i] = best; changed += 1
        lab = newlab
        if changed == 0: break
    return lab

# stitch parameters
GAP_MAX  = 55.0    # mm, max endpoint gap to bridge (bigger than intra-cluster RMAX)
COS_LINK = 0.80    # endpoint continuation + smoothness cone (~37 deg)
QLEVEL   = 0.55    # median-charge ratio between fragments (same dE/dx band)

def frag_features(P, Q, lab):
    feats = {}
    for c in sorted(set(lab[lab >= 0])):
        idx = np.where(lab == c)[0]
        pts = P[idx]; cen = pts.mean(0)
        # principal axis to find the two extremities
        vt = np.linalg.svd(pts - cen, full_matrices=False)[2]
        proj = (pts - cen) @ vt[0]
        ends = [pts[proj.argmin()], pts[proj.argmax()]]
        tans = []
        for e in ends:
            near = pts[np.argsort(np.linalg.norm(pts - e, axis=1))[:min(6, len(pts))]]
            a = np.linalg.svd(near - near.mean(0), full_matrices=False)[2][0]
            if a @ (e - cen) < 0: a = -a          # orient outward
            tans.append(a / (np.linalg.norm(a) + 1e-9))
        feats[c] = dict(ends=ends, tans=tans, medq=np.median(Q[idx]), size=len(idx))
    return feats

def stitch(P, Q, lab, gap_max=GAP_MAX, cos_link=COS_LINK, qlevel=QLEVEL):
    f = frag_features(P, Q, lab)
    cs = list(f.keys())
    parent = {c: c for c in cs}
    def find(a):
        while parent[a] != a: parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b): parent[find(a)] = find(b)
    for ia in range(len(cs)):
        for ib in range(ia + 1, len(cs)):
            A, B = f[cs[ia]], f[cs[ib]]
            qr = min(A['medq'], B['medq']) / (max(A['medq'], B['medq']) + 1e-9)
            if qr < qlevel:
                continue                          # different dE/dx band -> don't stitch
            linked = False
            for ea, da in zip(A['ends'], A['tans']):
                for eb, db in zip(B['ends'], B['tans']):
                    g = eb - ea; dist = np.linalg.norm(g)
                    if dist > gap_max or dist < 1e-6:
                        continue
                    gh = g / dist
                    # A heads toward B, B heads toward A, and the two outward dirs oppose (smooth)
                    if (gh @ da) > cos_link and (-gh @ db) > cos_link and (da @ db) < -cos_link + 0.2:
                        linked = True
                if linked: break
            if linked:
                union(cs[ia], cs[ib])
    # relabel by root
    roots = {c: find(c) for c in cs}
    remap = {r: i for i, r in enumerate(sorted(set(roots.values())))}
    out = -np.ones(len(P), int)
    for c in cs:
        out[lab == c] = remap[roots[c]]
    return out

def combined(P, Q, qratio_split=0.4):
    lab = cluster(P, Q, qratio=qratio_split, min_hits=4)
    return stitch(P, Q, lab)

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    D = pd.read_csv("labeling/data/real_events.csv")
    if len(sys.argv) > 1 and any(c.isdigit() for c in sys.argv[1]):
        evs = [int(e) for e in sys.argv[1].split(",")]; out = "diagnostics/stitch_events.png"
    else:
        sz = D.groupby('event').size()
        cand = np.array(sz[(sz > 100) & (sz < 1200)].index.to_numpy(), copy=True)
        np.random.seed(3); np.random.shuffle(cand); evs = list(cand[:16]); out = "diagnostics/stitch_gallery.png"
    ncol = 2 if len(evs) <= 6 else 4
    nrow = int(np.ceil(len(evs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4*ncol, 3.1*nrow)); axes = np.atleast_1d(axes).flatten()
    for i, ev in enumerate(evs):
        g = D[D.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        lab = combined(np.stack([x, y, z], 1), q)
        ncl = len(set(lab[lab >= 0])); ng = int((lab < 0).sum())
        axes[i].scatter(x[lab < 0], y[lab < 0], s=3, c='lightgray')
        cols = cm.tab20.colors
        for j2, cl in enumerate(sorted(set(lab[lab >= 0]))):
            m = lab == cl; axes[i].scatter(x[m], y[m], s=7, color=cols[j2 % 20])
        axes[i].set_title(f"evt {ev}: {ncl} tracks, {ng} gray", fontsize=8); axes[i].tick_params(labelsize=6)
    for j in range(len(evs), len(axes)): axes[j].axis('off')
    fig.suptitle("Splitter (dircluster) -> Stitcher, real data", fontsize=11)
    plt.tight_layout(); plt.savefig(out, dpi=95); print("wrote", out)
