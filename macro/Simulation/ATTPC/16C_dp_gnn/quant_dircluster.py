#!/usr/bin/env python
"""Quantify direction-continuity clustering on SIM truth (merge-focused, stitch-aware).
Sweeps the dE/dx charge-continuity knob QRATIO (edge features computed once per event,
only the charge threshold changes) to show whether charge helps cut proton/17C at the vertex.
  ~/gnn_env/bin/python quant_dircluster.py [nEvents]
"""
import sys, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import homogeneity_completeness_v_measure
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
K, COS_SEG, COS_TAN, RMAX, MIN_HITS = 12, 0.78, 0.72, 32.0, 3
NEV = int(sys.argv[1]) if len(sys.argv) > 1 else 300

SMOOTH_DEDX = True   # use median neighbourhood charge (local dE/dx) instead of raw single-hit charge

def edge_feats(P, Q):
    n = len(P)
    nn = NearestNeighbors(n_neighbors=K+1).fit(P); _, idx = nn.kneighbors(P)
    nb = idx[:, 1:]
    d = P[nb] - P[nb].mean(1, keepdims=True)
    cov = np.einsum('nki,nkj->nij', d, d) / K
    _, v = np.linalg.eigh(cov); tang = v[:, :, 2]
    # local dE/dx = median charge over the hit's own neighbourhood (robust to single-pad noise)
    Qloc = np.median(Q[nb], axis=1) if SMOOTH_DEDX else Q
    src = np.repeat(np.arange(n), K); dst = nb.ravel()
    seg = P[dst] - P[src]; dl = np.linalg.norm(seg, axis=1) + 1e-9; su = seg / dl[:, None]
    ai = np.abs((su * tang[src]).sum(1)); aj = np.abs((su * tang[dst]).sum(1))
    at = np.abs((tang[src] * tang[dst]).sum(1))
    qr = np.minimum(Qloc[src], Qloc[dst]) / (np.maximum(Qloc[src], Qloc[dst]) + 1e-9)
    geo = (ai > COS_SEG) & (aj > COS_SEG) & (at > COS_TAN) & (dl < RMAX)
    return src, dst, geo, qr, n

def labels(src, dst, geo, qr, n, qratio):
    keep = geo & (qr > qratio)
    A = coo_matrix((np.ones(keep.sum()), (src[keep], dst[keep])), shape=(n, n)); A = A + A.T
    _, lab = connected_components(A, directed=False)
    out = -np.ones(n, int); nxt = 0
    for c in range(lab.max() + 1):
        m = lab == c
        if m.sum() >= MIN_HITS: out[m] = nxt; nxt += 1
    return out

def metrics(pred, tid, mh=8):
    tracks = [t for t in set(tid[tid >= 0]) if (tid == t).sum() >= mh]
    nmerge = 0
    for c in set(pred[pred >= 0]):
        tt = tid[(pred == c)]; big = [t for t in tracks if (tt == t).sum() >= mh]
        if len(big) >= 2: nmerge += 1
    frags = []; rec = 0
    for t in tracks:
        tm = tid == t
        cs = [c for c in set(pred[tm & (pred >= 0)]) if ((pred == c) & tm).sum() >= mh]
        frags.append(len(cs))
        clustered = (pred[tm] >= 0).mean()
        merged = any(len([t2 for t2 in tracks if (tid[pred == c] == t2).sum() >= mh]) >= 2 for c in cs)
        if clustered >= 0.6 and not merged: rec += 1
    nm = tid < 0
    noise_in = (pred[nm] >= 0).mean() if nm.sum() else 0.0
    h, _, _ = homogeneity_completeness_v_measure(tid, pred)
    return nmerge, len(tracks), frags, rec, noise_in, h

D = pd.read_parquet("data/sim_noisy.parquet")
np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids); gids = gids[:NEV]
sub = D[D.gid.isin(gids)]
E = []
for gid, g in sub.groupby('gid'):
    x, y, z, q = (g[c].to_numpy() for c in ['x', 'y', 'z', 'q'])
    tid = g['label'].to_numpy().copy(); tid[g['particle'].to_numpy() == 'noise'] = -1
    if len(x) < K + 1: continue
    E.append((edge_feats(np.stack([x, y, z], 1), q), tid))
print(f"events: {len(E)}  (multi-track: {sum(len(set(t[t>=0]))>=2 for _,t in E)})\n")
print(f"{'QRATIO':>6} {'ev_merge%':>9} {'homog':>6} {'frags/trk':>9} {'recov%':>7} {'noise_in%':>9} {'clus/ev':>8}")
for qratio in [0.0, 0.2, 0.3, 0.4, 0.5, 0.6]:
    evm = 0; H = []; F = []; rec = 0; ntr = 0; NI = []; NC = []
    for (src, dst, geo, qr, n), tid in E:
        pred = labels(src, dst, geo, qr, n, qratio)
        nmerge, nt, frags, r, ni, h = metrics(pred, tid)
        evm += (nmerge > 0); H.append(h); F += frags; rec += r; ntr += nt; NI.append(ni)
        NC.append(len(set(pred[pred >= 0])))
    print(f"{qratio:>6} {100*evm/len(E):>9.1f} {np.mean(H):>6.3f} {np.mean(F):>9.2f} "
          f"{100*rec/max(ntr,1):>7.1f} {100*np.mean(NI):>9.1f} {np.mean(NC):>8.2f}")
print("\nwant: low ev_merge%, high homog, high recov%, low noise_in%. frags/trk = stitcher workload.")
