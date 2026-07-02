#!/usr/bin/env python
"""Noise cleaning as a standalone per-hit filter (no clustering commitment).

A hit is SIGNAL if it has >= MIN_DEG neighbours that are direction- AND charge-continuous
with it (i.e. it sits on a locally linear, smooth-dE/dx structure). Noise is isolated /
off-line and gets degree ~0. Optimised to KEEP signal (incl. sparse track tails) while
removing scattered background. Feeds any downstream track-finder.

  signal_mask = denoise(P, Q)
  ~/gnn_env/bin/python denoise.py            # ROC on sim truth + real-data gallery
"""
import sys, numpy as np, pandas as pd
from sklearn.neighbors import NearestNeighbors

K, COS_SEG, COS_TAN, RMAX, QRATIO, MIN_DEG = 12, 0.78, 0.72, 32.0, 0.65, 1

def denoise(P, Q, k=K, cos_seg=COS_SEG, cos_tan=COS_TAN, rmax=RMAX, qratio=QRATIO,
            min_deg=MIN_DEG, smooth_q=True, return_deg=False):
    n = len(P)
    if n < k + 1:
        return (np.zeros(n, bool) if not return_deg else np.zeros(n))
    nn = NearestNeighbors(n_neighbors=k+1).fit(P); _, idx = nn.kneighbors(P)
    nb = idx[:, 1:]
    d = P[nb] - P[nb].mean(1, keepdims=True)
    cov = np.einsum('nki,nkj->nij', d, d) / k
    _, v = np.linalg.eigh(cov); tang = v[:, :, 2]
    Qeff = np.median(Q[nb], axis=1) if smooth_q else Q
    src = np.repeat(np.arange(n), k); dst = nb.ravel()
    seg = P[dst] - P[src]; dl = np.linalg.norm(seg, axis=1) + 1e-9; su = seg / dl[:, None]
    ai = np.abs((su * tang[src]).sum(1)); aj = np.abs((su * tang[dst]).sum(1))
    at = np.abs((tang[src] * tang[dst]).sum(1))
    qr = np.minimum(Qeff[src], Qeff[dst]) / (np.maximum(Qeff[src], Qeff[dst]) + 1e-9)
    keep = (ai > cos_seg) & (aj > cos_seg) & (at > cos_tan) & (dl < rmax) & (qr > qratio)
    deg = np.zeros(n)
    np.add.at(deg, src[keep], 1.0); np.add.at(deg, dst[keep], 1.0)   # undirected support
    if return_deg:
        return deg
    return deg >= min_deg

if __name__ == "__main__":
    D = pd.read_parquet("data/sim_noisy.parquet")
    np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids); gids = gids[:250]
    EV = []
    for gid, g in D[D.gid.isin(gids)].groupby('gid'):
        x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q']); part = g['particle'].to_numpy()
        if len(x) < 14: continue
        EV.append((np.stack([x,y,z],1), q, part))
    # ROC over min_deg: signal eff (per particle) vs noise rejection
    print(f"{'min_deg':>7}{'proton_eff':>11}{'17C_eff':>9}{'allsig_eff':>11}{'noise_rej':>10}")
    for md in [1, 2, 3, 4]:
        kept = {'proton':[0,0], '17C_recoil':[0,0], 'sig':[0,0], 'noise':[0,0]}
        for P, Q, part in EV:
            deg = denoise(P, Q, return_deg=True); sig = deg >= md
            for p in ['proton', '17C_recoil']:
                m = part == p; kept[p][0] += int(sig[m].sum()); kept[p][1] += int(m.sum())
            st = part != 'noise'; kept['sig'][0] += int(sig[st].sum()); kept['sig'][1] += int(st.sum())
            nm = part == 'noise'; kept['noise'][0] += int((~sig[nm]).sum()); kept['noise'][1] += int(nm.sum())
        f = lambda a: 100*a[0]/max(a[1],1)
        print(f"{md:>7}{f(kept['proton']):>10.1f}%{f(kept['17C_recoil']):>8.1f}%{f(kept['sig']):>10.1f}%{f(kept['noise']):>9.1f}%")
    print("\nwant: high signal eff (esp 17C), high noise rejection. min_deg trades them.")

    # real-data gallery: kept (color) vs removed (gray)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    R = pd.read_csv("labeling/data/real_events.csv")
    sz = R.groupby('event').size(); cand = np.array(sz[(sz>100)&(sz<1200)].index.to_numpy(), copy=True)
    np.random.seed(3); np.random.shuffle(cand); evs = cand[:9]
    fig, axes = plt.subplots(3, 3, figsize=(11, 10)); axes = axes.flatten()
    for i, ev in enumerate(evs):
        g = R[R.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        sig = denoise(np.stack([x,y,z],1), q, min_deg=2)
        axes[i].scatter(x[~sig], y[~sig], s=4, c='lightgray')
        axes[i].scatter(x[sig], y[sig], s=7, c='tab:blue')
        axes[i].set_title(f"evt{ev}: {sig.sum()}/{len(x)} kept", fontsize=8); axes[i].tick_params(labelsize=6)
    fig.suptitle("Noise cleaning (min_deg=2): blue=signal, gray=removed", fontsize=12)
    plt.tight_layout(); plt.savefig("diagnostics/denoise_gallery.png", dpi=95); print("wrote diagnostics/denoise_gallery.png")
