#!/usr/bin/env python
"""DIFFERENT APPROACH: edge-classification GNN for track separation.

Instead of embedding hits + DBSCAN (which MERGES tracks that touch at the vertex),
predict per k-NN edge whether the two hits are the SAME real track, then cut low-score
edges and take connected components. Vertex-sharing tracks split (their bridging edges
are cut); a conservative (high) cut threshold OVER-segments -> stitchable, never merges.
Trained on noisy-sim truth. Edge label = same POSITIVE track id (noise -1 cut from all).
  ~/gnn_env/bin/python train_edge.py [parquet] [EPOCHS] [K]
"""
import sys, numpy as np, pandas as pd, json, time
import torch, torch.nn as nn
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import knn_graph, GravNetConv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
torch.manual_seed(0); np.random.seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
PARQ  = sys.argv[1] if len(sys.argv) > 1 else "data/sim_noisy.parquet"
EPOCHS = int(sys.argv[2]) if len(sys.argv) > 2 else 25
K = int(sys.argv[3]) if len(sys.argv) > 3 else 10
BATCH, MAX_HITS = 16, 1500

def feats(x, y, z, q):
    r = np.hypot(x, y); phi = np.arctan2(y, x)
    return np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0, r/250.0, phi/np.pi], 1)

def to_data(df):
    out = []
    for gid, g in df.groupby('gid'):
        x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        tid = g['label'].to_numpy().copy()
        tid[g['particle'].to_numpy() == 'noise'] = -1     # noise -> -1 (cut from everything)
        if len(x) < 12: continue
        if len(x) > MAX_HITS:
            idx = np.random.choice(len(x), MAX_HITS, replace=False)
            x, y, z, q, tid = x[idx], y[idx], z[idx], q[idx], tid[idx]
        pos = torch.tensor(np.stack([x, y, z], 1), dtype=torch.float)
        ei = knn_graph(pos, k=K, loop=False)
        t = torch.tensor(tid, dtype=torch.long)
        src, dst = ei
        # same-track edge ONLY if both endpoints share the same POSITIVE track id
        same = ((t[src] == t[dst]) & (t[src] >= 0)).float()
        out.append(Data(x=torch.tensor(feats(x,y,z,q), dtype=torch.float), edge_index=ei, y=same,
                        pos=pos, tid=t, num_nodes=len(x)))
    return out

class EdgeNet(nn.Module):
    def __init__(s, cin=6, h=64, nl=4, k=20):
        super().__init__()
        s.enc = nn.Sequential(nn.Linear(cin,h), nn.ReLU(), nn.Linear(h,h), nn.ReLU())
        s.gn = nn.ModuleList([GravNetConv(h,h,space_dimensions=4,propagate_dimensions=h,k=k) for _ in range(nl)])
        s.ln = nn.ModuleList([nn.LayerNorm(h) for _ in range(nl)])
        s.emlp = nn.Sequential(nn.Linear(2*h,h), nn.ReLU(), nn.Linear(h,h//2), nn.ReLU(), nn.Linear(h//2,1))
    def forward(s, x, ei, batch):
        h = s.enc(x)
        for g, n in zip(s.gn, s.ln): h = n(h + g(h, batch))
        src, dst = ei
        return s.emlp(torch.cat([h[src], h[dst]], -1)).squeeze(-1)

print(f"[{DEV}] loading {PARQ}")
D = pd.read_parquet(PARQ)
gids = D.gid.unique(); np.random.shuffle(gids)
n = len(gids); ntr = int(0.8*n); nva = int(0.1*n)
tr = to_data(D[D.gid.isin(gids[:ntr])]); va = to_data(D[D.gid.isin(gids[ntr:ntr+nva])])
te = to_data(D[D.gid.isin(gids[ntr+nva:])])
posf = np.mean([d.y.mean().item() for d in tr])
print(f"events: {len(tr)}/{len(va)}/{len(te)}   same-edge fraction {posf:.3f}")

net = EdgeNet().to(DEV)
opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-5)
# emphasise the rarer 'different' edges (the cuts) via pos_weight on positives
pw = torch.tensor([(1-posf)/max(posf,1e-3)], device=DEV)
lossf = nn.BCEWithLogitsLoss(pos_weight=pw)
tl = DataLoader(tr, batch_size=BATCH, shuffle=True); vl = DataLoader(va, batch_size=BATCH)
t0 = time.time(); best = 1e9
for ep in range(EPOCHS):
    net.train(); s = 0
    for d in tl:
        d = d.to(DEV); opt.zero_grad(); out = net(d.x, d.edge_index, d.batch)
        loss = lossf(out, d.y); loss.backward(); opt.step(); s += loss.item()
    net.eval(); vs = 0; tp = fp = tn = fn = 0
    with torch.no_grad():
        for d in vl:
            d = d.to(DEV); out = net(d.x, d.edge_index, d.batch); vs += lossf(out, d.y).item()
            p = (out > 0).float()
            tp += ((p==1)&(d.y==1)).sum().item(); fp += ((p==1)&(d.y==0)).sum().item()
            tn += ((p==0)&(d.y==0)).sum().item(); fn += ((p==0)&(d.y==1)).sum().item()
    vs /= max(len(vl),1)
    if vs < best: best = vs; torch.save(net.state_dict(), "data/edge_sim_noisy.pt")
    if ep % 3 == 0 or ep == EPOCHS-1:
        prec = tp/max(tp+fp,1); rec = tp/max(tp+fn,1)
        print(f"ep {ep:2d} tr {s/len(tl):.3f} va {vs:.3f} best {best:.3f} | same-edge P {prec:.3f} R {rec:.3f} "
              f"cut-edge acc {tn/max(tn+fp,1):.3f} ({time.time()-t0:.0f}s)")
net.load_state_dict(torch.load("data/edge_sim_noisy.pt", map_location=DEV)); net.eval()

# --- eval: connected components at cut threshold -> merge/fragment metrics on held-out sim ---
def cc(edge_index, score, n, thr):
    keep = score > thr
    s, d = edge_index[0][keep], edge_index[1][keep]
    A = coo_matrix((np.ones(len(s)), (s, d)), shape=(n, n))
    A = A + A.T
    ncomp, lab = connected_components(A, directed=False)
    return lab

def merges(lab, tid, min_hits=8):
    m = 0
    for c in set(lab):
        tt = tid[(lab == c) & (tid >= 0)]
        big = [t for t in set(tt) if np.sum(tt == t) >= min_hits]
        if len(big) >= 2: m += 1
    return m

print(f"\n{'thr':>5} {'comps/ev':>9} {'ev w/merge%':>11} {'sameEdgeR':>10}")
for thr in [0.0, 0.5, 1.0, 2.0, 3.0]:
    comps = []; evm = 0
    for d in te:
        with torch.no_grad():
            out = net(d.x.to(DEV), d.edge_index.to(DEV),
                      torch.zeros(d.num_nodes, dtype=torch.long, device=DEV)).cpu().numpy()
        lab = cc(d.edge_index.numpy(), out, d.num_nodes, thr)
        tid = d.tid.numpy()
        # count only components with >=8 real-track hits as tracks
        ncomp_real = len(set(lab[(tid >= 0)]))
        comps.append(ncomp_real); evm += (merges(lab, tid) > 0)
    print(f"{thr:>5} {np.mean(comps):>9.2f} {100*evm/len(te):>11.1f}")
json.dump({'epochs': EPOCHS, 'K': K, 'n_train': len(tr)}, open("data/metrics_edge.json", "w"))
print("saved data/edge_sim_noisy.pt")
