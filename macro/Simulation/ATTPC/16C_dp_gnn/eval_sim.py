#!/usr/bin/env python
"""Evaluate embed_sim_noisy.pt on held-out SIM (truth labels), with metrics that reward
over-segmentation and punish MERGES (under-segmentation), since fragments are stitchable.
  ~/gnn_env/bin/python eval_sim.py
"""
import numpy as np, pandas as pd
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure
np.random.seed(0); torch.manual_seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, EMB_DIM = 4, 4

def mlp(ch):
    L = []
    for i in range(len(ch)-1): L += [nn.Linear(ch[i],ch[i+1]), nn.BatchNorm1d(ch[i+1]), nn.ReLU()]
    return nn.Sequential(*L)
class Net(nn.Module):
    def __init__(s, k=20):
        super().__init__()
        s.ec1 = DynamicEdgeConv(mlp([2*IN_DIM,64,64]), k, 'max'); s.ec2 = DynamicEdgeConv(mlp([2*64,64,64]), k, 'max')
        s.ec3 = DynamicEdgeConv(mlp([2*64,128]), k, 'max')
        s.head = nn.Sequential(nn.Linear(256,128), nn.BatchNorm1d(128), nn.ReLU(), nn.Linear(128,EMB_DIM))
    def forward(s, x, b):
        x1 = s.ec1(x,b); x2 = s.ec2(x1,b); x3 = s.ec3(x2,b); return s.head(torch.cat([x1,x2,x3],1))

net = Net().to(DEV); net.load_state_dict(torch.load("data/embed_sim_noisy.pt", map_location=DEV)); net.eval()

D = pd.read_parquet("data/sim_noisy.parquet")
gids = D.gid.unique(); np.random.shuffle(gids)   # same seed/order as train_sim.py
n = len(gids); ntr = int(0.8*n); nva = int(0.1*n); te_g = gids[ntr+nva:]
test = D[D.gid.isin(te_g)]

# embed each test event
events = []
for gid, g in test.groupby('gid'):
    x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    if len(x) < 10: continue
    truth = g['label'].to_numpy()
    is_track = (g['particle'].to_numpy() != 'noise')   # real tracks (proton/17C), exclude noise-haze label
    feat = np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0], 1)
    with torch.no_grad():
        emb = net(torch.tensor(feat, dtype=torch.float, device=DEV),
                  torch.zeros(len(x), dtype=torch.long, device=DEV)).cpu().numpy()
    events.append((emb, truth, is_track))
print(f"held-out sim events: {len(events)}")

def merge_metrics(pred, truth, is_track, min_hits=8):
    """A predicted cluster is a MERGE if it contains >=2 distinct real-track truth ids each with
    >=min_hits hits. Returns (#merged_clusters, #real_tracks, #tracks_recovered_pure)."""
    merged = 0
    for c in set(pred[pred >= 0]):
        m = pred == c
        tl = truth[m & is_track]
        big = [t for t in set(tl) if np.sum(tl == t) >= min_hits]
        if len(big) >= 2: merged += 1
    # recovery: a truth track is 'recovered' if some predicted cluster holds >=60% of its hits AND
    # that cluster is >=80% pure for this track (over its track hits)
    rec = 0; ntr_ = 0
    for t in set(truth[is_track]):
        tm = (truth == t) & is_track
        if tm.sum() < min_hits: continue
        ntr_ += 1
        best = 0
        for c in set(pred[pred >= 0]):
            cm = pred == c
            inter = (cm & tm).sum()
            frac_track = inter / tm.sum()
            frac_pure = inter / max((cm & is_track).sum(), 1)
            if frac_track >= 0.6 and frac_pure >= 0.8: best = 1
        rec += best
    return merged, ntr_, rec

print(f"\n{'eps':>5} {'ARI':>6} {'homog':>6} {'compl':>6} {'clus/ev':>8} {'merged/ev':>10} {'ev w/merge%':>11} {'track_rec%':>10}")
for eps in [0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7]:
    aris = []; homs = []; comps = []; nclus = []; nmerge = []; evmerge = 0
    rec_tot = 0; trk_tot = 0
    for emb, truth, is_track in events:
        pred = DBSCAN(eps=eps, min_samples=6).fit_predict(emb)
        aris.append(adjusted_rand_score(truth, pred))
        h, c, v = homogeneity_completeness_v_measure(truth, pred)
        homs.append(h); comps.append(c)
        nclus.append(len(set(pred[pred >= 0])))
        mg, ntr_, rec = merge_metrics(pred, truth, is_track)
        nmerge.append(mg); evmerge += (mg > 0); rec_tot += rec; trk_tot += ntr_
    print(f"{eps:>5} {np.mean(aris):>6.3f} {np.mean(homs):>6.3f} {np.mean(comps):>6.3f} "
          f"{np.mean(nclus):>8.2f} {np.mean(nmerge):>10.2f} {100*evmerge/len(events):>11.1f} "
          f"{100*rec_tot/max(trk_tot,1):>10.1f}")
print("\nover-seg is fine (stitchable) -> want high homog, low merged/ev & ev-w-merge, high track_rec")
