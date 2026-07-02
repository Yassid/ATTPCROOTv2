#!/usr/bin/env python
"""Galleries of the sim-trained DGCNN separating tracks.
  real:  ~/gnn_env/bin/python gallery.py real [eps] [nEvents] [out.png]
  sim :  ~/gnn_env/bin/python gallery.py sim  [eps] [nEvents] [out.png]   (truth vs pred pairs)
"""
import sys, numpy as np, pandas as pd
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, EMB_DIM = 4, 4
MODE = sys.argv[1] if len(sys.argv) > 1 else "real"
EPS  = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2
NEV  = int(sys.argv[3]) if len(sys.argv) > 3 else 16
OUT  = sys.argv[4] if len(sys.argv) > 4 else f"diagnostics/gallery_{MODE}.png"
MINSAMP = int(sys.argv[5]) if len(sys.argv) > 5 else 3   # lower keeps sparse track segments (stitchable)
EVLIST = [int(e) for e in sys.argv[6].split(",")] if len(sys.argv) > 6 else None  # explicit event ids

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

def embed(x, y, z, q):
    feat = np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0], 1)
    with torch.no_grad():
        return net(torch.tensor(feat, dtype=torch.float, device=DEV),
                   torch.zeros(len(x), dtype=torch.long, device=DEV)).cpu().numpy()

def draw(ax, x, y, lab, title):
    noise = lab < 0
    ax.scatter(x[noise], y[noise], s=3, c='lightgray')
    cols = cm.tab10.colors
    for i, cl in enumerate(sorted(set(lab[lab >= 0]))):
        m = lab == cl
        ax.scatter(x[m], y[m], s=6, color=cols[i % 10])
    ax.set_title(title, fontsize=8); ax.set_xlabel("x[mm]", fontsize=7); ax.set_ylabel("y[mm]", fontsize=7)
    ax.tick_params(labelsize=6)

if MODE == "real":
    D = pd.read_csv("labeling/data/real_events.csv")
    sizes = D.groupby('event').size()
    if EVLIST is not None:
        cand = np.array(EVLIST)
    else:
        cand = np.array(sizes[(sizes > 100) & (sizes < 1200)].index.to_numpy(), copy=True)
        np.random.seed(3); np.random.shuffle(cand); cand = cand[:NEV]
    ncol = 4; nrow = int(np.ceil(len(cand)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2*ncol, 3.0*nrow)); axes = axes.flatten()
    for i, ev in enumerate(cand):
        g = D[D.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        pred = DBSCAN(eps=EPS, min_samples=MINSAMP).fit_predict(embed(x, y, z, q))
        ncl = len(set(pred[pred >= 0]))
        draw(axes[i], x, y, pred, f"real {ev}: {ncl} clusters, {len(x)} hits")
    for j in range(len(cand), len(axes)): axes[j].axis('off')
    fig.suptitle(f"SIM-trained DGCNN on REAL data (eps={EPS}, gray=noise)", fontsize=11)
else:
    D = pd.read_parquet("data/sim_noisy.parquet")
    np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids)
    n = len(gids); te = gids[int(0.9*n):]        # held-out
    np.random.seed(5); sel = np.array(te, copy=True); np.random.shuffle(sel); sel = sel[:NEV]
    ncol = 4; nrow = int(np.ceil(len(sel)/ncol))
    fig, axes = plt.subplots(nrow, ncol*2, figsize=(3.0*ncol*2, 2.8*nrow))
    for i, ev in enumerate(sel):
        g = D[D.gid == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        truth = g['label'].to_numpy()
        pred = DBSCAN(eps=EPS, min_samples=MINSAMP).fit_predict(embed(x, y, z, q))
        r, c = divmod(i, ncol)
        draw(axes[r, 2*c],   x, y, truth, f"TRUTH {ev}")
        draw(axes[r, 2*c+1], x, y, pred,  f"PRED {len(set(pred[pred>=0]))}cl")
    fig.suptitle(f"Held-out SIM: TRUTH vs PRED (eps={EPS})", fontsize=11)
plt.tight_layout(); plt.savefig(OUT, dpi=95)
print("wrote", OUT)
