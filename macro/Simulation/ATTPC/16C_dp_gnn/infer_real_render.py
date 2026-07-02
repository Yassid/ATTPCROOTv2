#!/usr/bin/env python
"""Apply the SIM-trained DGCNN to REAL events and render its cluster separation.
The real test of sim->real transfer: does a net trained only on noisy simulation
separate tracks in real data? (Visual > metrics.)
  ~/gnn_env/bin/python infer_real_render.py [real_events.csv] [model.pt] [eps] [nEvents] [out.png]
"""
import sys, numpy as np, pandas as pd
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, EMB_DIM = 4, 4
REAL = sys.argv[1] if len(sys.argv) > 1 else "labeling/data/real_events.csv"
MODEL = sys.argv[2] if len(sys.argv) > 2 else "data/embed_sim_noisy.pt"
EPS = float(sys.argv[3]) if len(sys.argv) > 3 else 0.7
NEV = int(sys.argv[4]) if len(sys.argv) > 4 else 6
OUT = sys.argv[5] if len(sys.argv) > 5 else "diagnostics/sim2real_clusters.png"

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

net = Net().to(DEV); net.load_state_dict(torch.load(MODEL, map_location=DEV)); net.eval()
D = pd.read_csv(REAL)
# pick events with a reasonable number of hits
sizes = D.groupby('event').size()
cand = np.array(sizes[(sizes > 80) & (sizes < 1500)].index.to_numpy(), copy=True)
np.random.seed(1); np.random.shuffle(cand)
cand = cand[:NEV]

fig, axes = plt.subplots(len(cand), 2, figsize=(9, 3.2*len(cand)))
if len(cand) == 1: axes = axes[None, :]
for r, ev in enumerate(cand):
    g = D[D.event == ev]
    x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    feat = np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0], 1)
    with torch.no_grad():
        emb = net(torch.tensor(feat, dtype=torch.float, device=DEV),
                  torch.zeros(len(x), dtype=torch.long, device=DEV)).cpu().numpy()
    pred = DBSCAN(eps=EPS, min_samples=6).fit_predict(emb)
    ncl = len(set(pred[pred >= 0]))
    for c, (ax, aa, bb, xl, yl) in enumerate([(axes[r,0], x, y, 'x[mm]', 'y[mm]'),
                                              (axes[r,1], z, x, 'z[mm]', 'x[mm]')]):
        noise = pred < 0
        ax.scatter(aa[noise], bb[noise], s=3, c='lightgray')
        for cl in sorted(set(pred[pred >= 0])):
            m = pred == cl
            ax.scatter(aa[m], bb[m], s=5)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        if c == 0: ax.set_title(f"real evt {ev}: {ncl} clusters ({len(x)} hits, gray=noise)")
        else: ax.set_title("z-x")
plt.tight_layout(); plt.savefig(OUT, dpi=90)
print("wrote", OUT, "eps", EPS)
