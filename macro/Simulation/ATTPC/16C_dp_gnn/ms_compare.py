#!/usr/bin/env python
"""Show DBSCAN min_samples effect on sparse-track recovery for chosen real events.
Columns = min_samples [6,3,2]; rows = events.  ~/gnn_env/bin/python ms_compare.py ev1,ev2,...
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
EVS = [int(e) for e in (sys.argv[1] if len(sys.argv) > 1 else "2272,631,2078,1232").split(",")]
EPS = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2
MSS = [6, 3, 2]

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
D = pd.read_csv("labeling/data/real_events.csv")

def emb(x, y, z, q):
    f = np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0], 1)
    with torch.no_grad():
        return net(torch.tensor(f, dtype=torch.float, device=DEV),
                   torch.zeros(len(x), dtype=torch.long, device=DEV)).cpu().numpy()

fig, axes = plt.subplots(len(EVS), len(MSS), figsize=(3.3*len(MSS), 3.0*len(EVS)))
if len(EVS) == 1: axes = axes[None, :]
cols = cm.tab10.colors
for r, ev in enumerate(EVS):
    g = D[D.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    e = emb(x, y, z, q)
    for c, ms in enumerate(MSS):
        pred = DBSCAN(eps=EPS, min_samples=ms).fit_predict(e)
        ngray = int((pred < 0).sum()); ncl = len(set(pred[pred >= 0]))
        ax = axes[r, c]
        m0 = pred < 0
        ax.scatter(x[m0], y[m0], s=4, c='lightgray')
        for i, cl in enumerate(sorted(set(pred[pred >= 0]))):
            m = pred == cl; ax.scatter(x[m], y[m], s=7, color=cols[i % 10])
        ax.set_title(f"evt {ev}  ms={ms}: {ncl}cl, {ngray} gray", fontsize=8)
        ax.tick_params(labelsize=6)
plt.tight_layout(); plt.savefig("diagnostics/ms_compare.png", dpi=95)
print("wrote diagnostics/ms_compare.png")
