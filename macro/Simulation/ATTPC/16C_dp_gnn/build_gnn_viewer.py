#!/usr/bin/env python
# Run the HDBSCAN-trained GNN on run_0305 and build matched CSVs for the GNN-vs-Spyral viewer.
#   GNN clustering (left) vs Spyral (right), same global GET event #. ATTPCROOT z flipped to Spyral frame.
# Out: data/cmp_gnn_0305.csv (+ reuses data/cmp_spyral_0305.csv). Also a static preview.
import json
import pandas as pd, numpy as np
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RUN, OFF, ZFLIP = 305, 34675, 1137.0
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM = EMB_DIM = 4; K = 20
EPS = json.load(open("sup/metrics_attpc.json")).get("best_eps", 0.7)

def mlp(ch):
    L = []
    for i in range(len(ch)-1):
        L += [nn.Linear(ch[i], ch[i+1]), nn.BatchNorm1d(ch[i+1]), nn.ReLU()]
    return nn.Sequential(*L)
class Net(nn.Module):
    def __init__(s, k=20):
        super().__init__()
        s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max'); s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max'); s.head=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU(),nn.Linear(128,EMB_DIM))
    def forward(s,x,b):
        x1=s.ec1(x,b); x2=s.ec2(x1,b); x3=s.ec3(x2,b); return s.head(torch.cat([x1,x2,x3],1))

a = pd.read_parquet("data/attpc_all.parquet"); a = a[a.run == RUN].copy()
m = Net(K).to(DEV); m.load_state_dict(torch.load("sup/embed_attpc.pt", map_location=DEV)); m.eval()
print(f"GNN inference (eps={EPS}) on {a.gid.nunique()} run_0{RUN} events ...")
rows = []
for gid, g in a.groupby('gid'):
    x=g.x.values; y=g.y.values; z=g.z.values; q=g.q.values; n=len(x); lab=np.full(n,-1)
    if n>=10:
        feat=np.stack([x/250,y/250,(z-600)/450,np.log1p(np.clip(q,0,None))/9.0],1)
        with torch.no_grad():
            emb=m(torch.tensor(feat,dtype=torch.float,device=DEV),torch.zeros(n,dtype=torch.long,device=DEV)).cpu().numpy()
        lab=DBSCAN(eps=EPS,min_samples=6).fit_predict(emb)
    gg=g.copy(); gg['gnn']=lab; gg['glob']=OFF+gg['event']; rows.append(gg)
gn=pd.concat(rows)

sp=pd.read_parquet("data/spyral_all.parquet"); sp=sp[sp.run==RUN]
common=sorted(set(gn['glob'])&set(sp.event))
go=gn[gn['glob'].isin(common)].copy(); go['z']=ZFLIP-go['z']
go[['glob','x','y','z','q','gnn']].rename(columns={'glob':'event','gnn':'cluster'}).to_csv("data/cmp_gnn_0305.csv",index=False)
spo=sp[sp.event.isin(common)][['event','x','y','z','q','label']].rename(columns={'label':'cluster'})
spo.to_csv("data/cmp_spyral_0305.csv",index=False)
print(f"{len(common)} matched events -> data/cmp_gnn_0305.csv / data/cmp_spyral_0305.csv")

# static preview on the crossing events (entry idx -> global)
PAL=plt.cm.tab10.colors; col=lambda c:(0.7,0.7,0.7,0.3) if c<0 else PAL[int(c)%10]
CROSS=[OFF+e for e in [60,225,309,331,354]]
g2=pd.read_csv("data/cmp_gnn_0305.csv"); s2=pd.read_csv("data/cmp_spyral_0305.csv")
ev=[e for e in CROSS if e in set(g2.event)][:5]
fig,ax=plt.subplots(len(ev),2,figsize=(13,3.3*len(ev)))
for r,e in enumerate(ev):
    gg=g2[g2.event==e]; ss=s2[s2.event==e]
    for cl in sorted(gg.cluster.unique()): mk=gg.cluster==cl; ax[r,0].scatter(gg.z[mk],gg.y[mk],s=8,color=col(cl))
    for cl in sorted(ss.cluster.unique()): mk=ss.cluster==cl; ax[r,1].scatter(ss.z[mk],ss.y[mk],s=8,color=col(cl))
    ax[r,0].set_title(f"GNN ev{e} ({gg[gg.cluster>=0].cluster.nunique()} cl)",fontsize=9)
    ax[r,1].set_title(f"Spyral ev{e} ({ss.cluster.nunique()} cl)",fontsize=9)
    for c in range(2): ax[r,c].set_xlabel('z'); ax[r,c].set_ylabel('y')
plt.tight_layout(); plt.savefig("gnn_separation_preview.png",dpi=80); print("-> gnn_separation_preview.png")
