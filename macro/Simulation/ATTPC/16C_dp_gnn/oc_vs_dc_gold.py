#!/usr/bin/env python
"""Illustrate the OC sim->real failure: GOLD hand labels | OC pred | dircluster pred on real
2-track gold events. OC merges the two tracks; dircluster keeps them separate."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
from dircluster import cluster as dircluster
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
DEV='cuda' if torch.cuda.is_available() else 'cpu'; IN_DIM=4
def mlp(ch):
    L=[]
    for i in range(len(ch)-1): L+=[nn.Linear(ch[i],ch[i+1]),nn.BatchNorm1d(ch[i+1]),nn.ReLU()]
    return nn.Sequential(*L)
class OCNet(nn.Module):
    def __init__(s,lat=3,k=20):
        super().__init__(); s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max');s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max'); s.trunk=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU())
        s.coord=nn.Linear(128,lat); s.beta=nn.Linear(128,1)
    def forward(s,x,b):
        x1=s.ec1(x,b);x2=s.ec2(x1,b);x3=s.ec3(x2,b);h=s.trunk(torch.cat([x1,x2,x3],1)); return s.coord(h),torch.sigmoid(s.beta(h)).squeeze(-1)
net=OCNet().to(DEV); net.load_state_dict(torch.load('data/oc_net.pt',map_location=DEV)); net.eval()
G=pd.read_parquet('labeling/data/labels.parquet'); G=G[G.reviewed]
# pick clear 2-track events (both labels >=20 hits)
cand=[]
for ev in sorted(G.event.unique()):
    g=G[G.event==ev]; vc=g.label.value_counts()
    if (vc>=20).sum()>=2 and len(g)<600: cand.append(ev)
np.random.seed(1); np.random.shuffle(cand); evs=cand[:6]
def ocp(x,y,z,q):
    f=np.stack([x/250.,y/250.,(z-600.)/450.,np.log1p(np.clip(q,0,None))/9.],1)
    with torch.no_grad(): co,_=net(torch.tensor(f,dtype=torch.float,device=DEV),torch.zeros(len(x),dtype=torch.long,device=DEV))
    return DBSCAN(eps=0.03,min_samples=4).fit_predict(co.cpu().numpy())
def draw(ax,x,y,lab,title):
    ax.scatter(x[lab<0],y[lab<0],s=4,c='lightgray')
    for i,cl in enumerate(sorted(set(lab[lab>=0]))): m=lab==cl; ax.scatter(x[m],y[m],s=8,color=cm.tab10.colors[i%10])
    ax.set_title(title,fontsize=8); ax.tick_params(labelsize=6)
fig,axes=plt.subplots(len(evs),3,figsize=(10,3*len(evs)))
for r,ev in enumerate(evs):
    g=G[G.event==ev]; x,y,z,q=(g[c].to_numpy() for c in ['x','y','z','q']); gold=g.label.to_numpy()
    op=ocp(x,y,z,q); dp=dircluster(np.stack([x,y,z],1),q,qratio=0.65,min_hits=4)
    draw(axes[r,0],x,y,gold,f"GOLD evt{ev} ({len(set(gold))} tr)")
    draw(axes[r,1],x,y,op,f"OC ({len(set(op[op>=0]))} cl)")
    draw(axes[r,2],x,y,dp,f"dircluster ({len(set(dp[dp>=0]))} cl)")
fig.suptitle("Real gold events: GOLD | Object Condensation (merges) | dircluster (separates)",fontsize=11)
plt.tight_layout(); plt.savefig("diagnostics/oc_vs_dc_gold.png",dpi=95); print("wrote diagnostics/oc_vs_dc_gold.png")
