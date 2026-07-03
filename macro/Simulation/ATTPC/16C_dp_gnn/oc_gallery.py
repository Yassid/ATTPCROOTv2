#!/usr/bin/env python
"""Apply the sim-trained Object Condensation net to REAL events; render DBSCAN-on-coords clusters."""
import sys, numpy as np, pandas as pd, torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import cm
DEV='cuda' if torch.cuda.is_available() else 'cpu'; IN_DIM=4
EPS=float(sys.argv[1]) if len(sys.argv)>1 else 0.07
def mlp(ch):
    L=[]
    for i in range(len(ch)-1): L+=[nn.Linear(ch[i],ch[i+1]),nn.BatchNorm1d(ch[i+1]),nn.ReLU()]
    return nn.Sequential(*L)
class OCNet(nn.Module):
    def __init__(s,lat=3,k=20):
        super().__init__()
        s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max');s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max')
        s.trunk=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU())
        s.coord=nn.Linear(128,lat); s.beta=nn.Linear(128,1)
    def forward(s,x,b):
        x1=s.ec1(x,b);x2=s.ec2(x1,b);x3=s.ec3(x2,b);h=s.trunk(torch.cat([x1,x2,x3],1))
        return s.coord(h), torch.sigmoid(s.beta(h)).squeeze(-1)
net=OCNet(3).to(DEV); net.load_state_dict(torch.load("data/oc_net.pt",map_location=DEV)); net.eval()
R=pd.read_csv("labeling/data/real_events.csv"); sz=R.groupby('event').size()
cand=np.array(sz[(sz>100)&(sz<1000)].index.to_numpy(),copy=True); np.random.seed(3); np.random.shuffle(cand); evs=cand[:16]
fig,ax=plt.subplots(4,4,figsize=(15,15)); ax=ax.flatten()
for i,ev in enumerate(evs):
    g=R[R.event==ev]; x,y,z,q=(g[c].to_numpy() for c in ['x','y','z','q'])
    f=np.stack([x/250.,y/250.,(z-600.)/450.,np.log1p(np.clip(q,0,None))/9.],1)
    with torch.no_grad(): co,be=net(torch.tensor(f,dtype=torch.float,device=DEV),torch.zeros(len(x),dtype=torch.long,device=DEV))
    pred=DBSCAN(eps=EPS,min_samples=4).fit_predict(co.cpu().numpy())
    ax[i].scatter(x[pred<0],y[pred<0],s=4,c='lightgray')
    for j,cl in enumerate(sorted(set(pred[pred>=0]))): m=pred==cl; ax[i].scatter(x[m],y[m],s=7,color=cm.tab10.colors[j%10])
    ax[i].set_title(f"evt{ev}: {len(set(pred[pred>=0]))} tracks, {len(x)} hits",fontsize=8); ax[i].tick_params(labelsize=6)
fig.suptitle(f"Object Condensation (sim-trained) on REAL data, DBSCAN eps={EPS}",fontsize=13)
plt.tight_layout(); plt.savefig("diagnostics/oc_real.png",dpi=90); print("wrote diagnostics/oc_real.png")
