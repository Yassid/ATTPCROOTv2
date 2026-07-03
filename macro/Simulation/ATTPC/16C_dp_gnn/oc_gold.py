#!/usr/bin/env python
"""Validate Object Condensation (oc_net.pt) on the REAL hand-labeled gold set, vs dircluster.
Same events, same metrics (event-merge% of 2 hand-drawn tracks, ARI, homogeneity)."""
import sys, numpy as np, pandas as pd, torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure
from dircluster import cluster as dircluster
DEV='cuda' if torch.cuda.is_available() else 'cpu'; IN_DIM=4; LAT=int(sys.argv[1]) if len(sys.argv)>1 else 3
NET=sys.argv[2] if len(sys.argv)>2 else "data/oc_net.pt"; EPS=float(sys.argv[3]) if len(sys.argv)>3 else 0.05
def mlp(ch):
    L=[]
    for i in range(len(ch)-1): L+=[nn.Linear(ch[i],ch[i+1]),nn.BatchNorm1d(ch[i+1]),nn.ReLU()]
    return nn.Sequential(*L)
class OCNet(nn.Module):
    def __init__(s,lat,k=20):
        super().__init__()
        s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max');s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max')
        s.trunk=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU())
        s.coord=nn.Linear(128,lat); s.beta=nn.Linear(128,1)
    def forward(s,x,b):
        x1=s.ec1(x,b);x2=s.ec2(x1,b);x3=s.ec3(x2,b);h=s.trunk(torch.cat([x1,x2,x3],1))
        return s.coord(h), torch.sigmoid(s.beta(h)).squeeze(-1)
net=OCNet(LAT).to(DEV); net.load_state_dict(torch.load(NET,map_location=DEV)); net.eval()
G=pd.read_parquet("labeling/data/labels.parquet"); G=G[G.reviewed]; evs=sorted(G.event.unique()); MH=8
def metrics(pred, gold):
    tracks=[t for t in set(gold) if (gold==t).sum()>=MH]
    merged=any(len([t for t in tracks if (gold[pred==c]==t).sum()>=MH])>=2 for c in set(pred[pred>=0]))
    return int(merged), adjusted_rand_score(gold,pred), homogeneity_completeness_v_measure(gold,pred)[0]
def oc_pred(x,y,z,q):
    f=np.stack([x/250.,y/250.,(z-600.)/450.,np.log1p(np.clip(q,0,None))/9.],1)
    with torch.no_grad(): co,be=net(torch.tensor(f,dtype=torch.float,device=DEV),torch.zeros(len(x),dtype=torch.long,device=DEV))
    return DBSCAN(eps=EPS,min_samples=4).fit_predict(co.cpu().numpy())
def dc_pred(x,y,z,q): return dircluster(np.stack([x,y,z],1),q,qratio=0.65,min_hits=4)
def run(fn):
    M=[];A=[];H=[]
    for ev in evs:
        g=G[G.event==ev]; x,y,z,q=(g[c].to_numpy() for c in ['x','y','z','q']); gold=g['label'].to_numpy()
        if len(x)<14: continue
        pred=fn(x,y,z,q); m,a,h=metrics(pred,gold); M.append(m);A.append(a);H.append(h)
    return 100*np.mean(M), np.mean(A), np.mean(H), len(M)
oc=run(oc_pred); dc=run(dc_pred)
print(f"gold reviewed events: {oc[3]}   (OC net={NET} eps={EPS})")
print(f"{'method':<14}{'merge%':>8}{'ARI':>8}{'homog':>8}")
print(f"{'ObjCondense':<14}{oc[0]:>7.1f}%{oc[1]:>8.3f}{oc[2]:>8.3f}")
print(f"{'dircluster':<14}{dc[0]:>7.1f}%{dc[1]:>8.3f}{dc[2]:>8.3f}")
