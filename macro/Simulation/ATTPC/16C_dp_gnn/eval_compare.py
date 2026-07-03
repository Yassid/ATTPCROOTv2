#!/usr/bin/env python
"""Head-to-head on held-out sim: Object Condensation (oc_net.pt) vs dircluster, same events,
same metrics (event-merge%, track recovery%, per-particle kept%). Uses full events."""
import numpy as np, pandas as pd, torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from dircluster import cluster as dircluster
np.random.seed(0); torch.manual_seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, LAT = 4, 3

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

def cluster_oc(coord, beta, t_beta=0.1, t_d=0.5):
    n=len(beta); lab=-np.ones(n,int); un=np.ones(n,bool); cid=0
    for i in np.argsort(-beta):
        if beta[i]<t_beta: break
        if not un[i]: continue
        sel=un&(np.linalg.norm(coord-coord[i],axis=1)<t_d); lab[sel]=cid; un[sel]=False; cid+=1
    return lab

def metrics(pred, tid, part, mh=8):
    tracks=[t for t in set(tid[tid>=0]) if (tid==t).sum()>=mh]
    merged=any(len([t for t in tracks if (tid[pred==c]==t).sum()>=mh])>=2 for c in set(pred[pred>=0]))
    rec=sum(any(((pred==c)&(tid==t)).sum()/(tid==t).sum()>=.6 and ((pred==c)&(tid==t)).sum()/(pred==c).sum()>=.6
                for c in set(pred[(tid==t)&(pred>=0)])) for t in tracks)
    kp={p:[int((pred[part==p]>=0).sum()),int((part==p).sum())] for p in ['proton','17C_recoil']}
    return int(merged),len(tracks),rec,kp

D=pd.read_parquet("data/sim_noisy.parquet")
g=D.gid.unique(); np.random.shuffle(g); n=len(g); te_g=g[int(0.8*n)+int(0.1*n):]
net=OCNet(LAT).to(DEV); net.load_state_dict(torch.load("data/oc_net.pt",map_location=DEV)); net.eval()

def run(fn):
    M=0;rec=0;ntr=0;nc=[];kp={'proton':[0,0],'17C_recoil':[0,0]}
    for gid in te_g:
        gg=D[D.gid==gid]; x,y,z,q=(gg[c].to_numpy() for c in ['x','y','z','q']); part=gg['particle'].to_numpy()
        tid=gg['label'].to_numpy().copy(); tid[part=='noise']=-1
        if len(x)<14: continue
        pred=fn(x,y,z,q)
        m,nt,r,kpp=metrics(pred,tid,part); M+=m;rec+=r;ntr+=nt;nc.append(len(set(pred[pred>=0])))
        for p in kp: kp[p][0]+=kpp[p][0]; kp[p][1]+=kpp[p][1]
    N=sum(1 for gid in te_g if (D.gid==gid).sum()>=14)
    return dict(merge=100*M/N, rec=100*rec/max(ntr,1), clus=np.mean(nc),
                proton=100*kp['proton'][0]/max(kp['proton'][1],1), c17=100*kp['17C_recoil'][0]/max(kp['17C_recoil'][1],1))

from sklearn.cluster import DBSCAN
OC_EPS = 0.05
def oc_fn(x,y,z,q):
    feat=np.stack([x/250.,y/250.,(z-600.)/450.,np.log1p(np.clip(q,0,None))/9.],1)
    with torch.no_grad():
        co,be=net(torch.tensor(feat,dtype=torch.float,device=DEV),torch.zeros(len(x),dtype=torch.long,device=DEV))
    return DBSCAN(eps=OC_EPS, min_samples=4).fit_predict(co.cpu().numpy())  # cluster the latent coords
def dc_fn(x,y,z,q): return dircluster(np.stack([x,y,z],1), q, qratio=0.65, min_hits=4)

OC=run(oc_fn); DC=run(dc_fn)
print(f"held-out sim events: {len(te_g)}")
print(f"{'method':<14}{'merge%':>8}{'recovery%':>10}{'proton%':>9}{'17C%':>7}{'clus/ev':>9}")
print(f"{'ObjCondense':<14}{OC['merge']:>7.1f}%{OC['rec']:>9.1f}%{OC['proton']:>8.1f}%{OC['c17']:>6.1f}%{OC['clus']:>9.2f}")
print(f"{'dircluster':<14}{DC['merge']:>7.1f}%{DC['rec']:>9.1f}%{DC['proton']:>8.1f}%{DC['c17']:>6.1f}%{DC['clus']:>9.2f}")
