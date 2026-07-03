#!/usr/bin/env python
"""Train Object Condensation on the REAL gold labels (no sim->real gap), 5-fold CV.
Tests whether real-trained OC separates real tracks (vs sim-trained OC which merged 57-74%,
vs geometric dircluster 1.4%). Only 73 gold events -> CV + heavy augmentation caution.
  ~/gnn_env/bin/python train_oc_gold.py [EPOCHS] [LATENT]
"""
import sys, numpy as np, pandas as pd, time
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure
torch.manual_seed(0); np.random.seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
EPOCHS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
LAT = int(sys.argv[2]) if len(sys.argv) > 2 else 3
IN_DIM, K, Q_MIN, S_REP, S_BETA = 4, 20, 0.34, 0.6, 0.1

G = pd.read_parquet("labeling/data/labels.parquet"); G = G[G.reviewed]
evs = [e for e in sorted(G.event.unique()) if (G.event == e).sum() >= 14]
def ev_data(ev):
    g = G[G.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    lab = g['label'].to_numpy()                      # gold tracks 0/1/2 (no explicit noise)
    feat = np.stack([x/250., y/250., (z-600.)/450., np.log1p(np.clip(q,0,None))/9.], 1)
    return Data(x=torch.tensor(feat, dtype=torch.float), y=torch.tensor(lab, dtype=torch.long),
                pos=torch.tensor(np.stack([x,y,z],1), dtype=torch.float), gid=int(ev))

def mlp(ch):
    L=[]
    for i in range(len(ch)-1): L+=[nn.Linear(ch[i],ch[i+1]),nn.BatchNorm1d(ch[i+1]),nn.ReLU()]
    return nn.Sequential(*L)
class OCNet(nn.Module):
    def __init__(s,lat,k=20):
        super().__init__(); s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max');s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max'); s.trunk=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU())
        s.coord=nn.Linear(128,lat); s.beta=nn.Linear(128,1)
    def forward(s,x,b):
        x1=s.ec1(x,b);x2=s.ec2(x1,b);x3=s.ec3(x2,b);h=s.trunk(torch.cat([x1,x2,x3],1)); return s.coord(h),torch.sigmoid(s.beta(h)).squeeze(-1)

def oc_loss(coord, beta, y, batch):
    q = torch.atanh(beta.clamp(1e-6,1-1e-6))**2 + Q_MIN; Lv=Lb=0.0; nev=0
    for b in batch.unique():
        m=batch==b; co,be,la,qq=coord[m],beta[m],y[m],q[m]; N=co.shape[0]
        objs=[int(t) for t in la.unique().tolist() if t>=0]; Lv_ev=0.0
        for t in objs:
            tm=la==t; ti=torch.where(tm)[0]; cp=ti[torch.argmax(be[ti])]; qcp,ccp=qq[cp],co[cp]
            Lv_ev=Lv_ev+(qq[tm]*qcp*((co[tm]-ccp)**2).sum(1)).sum()
            Lv_ev=Lv_ev+S_REP*(qq[~tm]*qcp*torch.clamp(1-torch.norm(co[~tm]-ccp,dim=1),min=0)).sum()
        Lv=Lv+Lv_ev/max(N,1)
        Lb=Lb+sum((1-be[torch.where(la==t)[0][torch.argmax(be[torch.where(la==t)[0]])]]) for t in objs)/max(len(objs),1)
        nev+=1
    return (Lv+S_BETA*Lb)/max(nev,1)

MH=8
def metrics(pred, gold):
    tr=[t for t in set(gold) if (gold==t).sum()>=MH]
    mg=any(len([t for t in tr if (gold[pred==c]==t).sum()>=MH])>=2 for c in set(pred[pred>=0]))
    return int(mg), adjusted_rand_score(gold,pred), homogeneity_completeness_v_measure(gold,pred)[0]

# 5-fold CV
np.random.seed(0); order=np.array(evs); np.random.shuffle(order); folds=np.array_split(order,5)
all_pred=[]  # (coords, gold) for held-out, all folds
t0=time.time()
for fi in range(5):
    te_ev=set(folds[fi].tolist()); tr_ev=[e for e in evs if e not in te_ev]
    tr=[ev_data(e) for e in tr_ev]; te=[ev_data(e) for e in folds[fi]]
    net=OCNet(LAT).to(DEV); opt=torch.optim.Adam(net.parameters(),lr=1e-3,weight_decay=1e-4)
    tl=DataLoader(tr,batch_size=8,shuffle=True)
    for ep in range(EPOCHS):
        net.train()
        for d in tl: d=d.to(DEV); opt.zero_grad(); co,be=net(d.x,d.batch); l=oc_loss(co,be,d.y,d.batch); l.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        for d in te:
            co,_=net(d.x.to(DEV),torch.zeros(d.x.shape[0],dtype=torch.long,device=DEV))
            all_pred.append((co.cpu().numpy(), d.y.numpy()))
    print(f"fold {fi} done ({time.time()-t0:.0f}s)")

print(f"\nheld-out gold events (5-fold CV): {len(all_pred)}")
print(f"{'eps':>6}{'merge%':>8}{'ARI':>8}{'homog':>8}")
for eps in [0.02,0.03,0.05,0.08,0.12]:
    M=[];A=[];H=[]
    for co,gold in all_pred:
        pred=DBSCAN(eps=eps,min_samples=4).fit_predict(co); m,a,h=metrics(pred,gold); M.append(m);A.append(a);H.append(h)
    print(f"{eps:>6}{100*np.mean(M):>7.1f}%{np.mean(A):>8.3f}{np.mean(H):>8.3f}")
print("compare: sim-trained OC merge 57-74%, dircluster merge 1.4% ARI0.45 homog0.73")
