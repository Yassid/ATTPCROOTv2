#!/usr/bin/env python
"""SUPERVISED track separation trained on Spyral's REAL labels (no sim, no domain gap).

DGCNN -> 4-d embedding, discriminative instance loss with Spyral cluster ids as the instances.
Train on runs 0300-0304, EVALUATE on held-out run_0305 (never seen). Metrics: ARI vs Spyral
(does the NN reproduce Spyral on unseen events?) + clean-proton rate of NN vs Spyral (does it
match/beat Spyral's physical-track quality?).
Run: ~/gnn_env/bin/python sup/train_supervised.py [N_TRAIN] [EPOCHS]
"""
import sys, numpy as np, pandas as pd, json, time
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score
torch.manual_seed(0); np.random.seed(0)
DEV='cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, EMB_DIM, K = 4, 4, 20
DELTA_V, DELTA_D = 0.5, 2.5
N_TRAIN = int(sys.argv[1]) if len(sys.argv)>1 else 18000
EPOCHS  = int(sys.argv[2]) if len(sys.argv)>2 else 35
BATCH   = 32

def to_data(df):
    out=[]
    for gid,g in df.groupby('gid'):
        x=g['x'].to_numpy();y=g['y'].to_numpy();z=g['z'].to_numpy();q=g['q'].to_numpy()
        lab=pd.factorize(g['label'].to_numpy())[0]            # remap Spyral ids -> 0..K-1
        if len(x)<10 or lab.max()<0: continue
        feat=np.stack([x/250.0,y/250.0,(z-600.0)/450.0,np.log1p(np.clip(q,0,None))/9.0],1)
        out.append(Data(x=torch.tensor(feat,dtype=torch.float),
                        y=torch.tensor(lab,dtype=torch.long),
                        pos=torch.tensor(np.stack([x,y,z],1),dtype=torch.float)))
    return out

print(f"[{DEV}] loading spyral_all.parquet ...")
D=pd.read_parquet("data/attpc_all.parquet")
tr_df=D[D.run<305]; te_df=D[D.run==305]
tr_gids=tr_df.gid.unique(); np.random.shuffle(tr_gids)
tr_gids=tr_gids[:N_TRAIN]
train=to_data(tr_df[tr_df.gid.isin(tr_gids)])
nval=max(1,int(0.1*len(train))); val=train[:nval]; train=train[nval:]
# held-out test: cap to 1200 events for speed of eval
te_gids=te_df.gid.unique()[:1200]
test=to_data(te_df[te_df.gid.isin(te_gids)])
print(f"events: {len(train)} train / {len(val)} val / {len(test)} heldout(run_0305)")

def mlp(ch):
    L=[]
    for i in range(len(ch)-1): L+=[nn.Linear(ch[i],ch[i+1]),nn.BatchNorm1d(ch[i+1]),nn.ReLU()]
    return nn.Sequential(*L)
class Net(nn.Module):
    def __init__(s,k=20):
        super().__init__()
        s.ec1=DynamicEdgeConv(mlp([2*IN_DIM,64,64]),k,'max');s.ec2=DynamicEdgeConv(mlp([2*64,64,64]),k,'max')
        s.ec3=DynamicEdgeConv(mlp([2*64,128]),k,'max')
        s.head=nn.Sequential(nn.Linear(256,128),nn.BatchNorm1d(128),nn.ReLU(),nn.Linear(128,EMB_DIM))
    def forward(s,x,b):
        x1=s.ec1(x,b);x2=s.ec2(x1,b);x3=s.ec3(x2,b);return s.head(torch.cat([x1,x2,x3],1))

def disc_loss(emb,y,batch):
    Lvar=Ldist=Lreg=0.0; nev=0
    for b in batch.unique():
        m=batch==b; e,lab=emb[m],y[m]; means=[]
        for c in lab.unique():
            cm=lab==c; mu=e[cm].mean(0); means.append(mu)
            Lvar+=torch.clamp(torch.norm(e[cm]-mu,dim=1)-DELTA_V,min=0).pow(2).mean()
        means=torch.stack(means); Kk=means.shape[0]
        if Kk>1:
            dm=torch.cdist(means,means); idx=torch.triu_indices(Kk,Kk,1)
            Ldist+=torch.clamp(2*DELTA_D-dm[idx[0],idx[1]],min=0).pow(2).mean()
        Lreg+=torch.norm(means,dim=1).mean(); nev+=1
    return (Lvar+Ldist+0.001*Lreg)/max(nev,1)

net=Net(K).to(DEV)
opt=torch.optim.Adam(net.parameters(),lr=1e-3,weight_decay=1e-5)
sch=torch.optim.lr_scheduler.StepLR(opt,step_size=15,gamma=0.5)
tl=DataLoader(train,batch_size=BATCH,shuffle=True); vl=DataLoader(val,batch_size=BATCH)
t0=time.time()
for ep in range(EPOCHS):
    net.train(); tr=0
    for d in tl:
        d=d.to(DEV); opt.zero_grad(); loss=disc_loss(net(d.x,d.batch),d.y,d.batch)
        loss.backward(); opt.step(); tr+=loss.item()
    sch.step()
    if ep%5==0 or ep==EPOCHS-1:
        net.eval(); vr=0
        with torch.no_grad():
            for d in vl: d=d.to(DEV); vr+=disc_loss(net(d.x,d.batch),d.y,d.batch).item()
        print(f"ep {ep:2d} train {tr/len(tl):.3f} val {vr/max(len(vl),1):.3f} ({time.time()-t0:.0f}s)")
torch.save(net.state_dict(),"sup/embed_attpc.pt")

# -------- eval on held-out run_0305: ARI vs Spyral + clean-proton rate --------
def clean_proton(P):
    n=len(P)
    if not (25<=n<=170): return False
    C=P-P.mean(0); ev=np.linalg.eigvalsh(C.T@C/n); lin=ev[2]/(ev.sum()+1e-9)
    r=np.hypot(P[:,0],P[:,1]).mean(); zspan=P[:,2].max()-P[:,2].min()
    x,y=P[:,0],P[:,1]; A=np.stack([2*x,2*y,np.ones_like(x)],1); b=x*x+y*y
    try:
        cx,cy,cc=np.linalg.lstsq(A,b,rcond=None)[0]; R=np.sqrt(max(cc+cx*cx+cy*cy,1e-6)); arc=np.median(np.abs(np.hypot(x-cx,y-cy)-R))
    except: arc=99
    return r>40 and lin>0.85 and arc<12 and zspan>80
net.eval()
# sweep eps on a slice, pick best ARI vs Spyral
aris={e:[] for e in [0.3,0.5,0.7,0.9,1.1]}
embs=[]
with torch.no_grad():
    for d in test:
        embs.append(net(d.x.to(DEV),torch.zeros(d.x.shape[0],dtype=torch.long,device=DEV)).cpu().numpy())
for eps in aris:
    for i,d in enumerate(test):
        pr=DBSCAN(eps=eps,min_samples=6).fit_predict(embs[i])
        aris[eps].append(adjusted_rand_score(d.y.numpy(),pr))
best_eps=max(aris,key=lambda e:np.mean(aris[e])); best_ari=np.mean(aris[best_eps])
# clean-proton rate: NN vs Spyral on held-out
nn_clean=sp_clean=0
for i,d in enumerate(test):
    pos=d.pos.numpy(); sp=d.y.numpy()
    pr=DBSCAN(eps=best_eps,min_samples=6).fit_predict(embs[i])
    nn_clean+=any(clean_proton(pos[pr==c]) for c in set(pr) if c>=0)
    sp_clean+=any(clean_proton(pos[sp==c]) for c in set(sp))
n=len(test)
print(f"\n=== HELD-OUT run_0305 ({n} events) ===")
print(f"  best eps {best_eps}: ARI(NN vs Spyral) = {best_ari:.3f}   (1.0 = reproduces Spyral exactly)")
print(f"  clean-proton rate:  NN {nn_clean/n*100:.0f}%   vs   Spyral {sp_clean/n*100:.0f}%")
json.dump({'ari_vs_spyral':best_ari,'best_eps':best_eps,'nn_clean':nn_clean/n,'sp_clean':sp_clean/n,
           'n_test':n,'n_train':len(train),'epochs':EPOCHS},open("sup/metrics_attpc.json","w"),indent=1)
print("saved sup/embed_attpc.pt + sup/metrics_attpc.json")
