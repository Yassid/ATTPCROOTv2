#!/usr/bin/env python
"""Train DGCNN track-separator on the NOISY SIMULATION (perfect truth labels).

Discriminative instance loss; each per-event label (proton, 17C, noise-haze) is an
instance. Noise is its own instance -> repelled from tracks -> dropped at inference.
Split events randomly (sim has no runs). Eval on held-out sim (ARI vs truth).
  ~/gnn_env/bin/python train_sim.py [parquet] [N_TRAIN] [EPOCHS]
"""
import sys, numpy as np, pandas as pd, json, time
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score
torch.manual_seed(0); np.random.seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
IN_DIM, EMB_DIM, K = 4, 4, 20
DELTA_V, DELTA_D = 0.5, 2.5
PARQ   = sys.argv[1] if len(sys.argv) > 1 else "data/sim_noisy.parquet"
N_TRAIN = int(sys.argv[2]) if len(sys.argv) > 2 else 100000
EPOCHS  = int(sys.argv[3]) if len(sys.argv) > 3 else 30
BATCH, MAX_HITS = 32, 800

def to_data(df):
    out = []
    for gid, g in df.groupby('gid'):
        x = g['x'].to_numpy(); y = g['y'].to_numpy(); z = g['z'].to_numpy(); q = g['q'].to_numpy()
        lab = pd.factorize(g['label'].to_numpy())[0]
        if len(x) < 10 or lab.max() < 0: continue
        if len(x) > MAX_HITS:
            idx = np.random.choice(len(x), MAX_HITS, replace=False)
            x, y, z, q, lab = x[idx], y[idx], z[idx], q[idx], lab[idx]
        feat = np.stack([x/250.0, y/250.0, (z-600.0)/450.0, np.log1p(np.clip(q,0,None))/9.0], 1)
        out.append(Data(x=torch.tensor(feat, dtype=torch.float),
                        y=torch.tensor(lab, dtype=torch.long),
                        pos=torch.tensor(np.stack([x,y,z],1), dtype=torch.float)))
    return out

print(f"[{DEV}] loading {PARQ}")
D = pd.read_parquet(PARQ)
gids = D.gid.unique(); np.random.shuffle(gids)
gids = gids[:N_TRAIN]
n = len(gids); ntr = int(0.8*n); nva = int(0.1*n)
tr_g, va_g, te_g = gids[:ntr], gids[ntr:ntr+nva], gids[ntr+nva:]
train = to_data(D[D.gid.isin(tr_g)]); val = to_data(D[D.gid.isin(va_g)]); test = to_data(D[D.gid.isin(te_g)])
print(f"events: {len(train)} train / {len(val)} val / {len(test)} test")

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

def disc_loss(emb, y, batch):
    Lvar = Ldist = Lreg = 0.0; nev = 0
    for b in batch.unique():
        m = batch == b; e, lab = emb[m], y[m]; means = []
        for c in lab.unique():
            cm = lab == c; mu = e[cm].mean(0); means.append(mu)
            Lvar += torch.clamp(torch.norm(e[cm]-mu, dim=1)-DELTA_V, min=0).pow(2).mean()
        means = torch.stack(means); Kk = means.shape[0]
        if Kk > 1:
            dm = torch.cdist(means, means); idx = torch.triu_indices(Kk, Kk, 1)
            Ldist += torch.clamp(2*DELTA_D-dm[idx[0],idx[1]], min=0).pow(2).mean()
        Lreg += torch.norm(means, dim=1).mean(); nev += 1
    return (Lvar+Ldist+0.001*Lreg)/max(nev,1)

net = Net(K).to(DEV)
opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-5)
sch = torch.optim.lr_scheduler.StepLR(opt, step_size=12, gamma=0.5)
tl = DataLoader(train, batch_size=BATCH, shuffle=True)
vl = DataLoader(val, batch_size=BATCH)
t0 = time.time(); best = 1e9
for ep in range(EPOCHS):
    net.train(); tr = 0
    for d in tl:
        d = d.to(DEV); opt.zero_grad(); loss = disc_loss(net(d.x, d.batch), d.y, d.batch)
        loss.backward(); opt.step(); tr += loss.item()
    sch.step()
    net.eval(); vr = 0
    with torch.no_grad():
        for d in vl: d = d.to(DEV); vr += disc_loss(net(d.x, d.batch), d.y, d.batch).item()
    vr /= max(len(vl), 1)
    if vr < best: best = vr; torch.save(net.state_dict(), "data/embed_sim_noisy.pt")
    if ep % 3 == 0 or ep == EPOCHS-1:
        print(f"ep {ep:2d} train {tr/len(tl):.3f} val {vr:.3f} best {best:.3f} ({time.time()-t0:.0f}s)")
net.load_state_dict(torch.load("data/embed_sim_noisy.pt", map_location=DEV)); net.eval()

# held-out sim ARI vs truth (sweep DBSCAN eps)
aris = {e: [] for e in [0.3, 0.5, 0.7, 0.9, 1.1]}
embs = []
with torch.no_grad():
    for d in test:
        embs.append(net(d.x.to(DEV), torch.zeros(d.x.shape[0], dtype=torch.long, device=DEV)).cpu().numpy())
for eps in aris:
    for i, d in enumerate(test):
        pr = DBSCAN(eps=eps, min_samples=6).fit_predict(embs[i])
        aris[eps].append(adjusted_rand_score(d.y.numpy(), pr))
best_eps = max(aris, key=lambda e: np.mean(aris[e])); best_ari = np.mean(aris[best_eps])
print(f"\n=== HELD-OUT SIM ({len(test)} events) ===")
print(f"  best eps {best_eps}: ARI(pred vs truth) = {best_ari:.3f}")
json.dump({'ari_sim': best_ari, 'best_eps': best_eps, 'n_train': len(train), 'epochs': EPOCHS},
          open("data/metrics_sim_noisy.json", "w"), indent=1)
print("saved data/embed_sim_noisy.pt + data/metrics_sim_noisy.json")
