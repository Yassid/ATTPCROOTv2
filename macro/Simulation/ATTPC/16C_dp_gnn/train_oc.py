#!/usr/bin/env python
"""Object Condensation (arXiv:2309.16754 / Kieseler 2020) for AT-TPC track separation.

Per hit: condensation strength beta + clustering coords c (learned latent space). Tracks are
condensed around their max-beta hit (condensation point). Attractive pull to own CP + repulsive
push from other CPs separates overlapping tracks in latent space; beta term suppresses noise.
Trained on noisy-sim truth (proton/17C tracks + noise). Inference: greedy CP assignment.
  ~/gnn_env/bin/python train_oc.py [parquet] [EPOCHS] [LATENT]
"""
import sys, numpy as np, pandas as pd, json, time
import torch, torch.nn as nn
from torch_geometric.nn import DynamicEdgeConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
torch.manual_seed(0); np.random.seed(0)
DEV = 'cuda' if torch.cuda.is_available() else 'cpu'
PARQ   = sys.argv[1] if len(sys.argv) > 1 else "data/sim_noisy.parquet"
EPOCHS = int(sys.argv[2]) if len(sys.argv) > 2 else 30
LAT    = int(sys.argv[3]) if len(sys.argv) > 3 else 3
IN_DIM, K, BATCH, MAX_HITS = 4, 20, 16, 1000
S_BETA = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
OUT    = sys.argv[5] if len(sys.argv) > 5 else "data/oc_net.pt"
Q_MIN, S_REP, S_B = 0.34, 0.6, 0.09

def to_data(df):
    out = []
    for gid, g in df.groupby('gid'):
        x = g['x'].to_numpy(); y = g['y'].to_numpy(); z = g['z'].to_numpy(); q = g['q'].to_numpy()
        lab = g['label'].to_numpy().copy(); lab[g['particle'].to_numpy() == 'noise'] = -1  # noise = -1
        if len(x) < 12: continue
        if len(x) > MAX_HITS:
            idx = np.random.choice(len(x), MAX_HITS, replace=False)
            x, y, z, q, lab = x[idx], y[idx], z[idx], q[idx], lab[idx]
        feat = np.stack([x/250., y/250., (z-600.)/450., np.log1p(np.clip(q,0,None))/9.], 1)
        out.append(Data(x=torch.tensor(feat, dtype=torch.float), y=torch.tensor(lab, dtype=torch.long),
                        pos=torch.tensor(np.stack([x,y,z],1), dtype=torch.float)))
    return out

def mlp(ch):
    L = []
    for i in range(len(ch)-1): L += [nn.Linear(ch[i],ch[i+1]), nn.BatchNorm1d(ch[i+1]), nn.ReLU()]
    return nn.Sequential(*L)
class OCNet(nn.Module):
    def __init__(s, lat=3, k=20):
        super().__init__()
        s.ec1 = DynamicEdgeConv(mlp([2*IN_DIM,64,64]), k, 'max'); s.ec2 = DynamicEdgeConv(mlp([2*64,64,64]), k, 'max')
        s.ec3 = DynamicEdgeConv(mlp([2*64,128]), k, 'max')
        s.trunk = nn.Sequential(nn.Linear(256,128), nn.BatchNorm1d(128), nn.ReLU())
        s.coord = nn.Linear(128, lat)          # clustering coordinates
        s.beta = nn.Linear(128, 1)             # condensation strength (logit)
    def forward(s, x, b):
        x1 = s.ec1(x, b); x2 = s.ec2(x1, b); x3 = s.ec3(x2, b)
        h = s.trunk(torch.cat([x1, x2, x3], 1))
        return s.coord(h), torch.sigmoid(s.beta(h)).squeeze(-1)

def oc_loss(coord, beta, y, batch):
    beta_c = beta.clamp(1e-6, 1-1e-6)
    q = torch.atanh(beta_c)**2 + Q_MIN
    Lv = Lb = 0.0; nev = 0
    for b in batch.unique():
        m = batch == b; co, be, la, qq = coord[m], beta[m], y[m], q[m]
        N = co.shape[0]
        objs = [int(t) for t in la.unique().tolist() if t >= 0]
        noise = la < 0
        Lv_ev = 0.0
        for t in objs:
            tm = la == t
            ti = torch.where(tm)[0]
            cp = ti[torch.argmax(be[ti])]          # condensation point = max-beta hit in track
            qcp, ccp = qq[cp], co[cp]
            datt = ((co[tm] - ccp)**2).sum(1)                              # attract own hits
            Lv_ev = Lv_ev + (qq[tm] * qcp * datt).sum()
            drep = torch.clamp(1 - torch.norm(co[~tm] - ccp, dim=1), min=0)  # repel others
            Lv_ev = Lv_ev + S_REP * (qq[~tm] * qcp * drep).sum()
        Lv = Lv + Lv_ev / max(N, 1)
        # beta: CP->1, noise beta suppressed
        Lb_cp = sum((1 - be[torch.where(la == t)[0][torch.argmax(be[torch.where(la==t)[0]])]]) for t in objs)
        Lb_cp = Lb_cp / max(len(objs), 1) if objs else torch.tensor(0.0, device=DEV)
        Lb_noise = S_B * be[noise].sum() / (noise.sum() + 1e-9) if noise.any() else 0.0
        Lb = Lb + Lb_cp + Lb_noise
        nev += 1
    return (Lv + S_BETA * Lb) / max(nev, 1)

@torch.no_grad()
def cluster_oc(coord, beta, t_beta=0.1, t_d=0.5):
    n = len(beta); lab = -np.ones(n, int); un = np.ones(n, bool)
    order = np.argsort(-beta); cid = 0
    C = coord
    for i in order:
        if beta[i] < t_beta: break
        if not un[i]: continue
        d = np.linalg.norm(C - C[i], axis=1)
        sel = un & (d < t_d)
        lab[sel] = cid; un[sel] = False; cid += 1
    return lab

print(f"[{DEV}] loading {PARQ}  latent={LAT}")
D = pd.read_parquet(PARQ)
g = D.gid.unique(); np.random.shuffle(g); n = len(g); ntr, nva = int(.8*n), int(.1*n)
tr = to_data(D[D.gid.isin(g[:ntr])]); va = to_data(D[D.gid.isin(g[ntr:ntr+nva])]); te = to_data(D[D.gid.isin(g[ntr+nva:])])
print(f"events: {len(tr)}/{len(va)}/{len(te)}")
net = OCNet(LAT, K).to(DEV)
opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-5)
sch = torch.optim.lr_scheduler.StepLR(opt, 12, 0.5)
tl = DataLoader(tr, batch_size=BATCH, shuffle=True); vl = DataLoader(va, batch_size=BATCH)
t0 = time.time(); best = 1e9
for ep in range(EPOCHS):
    net.train(); s = 0
    for d in tl:
        d = d.to(DEV); opt.zero_grad(); co, be = net(d.x, d.batch)
        loss = oc_loss(co, be, d.y, d.batch); loss.backward(); opt.step(); s += loss.item()
    sch.step(); net.eval(); vs = 0
    with torch.no_grad():
        for d in vl: d = d.to(DEV); co, be = net(d.x, d.batch); vs += oc_loss(co, be, d.y, d.batch).item()
    vs /= max(len(vl), 1)
    if vs < best: best = vs; torch.save(net.state_dict(), OUT)
    if ep % 3 == 0 or ep == EPOCHS-1:
        print(f"ep {ep:2d} train {s/len(tl):.3f} val {vs:.3f} best {best:.3f} ({time.time()-t0:.0f}s)")
net.load_state_dict(torch.load(OUT, map_location=DEV)); net.eval()

# ---- eval on held-out sim: merge / recovery / per-particle (sweep t_d) ----
from sklearn.metrics import homogeneity_completeness_v_measure
def metrics(pred, tid, part, mh=8):
    tracks = [t for t in set(tid[tid>=0]) if (tid==t).sum()>=mh]
    merged = any(len([t for t in tracks if (tid[pred==c]==t).sum()>=mh])>=2 for c in set(pred[pred>=0]))
    rec = 0
    for t in tracks:
        tm = tid==t
        if any(((pred==c)&tm).sum()/tm.sum()>=.6 and ((pred==c)&tm).sum()/(pred==c).sum()>=.6 for c in set(pred[tm&(pred>=0)])): rec += 1
    keptp = {} if part is None else {p: [int((pred[part==p]>=0).sum()), int((part==p).sum())] for p in ['proton','17C_recoil']}
    return int(merged), len(tracks), rec, keptp
EV = []  # (coords, beta, truth_label, particle-array)
with torch.no_grad():
    for d in te:
        co, be = net(d.x.to(DEV), torch.zeros(d.x.shape[0], dtype=torch.long, device=DEV))
        EV.append((co.cpu().numpy(), be.cpu().numpy(), d.y.numpy()))
from sklearn.cluster import DBSCAN
print(f"\n{'eps':>6}{'merge%':>8}{'recov%':>8}{'clus/ev':>9}")
res = []
for eps in [0.04, 0.05, 0.06, 0.07, 0.08]:
    M=0; rec=0; ntr_=0; nc=[]
    for co, be, tid in EV:
        pred = DBSCAN(eps=eps, min_samples=4).fit_predict(co)
        m, nt, r, _ = metrics(pred, tid, None); M += m; rec += r; ntr_ += nt; nc.append(len(set(pred[pred>=0])))
    mr = 100*M/len(EV); rc = 100*rec/max(ntr_,1); res.append((eps, mr, rc))
    print(f"{eps:>6}{mr:>7.1f}%{rc:>7.1f}%{np.mean(nc):>9.2f}")
ok = [r for r in res if r[1] <= 3.0]; best = max(ok or res, key=lambda r: r[2])
print(f"SUMMARY lat={LAT} sbeta={S_BETA} out={OUT} :: eps={best[0]} merge={best[1]:.1f}% recovery={best[2]:.1f}%")
json.dump({'epochs':EPOCHS,'latent':LAT}, open("data/oc_metrics.json","w"))
print("saved data/oc_net.pt")
