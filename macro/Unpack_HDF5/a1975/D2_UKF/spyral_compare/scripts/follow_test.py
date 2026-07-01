#!/usr/bin/env python
import sys, numpy as np, polars as pl
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "scripts")
from follower import follow_event

H = pl.read_csv("/home/yassid/spyral_d2/follower_hits.csv")
SPY = "/home/yassid/spyral_d2/spyral_clouds.root"
import uproot
sp = uproot.open(SPY)['spyral'].arrays(['event','kind','lab'], library='np')

def spyral_maxcluster(ev):
    m = (sp['event'] == ev) & (sp['kind'] == 0)
    if not m.any(): return 0
    labs, cnts = np.unique(sp['lab'][m], return_counts=True)
    return int(cnts.max())

def ev_arrays(ev):
    d = H.filter(pl.col('event') == ev)
    return dict(x=d['x'].to_numpy(), y=d['y'].to_numpy(), z=d['z'].to_numpy(),
                praTrk=d['praTrk'].to_numpy())

def plot(ev, ax_xy, ax_zy):
    a = ev_arrays(ev); r = follow_event(a)
    x, y, z, pra = a['x'], a['y'], a['z'], a['praTrk']
    seed = np.array(r.get('seed_idx', []), dtype=int)
    grab = np.array(r.get('grabbed_idx', []), dtype=int)
    other = np.ones(len(x), bool); other[seed] = False
    if len(grab): other[grab] = False
    for ax, (h, v) in [(ax_xy, (x, y)), (ax_zy, (z, y))]:
        ax.scatter(h[other], v[other], s=3, c='lightgray')
        ax.scatter(h[seed], v[seed], s=8, c='tab:blue', label=f'PRA seed ({r["seed"]})')
        if len(grab):
            ax.scatter(h[grab], v[grab], s=8, c='red', label=f'KF recovered (+{r["grabbed"]})')
    smax = spyral_maxcluster(ev)
    ax_xy.set_title(f"ev {ev}: PRA {r['seed']} -> KF {r['extended']}  (Spyral {smax})")
    ax_xy.legend(fontsize=7); ax_xy.set_xlabel('x'); ax_xy.set_ylabel('y')
    ax_zy.set_xlabel('z drift'); ax_zy.set_ylabel('y')
    return r['seed'], r['extended'], smax

evs = [int(e) for e in (sys.argv[1:] or [49, 57, 15, 11090])]
fig, axes = plt.subplots(len(evs), 2, figsize=(13, 4.2*len(evs)))
if len(evs) == 1: axes = axes.reshape(1, 2)
for i, ev in enumerate(evs):
    s, e, sm = plot(ev, axes[i, 0], axes[i, 1])
    print(f"ev {ev}: seed={s} extended={e} spyral={sm}  recovered_frac PRA={s/max(sm,1):.2f} -> KF={e/max(sm,1):.2f}")
plt.tight_layout()
out = "plots/follower_cases.png"
plt.savefig(out, dpi=105); print("saved", out)
