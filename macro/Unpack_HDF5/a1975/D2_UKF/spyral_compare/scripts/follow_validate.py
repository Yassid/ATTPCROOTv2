#!/usr/bin/env python
"""Run the grow+bridge follower over the validation set and report recovery by
category vs default PRA, plus over-grow on clean controls. Uses Spyral's largest
cluster as the (approximate) truth size for the main track."""
import sys, numpy as np, polars as pl, uproot
sys.path.insert(0, "scripts")
import follower

H = pl.read_csv("/home/yassid/spyral_d2/follower_hits.csv")
sp = uproot.open("/home/yassid/spyral_d2/spyral_clouds.root")['spyral'].arrays(['event','kind','lab'], library='np')
import json
SET = json.load(open("/home/yassid/spyral_d2/follower_events.json"))

# spyral max-cluster per event (vectorized: group by event,lab -> count -> max per event)
spdf = pl.DataFrame({'event': sp['event'], 'kind': sp['kind'], 'lab': sp['lab']}) \
    .filter(pl.col('kind') == 0) \
    .group_by(['event', 'lab']).len() \
    .group_by('event').agg(pl.col('len').max().alias('smax'))
spm = {int(e): int(s) for e, s in zip(spdf['event'], spdf['smax'])}

# partition once (fast) instead of 334 filters
H_ev = {int(k[0]): d for k, d in H.partition_by('event', as_dict=True).items()}

par = {}
if len(sys.argv) > 1:
    for kv in sys.argv[1:]:
        k, v = kv.split('='); par[k] = float(v)

cat_of = {}
for c, evs in SET['sets'].items():
    for e in evs: cat_of[int(e)] = c

rows = []
nproc = 0
for ev, d in H_ev.items():
    sm = spm.get(ev, 0)
    if sm < 30:
        continue
    if d.height > 3000:          # skip beam-dominated huge clouds
        continue
    a = dict(x=d['x'].to_numpy(), y=d['y'].to_numpy(), z=d['z'].to_numpy(), praTrk=d['praTrk'].to_numpy())
    r = follower.follow_event(a, par)
    nproc += 1
    if nproc % 50 == 0:
        print(f"  ...{nproc} events", flush=True)
    if not r.get('ok'):
        continue
    rows.append(dict(event=ev, cat=cat_of.get(ev, '?'), seed=r['seed'], ext=r['extended'],
                     spyral=sm, pra_frac=r['seed']/sm, kf_frac=r['extended']/sm))

df = pl.DataFrame(rows)
print(f"params: {par or 'defaults'}   events: {df.height}\n")
print(f"{'category':14s} {'n':>4s} {'PRA frac':>9s} {'KF frac':>9s} {'recovered':>10s} {'overgrow>1.3':>12s}")
for c in ['severe', 'partial', 'match', 'clean_single']:
    s = df.filter(pl.col('cat') == c)
    if s.height == 0: continue
    pf = s['pra_frac'].median(); kf = s['kf_frac'].median()
    rec = (s['kf_frac'] - s['pra_frac']).mean()
    og = (s['kf_frac'] > 1.3).mean()
    print(f"{c:14s} {s.height:4d} {pf:9.2f} {kf:9.2f} {rec:+10.2f} {100*og:11.0f}%")
# case studies
print("\ncase studies:")
for ev in [49, 57, 15, 11090, 36]:
    s = df.filter(pl.col('event') == ev)
    if s.height: r = s.row(0, named=True); print(f"  ev{ev:5d} [{r['cat']:12s}] PRA {r['pra_frac']:.2f} -> KF {r['kf_frac']:.2f}")
df.write_csv("/home/yassid/spyral_d2/follow_validate.csv")
