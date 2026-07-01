#!/usr/bin/env python
"""Final composition verification: sim-only -> overlay-cc -> overlay-v2 vs real.
Reports per-event feature KS to real-rxn (reaction-like) and real-all, plus the
charge giveaway KS within each overlay. Saves a summary PNG."""
import polars as pl, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp
SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
real = pl.read_csv(f"{SC}/data/real_events.csv").rename({"q": "charge"})
sets = {
    'sim-only':   pl.read_parquet(f"{SC}/data/train.parquet"),
    'overlay-cc': pl.read_parquet(f"{SC}/data/overlay_cc.parquet"),
    'overlay-v2': pl.read_parquet(f"{SC}/data/overlay_v2.parquet"),
}
rr = real.with_columns((pl.col('x')**2 + pl.col('y')**2).sqrt().alias('r'))
rxn = []
for ev, g in rr.group_by('event'):
    off = g.filter(pl.col('r') >= 35.0)
    if off.height >= 5 and cKDTree(off.select(['x','y','z']).to_numpy()).query_ball_point(
            off.select(['x','y','z']).to_numpy(), r=12.0, return_length=True).max() >= 7:
        rxn.append(int(ev[0]))
real_rxn = real.filter(pl.col('event').is_in(rxn))

def pe(df):
    return df.group_by('event').agg([pl.len().alias('nhits'),
        (pl.col('z').max()-pl.col('z').min()).alias('zspan'),
        ((pl.col('x')**2+pl.col('y')**2).sqrt().max()).alias('rmax')])
E = {k: pe(v) for k, v in sets.items()}
E['real-rxn'] = pe(real_rxn); E['real-all'] = pe(real)
def ks(a, b): return ks_2samp(a, b).statistic

for refname in ['real-rxn', 'real-all']:
    ref = E[refname]
    print(f"\n=== COMPOSITION KS vs {refname} (med nhits {np.median(ref['nhits'].to_numpy()):.0f}, "
          f"zspan {np.median(ref['zspan'].to_numpy()):.0f}, rmax {np.median(ref['rmax'].to_numpy()):.0f}) ===")
    print(f"{'feature':7s}  {'sim-only':>8s} {'overlay-cc':>11s} {'overlay-v2':>11s}")
    for f in ['nhits', 'zspan', 'rmax']:
        r = ref[f].to_numpy()
        print(f"{f:7s}  {ks(E['sim-only'][f].to_numpy(),r):8.2f} "
              f"{ks(E['overlay-cc'][f].to_numpy(),r):11.2f} {ks(E['overlay-v2'][f].to_numpy(),r):11.2f}")

def giveaway(df):
    q = df.with_columns([pl.col('charge').rank().over('event').alias('rk'),
                         pl.len().over('event').alias('nn')]).with_columns(
                         ((pl.col('rk')-0.5)/pl.col('nn')).alias('qr'))
    return ks(q.filter(pl.col('label')<2)['qr'].to_numpy(), q.filter(pl.col('label')==2)['qr'].to_numpy())
print("\n=== CHARGE GIVEAWAY (KS sim-vs-noise qrank; lower=better) ===")
print(f"  overlay-cc: {giveaway(sets['overlay-cc']):.2f}   overlay-v2: {giveaway(sets['overlay-v2']):.2f}")

# plot nhits + zspan progression
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for i, f in enumerate(['nhits', 'zspan']):
    bins = np.linspace(0, np.percentile(E['real-rxn'][f].to_numpy(), 97), 45)
    for n, c in [('sim-only','C1'),('overlay-cc','C2'),('overlay-v2','C3'),('real-rxn','C0')]:
        ax[i].hist(E[n][f].to_numpy(), bins=bins, density=True, histtype='step', lw=2, label=n, color=c)
    ax[i].set_title(f); ax[i].set_xlabel(f); ax[i].legend()
plt.tight_layout(); out = f"{SC}/overlay_v2_verify.png"; plt.savefig(out, dpi=110); print(f"\nsaved {out}")
