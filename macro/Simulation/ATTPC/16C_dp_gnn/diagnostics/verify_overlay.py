#!/usr/bin/env python
"""Verify the composition gap is closed by the overlay, and the charge-consistency fix.
 (1) Composition: per-event n_hits / z-span / r_max progression  sim-only -> overlay -> real.
     Did adding real beam+noise pull the sim distribution onto real? (baseline KS ~0.28)
 (2) Charge giveaway: within an overlay event, is per-event charge-quantile able to tell
     sim-signal (label 0/1) from real beam/noise (label 2)? RAW overlay should leak
     (high KS), charge-matched (cc) should not (low KS).
Reference 'real_rxn' = real events WITH a dense off-axis track (reaction-like), the true
analog of an overlay event; 'real_all' = every real event (beam-only dominated).
"""
import polars as pl, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp

SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
real = pl.read_csv(f"{SC}/data/real_events.csv").rename({"q": "charge"})
simo = pl.read_parquet(f"{SC}/data/train.parquet")          # sim signal only
ovc  = pl.read_parquet(f"{SC}/data/overlay_cc.parquet")     # charge-matched
ovr  = pl.read_parquet(f"{SC}/data/overlay.parquet")        # raw charge

# ---- split real into reaction-like vs beam-only (same logic as the canvas selector) ----
R_OFFAXIS, DENS_R, DENS_THRESH = 35.0, 12.0, 7
rr = real.with_columns((pl.col('x')**2 + pl.col('y')**2).sqrt().alias('r'))
rxn_ids = []
for ev, g in rr.group_by('event'):
    off = g.filter(pl.col('r') >= R_OFFAXIS)
    if off.height >= 5:
        Q = off.select(['x', 'y', 'z']).to_numpy()
        if cKDTree(Q).query_ball_point(Q, r=DENS_R, return_length=True).max() >= DENS_THRESH:
            rxn_ids.append(int(ev[0]))
real_rxn = real.filter(pl.col('event').is_in(rxn_ids))
print(f"real reaction-like events: {len(rxn_ids)} / {real['event'].n_unique()}")

def per_event(df):
    return df.group_by('event').agg([
        pl.len().alias('nhits'),
        (pl.col('z').max() - pl.col('z').min()).alias('zspan'),
        ((pl.col('x')**2 + pl.col('y')**2).sqrt().max()).alias('rmax'),
    ])

E = {n: per_event(d) for n, d in
     [('sim-only', simo), ('overlay-cc', ovc), ('real-rxn', real_rxn), ('real-all', real)]}

def ks(a, b): return ks_2samp(a, b).statistic
print("\n=== COMPOSITION GAP (KS vs real-rxn ; lower = better) ===")
print(f"{'feature':8s}  {'sim-only':>22s}  {'overlay-cc':>22s}   KS(sim-only)  KS(overlay)")
ref = E['real-rxn']
for f in ['nhits', 'zspan', 'rmax']:
    so, oc, rk = E['sim-only'][f].to_numpy(), E['overlay-cc'][f].to_numpy(), ref[f].to_numpy()
    print(f"{f:8s}  med {np.median(so):7.0f} [{np.percentile(so,10):5.0f},{np.percentile(so,90):6.0f}]"
          f"  med {np.median(oc):7.0f} [{np.percentile(oc,10):5.0f},{np.percentile(oc,90):6.0f}]"
          f"     {ks(so,rk):.2f}         {ks(oc,rk):.2f}")
print(f"(real-rxn medians: nhits {np.median(ref['nhits'].to_numpy()):.0f}, "
      f"zspan {np.median(ref['zspan'].to_numpy()):.0f}, rmax {np.median(ref['rmax'].to_numpy()):.0f})")

# ---- (2) charge giveaway: per-event quantile rank, sim(0/1) vs noise(2) ----
def qrank(df):
    return df.with_columns([
        pl.col('charge').rank().over('event').alias('rk'),
        pl.len().over('event').alias('nn'),
    ]).with_columns(((pl.col('rk') - 0.5) / pl.col('nn')).alias('qr'))
print("\n=== CHARGE GIVEAWAY: KS( qrank[sim 0/1]  vs  qrank[noise 2] ) within overlay ===")
print("   high KS = per-event normalization SEPARATES sim from real noise (domain leak); low = fixed")
for name, ov in [('overlay RAW', ovr), ('overlay CC ', ovc)]:
    q = qrank(ov)
    qs = q.filter(pl.col('label') < 2)['qr'].to_numpy()
    qn = q.filter(pl.col('label') == 2)['qr'].to_numpy()
    print(f"  {name}:  KS = {ks(qs, qn):.3f}   (sim qrank med {np.median(qs):.2f}, noise med {np.median(qn):.2f})")

# ---- plots ----
fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))
for i, f in enumerate(['nhits', 'zspan']):
    bins = np.linspace(0, np.percentile(E['real-rxn'][f].to_numpy(), 98), 45)
    for n, c in [('sim-only', 'C1'), ('overlay-cc', 'C2'), ('real-rxn', 'C0')]:
        ax[i].hist(E[n][f].to_numpy(), bins=bins, density=True, alpha=.45, label=n, color=c)
    ax[i].set_title(f"{f} : sim-only -> overlay -> real"); ax[i].set_xlabel(f); ax[i].legend()
for i, (name, ov) in enumerate([('overlay RAW', ovr), ('overlay CC', ovc)]):
    q = qrank(ov)
    qs = q.filter(pl.col('label') < 2)['qr'].to_numpy()
    qn = q.filter(pl.col('label') == 2)['qr'].to_numpy()
    ax[2+i].hist(qs, bins=np.linspace(0,1,40), density=True, alpha=.5, label='sim 0/1', color='C3')
    ax[2+i].hist(qn, bins=np.linspace(0,1,40), density=True, alpha=.5, label='noise 2', color='C0')
    ax[2+i].set_title(f"{name}: per-event qrank by label\n(overlap=good)  KS={ks(qs,qn):.2f}")
    ax[2+i].set_xlabel('per-event charge quantile'); ax[2+i].legend()
plt.tight_layout()
out = f"{SC}/overlay_verify.png"; plt.savefig(out, dpi=110); print(f"\nsaved {out}")
