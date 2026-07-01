#!/usr/bin/env python
"""One-to-one Spyral vs ATTPCROOT-PRA clustering comparison for 16C(d,p) run_0016.
Focus: backward-track (theta>90 deg) identification. Event alignment is direct
(Spyral orig_event == PRA tree entry); angle conventions agree (median |dtheta|~1.3 deg).
"""
import polars as pl, glob, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/spyral_compare/plots"
BACK = 90.0          # backward threshold (deg)
PRA_MIN_HITS = 20    # substantial PRA track
PRA_BEAM_MAXTH = 165 # PRA tracks above this theta treated as beam-like (along -z)

# --- load Spyral ---
sp = pl.concat([pl.read_parquet(f) for f in sorted(
    glob.glob("/home/yassid/spyral_d2/workspace/Estimation/run_030*.parquet"))])
sp = sp.with_columns((pl.col('polar')*180/np.pi).alias('th'))

# --- load PRA ---
pra = pl.read_csv("/home/yassid/spyral_d2/pra_run0016.csv", infer_schema_length=0)
for c in ["geo_theta","geo_radius","zmean"]:
    pra = pra.with_columns(pl.col(c).cast(pl.Float64, strict=False))
for c in ["event","nhits"]:
    pra = pra.with_columns(pl.col(c).cast(pl.Int64, strict=False))
pra = pra.drop_nulls(subset=["geo_theta"])
pra_sub = pra.filter(pl.col('nhits') >= PRA_MIN_HITS)

# --- per-event aggregation ---
sp_ev = sp.group_by('orig_event').agg(
    pl.len().alias('n_sp'),
    (pl.col('th') > BACK).sum().alias('n_sp_back'),
    pl.col('th').max().alias('sp_thmax'),
).rename({'orig_event':'event'})

pra_ev = pra_sub.group_by('event').agg(
    pl.len().alias('n_pra'),
    ((pl.col('geo_theta') > BACK) & (pl.col('geo_theta') <= PRA_BEAM_MAXTH)).sum().alias('n_pra_back'),
    ((pl.col('geo_theta') > PRA_BEAM_MAXTH)).sum().alias('n_pra_beam'),
)

all_ev = sp_ev.join(pra_ev, on='event', how='full', coalesce=True).fill_null(0)
N = all_ev.height

sp_b = all_ev['n_sp_back'].to_numpy() > 0
pra_b = all_ev['n_pra_back'].to_numpy() > 0
both = (sp_b & pra_b).sum()
only_sp = (sp_b & ~pra_b).sum()
only_pra = (~sp_b & pra_b).sum()
neither = (~sp_b & ~pra_b).sum()

print("=== BACKWARD-TRACK IDENTIFICATION (theta>90, PRA non-beam) ===")
print(f"events considered (union): {N}")
print(f"  both find backward      : {both:6d}")
print(f"  only Spyral             : {only_sp:6d}")
print(f"  only ATTPCROOT-PRA      : {only_pra:6d}")
print(f"  neither                 : {neither:6d}")
print(f"Spyral backward events    : {sp_b.sum()}")
print(f"PRA    backward events    : {pra_b.sum()}  (+ {all_ev['n_pra_beam'].sum()} beam-like tracks theta>{PRA_BEAM_MAXTH})")

# --- matched single-track angle agreement ---
sp1 = sp.group_by('orig_event').agg(pl.len().alias('nsp'), pl.first('th').alias('sp_th')).rename({'orig_event':'event'})
pra1 = pra_sub.group_by('event').agg(pl.len().alias('npra'), pl.first('geo_theta').alias('pra_th'))
j = sp1.join(pra1, on='event', how='inner').filter((pl.col('nsp')==1)&(pl.col('npra')==1)).drop_nulls()
a = j['sp_th'].to_numpy(); b = j['pra_th'].to_numpy()
m = np.isfinite(a)&np.isfinite(b); a,b=a[m],b[m]

# ============ PLOTS ============
fig, ax = plt.subplots(2, 2, figsize=(13, 10))

# (1) multiplicity
ax[0,0].hist(all_ev['n_sp'].to_numpy(), bins=np.arange(-0.5,8), alpha=0.6, label=f'Spyral (mean {all_ev["n_sp"].mean():.2f})')
ax[0,0].hist(all_ev['n_pra'].to_numpy(), bins=np.arange(-0.5,8), alpha=0.6, label=f'PRA nhits>={PRA_MIN_HITS} (mean {all_ev["n_pra"].mean():.2f})')
ax[0,0].set_xlabel('clusters / tracks per event'); ax[0,0].set_ylabel('events'); ax[0,0].legend(); ax[0,0].set_title('Per-event multiplicity')

# (2) theta distributions
ax[0,1].hist(sp['th'].to_numpy(), bins=np.arange(0,181,5), alpha=0.6, label='Spyral polar')
ax[0,1].hist(pra_sub['geo_theta'].to_numpy(), bins=np.arange(0,181,5), alpha=0.6, label=f'PRA GeoTheta (nhits>={PRA_MIN_HITS})')
ax[0,1].axvline(90, color='k', ls='--'); ax[0,1].axvline(PRA_BEAM_MAXTH, color='r', ls=':', label='beam cut')
ax[0,1].set_xlabel('theta (deg)'); ax[0,1].set_ylabel('tracks'); ax[0,1].legend(); ax[0,1].set_title('Polar-angle spectra')

# (3) matched angle scatter
ax[1,0].hist2d(a, b, bins=60, range=[[0,180],[0,180]], cmap='viridis', cmin=1)
ax[1,0].plot([0,180],[0,180],'r-',lw=1)
ax[1,0].set_xlabel('Spyral polar (deg)'); ax[1,0].set_ylabel('PRA GeoTheta (deg)')
ax[1,0].set_title(f'Single-track events (n={len(a)}, med|d|={np.median(np.abs(a-b)):.1f} deg)')

# (4) backward confusion bar
cats=['both','only\nSpyral','only\nPRA','neither']; vals=[both,only_sp,only_pra,neither]
ax[1,1].bar(cats, vals, color=['green','tab:blue','tab:orange','gray'])
for i,v in enumerate(vals): ax[1,1].text(i, v, str(v), ha='center', va='bottom')
ax[1,1].set_ylabel('events'); ax[1,1].set_title('Backward-track identification (theta>90)')

plt.tight_layout()
plt.savefig(f"{OUT}/d2_clustering_compare.png", dpi=110)
print(f"\nsaved {OUT}/d2_clustering_compare.png")
PY = None
