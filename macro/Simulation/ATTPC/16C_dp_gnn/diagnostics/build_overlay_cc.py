#!/usr/bin/env python
"""Charge-CONSISTENT overlay: same vertex-conditioned beam+noise compositing as
build_overlay.py, but the sim signal charge is QUANTILE-MAPPED onto the real charge
distribution first, so sim-signal and real-beam/noise hits share one charge scale.
Without this, a per-event charge normalization would trivially separate sim(label 0/1)
from real(label 2) hits — a domain giveaway. Output overlay_cc.parquet.
Run: ~/Spyral/venv/bin/python build_overlay_cc.py
"""
import polars as pl, numpy as np
from scipy.spatial import cKDTree

SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
SIM_PARQUET = f"{SC}/data/train.parquet"
REAL_CSV    = f"{SC}/data/real_events.csv"
OUT         = f"{SC}/data/overlay_cc.parquet"
BEAM_ENTERS_HIGH_Z = True
VERTEX_Z_MIN, VERTEX_Z_MAX = 300.0, 850.0
R_OFFAXIS   = 35.0
DENS_R, DENS_THRESH = 12.0, 7
BEAM_AXIS_XY = (0.0, 0.0)
DET_Z_MIN, DET_Z_MAX = 154.0, 1035.0
SEED = 12345
rng = np.random.default_rng(SEED)

# ---- charge quantile-map: sim signal charge -> real charge distribution ----
real = pl.read_csv(REAL_CSV)
real_q = real.filter(pl.col('q') > 0)['q'].to_numpy()
real_q_sorted = np.sort(real_q)
sim_all = pl.read_parquet(SIM_PARQUET)
sim_q = sim_all.filter(pl.col('charge') > 0)['charge'].to_numpy()
sim_q_sorted = np.sort(sim_q)
def qmap(q):
    """map a sim charge to the real charge at the same quantile (histogram matching)."""
    quant = np.searchsorted(sim_q_sorted, q, side='right') / len(sim_q_sorted)
    quant = np.clip(quant, 0, 1)
    return np.interp(quant, np.linspace(0, 1, len(real_q_sorted)), real_q_sorted)
print(f"charge qmap: sim med {np.median(sim_q):.1f} -> real med {np.median(real_q):.1f}"
      f"  (check: qmap(simMed)={qmap(np.median(sim_q)):.1f})")

# ---- 1. real beam-only canvases ----
real = real.with_columns((pl.col('x')**2 + pl.col('y')**2).sqrt().alias('r'))
canvases = {}
for ev, g in real.group_by('event'):
    off = g.filter(pl.col('r') >= R_OFFAXIS)
    if off.height >= 5:
        Q = off.select(['x', 'y', 'z']).to_numpy()
        dens = cKDTree(Q).query_ball_point(Q, r=DENS_R, return_length=True)
        if dens.max() >= DENS_THRESH:
            continue
    canvases[int(ev[0])] = g.select(['x', 'y', 'z', 'q']).to_numpy()
canvas_ids = list(canvases)
print(f"beam-only canvases: {len(canvas_ids)} / {real['event'].n_unique()} real events")

# ---- 2. overlay sim products (charge-mapped) onto truncated canvases ----
sim = sim_all.with_columns((pl.col('px')**2 + pl.col('py')**2 + pl.col('pz')**2).sqrt().alias('pmag'))
rows = []; out_ev = 0
for ev, g in sim.group_by('event'):
    prot = g.filter(pl.col('particle') == 'proton')
    if prot.height == 0:
        continue
    v = prot.sort('pmag', descending=True).row(0, named=True)
    vx, vy, vz = v['x'], v['y'], v['z']
    tz = float(rng.uniform(VERTEX_Z_MIN, VERTEX_Z_MAX))
    dx, dy, dz = BEAM_AXIS_XY[0] - vx, BEAM_AXIS_XY[1] - vy, tz - vz
    for h in g.iter_rows(named=True):
        zz = h['z'] + dz
        if zz < DET_Z_MIN or zz > DET_Z_MAX:
            continue
        lbl = 0 if h['particle'] == 'proton' else 1
        rows.append((out_ev, h['x'] + dx, h['y'] + dy, zz, float(qmap(h['charge'])), lbl, 'sim'))
    cid = canvas_ids[rng.integers(len(canvas_ids))]
    C = canvases[cid]
    keep = C[:, 2] > tz if BEAM_ENTERS_HIGH_Z else C[:, 2] < tz
    for x, y, z, q in C[keep]:
        rows.append((out_ev, x, y, z, q, 2, 'data'))
    out_ev += 1

df = pl.DataFrame(rows, schema=['event', 'x', 'y', 'z', 'charge', 'label', 'source'], orient='row')
df.write_parquet(OUT)
print(f"wrote {OUT}: {df.height} hits, {df['event'].n_unique()} overlaid events")
print("label mix:", dict(zip(*[s.to_list() for s in df.group_by('label').len().sort('label')])))
print(f"hits/event: median {df.group_by('event').len()['len'].median():.0f}")
