#!/usr/bin/env python
"""Overlay v2: fix the composition gap found by verify_overlay.py.
Three fixes over build_overlay.py:
  (1) SYNTHETIC BEAM TRACK  - dense axial line entrance->vertex (README TODO #1); the
      beam-only canvases are nearly empty (median 20 hits) so they don't supply a beam.
  (2) REDUCED CLIPPING      - sample vertex z lower so backward protons stay in-detector
      (was losing 39% of sim signal hits to z-clipping).
  (3) CHARGE CONSISTENCY    - sim signal + synthetic beam charge mapped onto the real
      charge scale (quantile match) so per-event normalization isn't a domain giveaway.
Real sparse canvas is still added (untruncated) as electronic-noise haze.
Labels: 0=proton, 1=17C, 2=beam/noise.  Output overlay_v2.parquet.
"""
import polars as pl, numpy as np
SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
SIM_PARQUET = f"{SC}/data/train.parquet"
REAL_CSV    = f"{SC}/data/real_events.csv"
OUT         = f"{SC}/data/overlay_v2.parquet"
# vertex placement: keep backward protons (extend to HIGH z) inside detector
VERTEX_Z_MIN, VERTEX_Z_MAX = 250.0, 600.0
DET_Z_MIN, DET_Z_MAX = 154.0, 1035.0
BEAM_ENTRANCE = 1000.0      # beam enters at high z
BEAM_SPACING  = 6.0         # mm between synthetic beam hits
BEAM_XY_SIGMA = 3.0         # mm beam-spot transverse spread
BEAM_Q_QLO, BEAM_Q_QHI = 0.55, 0.95   # beam is high-dE/dx -> upper real-charge quantiles
SEED = 12345
rng = np.random.default_rng(SEED)

real = pl.read_csv(REAL_CSV)
real_q_sorted = np.sort(real.filter(pl.col('q') > 0)['q'].to_numpy())
def real_q_at(quant):                       # real charge at given quantile(s)
    return np.interp(np.clip(quant, 0, 1), np.linspace(0, 1, len(real_q_sorted)), real_q_sorted)
sim_all = pl.read_parquet(SIM_PARQUET)
sim_q_sorted = np.sort(sim_all.filter(pl.col('charge') > 0)['charge'].to_numpy())
def qmap(q):                                # sim charge -> real charge (histogram match)
    quant = np.searchsorted(sim_q_sorted, q, side='right') / len(sim_q_sorted)
    return real_q_at(quant)

# sparse real beam-only canvases for electronic-noise haze (NOT reaction events, whose
# real tracks would be mislabeled as noise). Same density selector as build_overlay.py.
from scipy.spatial import cKDTree
rr = real.with_columns((pl.col('x')**2 + pl.col('y')**2).sqrt().alias('r'))
canvases = []
for _, g in rr.group_by('event'):
    off = g.filter(pl.col('r') >= 35.0)
    if off.height >= 5:
        Q = off.select(['x', 'y', 'z']).to_numpy()
        if cKDTree(Q).query_ball_point(Q, r=12.0, return_length=True).max() >= 7:
            continue                                   # reaction-like -> not a noise canvas
    canvases.append(g.select(['x', 'y', 'z', 'q']).to_numpy())
print(f"noise canvases (beam-only): {len(canvases)}  (median {np.median([c.shape[0] for c in canvases]):.0f} hits)")

sim = sim_all.with_columns((pl.col('px')**2 + pl.col('py')**2 + pl.col('pz')**2).sqrt().alias('pmag'))
rows = []; out_ev = 0; clip_lost = 0; sim_tot = 0
for ev, g in sim.group_by('event'):
    prot = g.filter(pl.col('particle') == 'proton')
    if prot.height == 0:
        continue
    v = prot.sort('pmag', descending=True).row(0, named=True)
    tz = float(rng.uniform(VERTEX_Z_MIN, VERTEX_Z_MAX))
    dx, dy, dz = -v['x'], -v['y'], tz - v['z']
    for h in g.iter_rows(named=True):
        sim_tot += 1
        zz = h['z'] + dz
        if zz < DET_Z_MIN or zz > DET_Z_MAX:
            clip_lost += 1; continue
        lbl = 0 if h['particle'] == 'proton' else 1
        rows.append((out_ev, h['x'] + dx, h['y'] + dy, zz, float(qmap(h['charge'])), lbl, 'sim'))
    # (1) synthetic beam: dense axial line entrance -> vertex
    nb = max(1, int((BEAM_ENTRANCE - tz) / BEAM_SPACING))
    bz = np.linspace(tz, BEAM_ENTRANCE, nb)
    bx = rng.normal(0, BEAM_XY_SIGMA, nb); by = rng.normal(0, BEAM_XY_SIGMA, nb)
    bq = real_q_at(rng.uniform(BEAM_Q_QLO, BEAM_Q_QHI, nb))
    for x, y, z, q in zip(bx, by, bz, bq):
        rows.append((out_ev, float(x), float(y), float(z), float(q), 2, 'beam'))
    # (3 cont.) sparse real-noise haze
    C = canvases[rng.integers(len(canvases))]
    for x, y, z, q in C:
        rows.append((out_ev, float(x), float(y), float(z), float(q), 2, 'noise'))
    out_ev += 1

df = pl.DataFrame(rows, schema=['event', 'x', 'y', 'z', 'charge', 'label', 'source'], orient='row')
df.write_parquet(OUT)
print(f"sim clip loss: {100*clip_lost/max(1,sim_tot):.0f}%  (was 39%)")
print(f"wrote {OUT}: {df.height} hits, {df['event'].n_unique()} events")
print("label mix:", dict(zip(*[s.to_list() for s in df.group_by('label').len().sort('label')])))
print(f"hits/event: median {df.group_by('event').len()['len'].median():.0f}")
