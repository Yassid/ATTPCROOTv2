#!/usr/bin/env python
# Build event-matched CSVs for the ATTPCROOT-vs-Spyral viewer.
#   ATTPCROOT entry index -> Spyral global event # via an empirically-found offset (chunk start).
#   ATTPCROOT z (CalculateZGeo) is FLIPPED into Spyral's frame so tracks overlay.
# Out: data/cmp_attpc_<run>.csv, data/cmp_spyral_<run>.csv  (event=global #, x,y,z,q,cluster)
# Then view with:  root -l 'viewer/view_events.C("data/cmp_attpc_0305.csv","data/cmp_spyral_0305.csv")'
import sys
import pandas as pd, numpy as np

RUN = int(sys.argv[1]) if len(sys.argv) > 1 else 305
# ATTPCROOT entry i = GET event (OFFSET + i)  -- GET events are consecutive from the h5's first event.
# OFFSET = first GET event number of the chunk (run_0305 -> 34675). Pass as arg 2 for other runs.
OFFSET = int(sys.argv[2]) if len(sys.argv) > 2 else 34675
ZFLIP = 1137.0  # z_spyral ~ ZFLIP - z_attpc (CalculateZGeo vs Spyral have opposite drift orientation)

at = pd.read_csv(f"data/cloud_run_0{RUN}.csv")          # event=entry idx, cluster=triplclust label
sp = pd.read_parquet("data/spyral_all.parquet"); sp = sp[sp.run == RUN]

at = at.copy()
best_off = OFFSET
at['glob'] = best_off + at['event']
# verify alignment on CLUSTERED y-range (frame-independent, noise-free)
atc = at[at.cluster >= 0]
at_yr = atc.groupby('glob').y.agg(lambda s: s.max() - s.min())
sp_yr = sp.groupby('event').y.agg(lambda s: s.max() - s.min())
at['z'] = ZFLIP - at['z']                                # flip into Spyral frame
common = sorted(set(at['glob']) & set(sp.event))
ato = at[at['glob'].isin(common)][['glob', 'x', 'y', 'z', 'q', 'cluster']].rename(columns={'glob': 'event'})
spo = sp[sp.event.isin(common)][['event', 'x', 'y', 'z', 'q', 'label']].rename(columns={'label': 'cluster'})
ato.to_csv(f"data/cmp_attpc_0{RUN}.csv", index=False)
spo.to_csv(f"data/cmp_spyral_0{RUN}.csv", index=False)

# match-quality report on the matched set (clustered y-range)
dy = [abs(at_yr[g] - sp_yr[g]) for g in common if g in at_yr.index and g in sp_yr.index]
print(f"offset {best_off}: {len(common)} matched events -> data/cmp_attpc_0{RUN}.csv / data/cmp_spyral_0{RUN}.csv")
print(f"median |clustered y-range diff|: {np.median(dy):.1f} mm  (small = correct matching)")
