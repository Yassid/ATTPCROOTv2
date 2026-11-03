# Compare the two codes' SIMPLE estimators (AtSpyralPID vs Spyral's Estimation phase) on the same
# events.  The earlier +0.45 deg was fitted theta (GENFIT) vs fitted theta (InterpSolver).  If the
# simple estimators agree, the offset is produced by the FITTERS; if they disagree by the same
# amount, it is upstream -- in the hit clouds or the geometry.
import glob, numpy as np, polars as pl
fs = sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/Estimation/*.parquet'))
print(f"Spyral Estimation files: {len(fs)}")
cols = pl.read_parquet(fs[0]).columns
print("columns:", cols)
run_of = lambda p: int(p.split('run_')[1].split('.')[0])
sp = pl.concat([pl.read_parquet(f).with_columns(pl.lit(run_of(f)).alias('run')) for f in fs],
               how='diagonal')
ev = 'orig_event' if 'orig_event' in cols else 'event'
sp = sp.select(['run', ev, 'polar', 'vertex_z', 'brho']).rename({ev:'event'})
# one track per event, the same convention as before: keep the largest brho
sp = sp.sort('brho', descending=True).unique(subset=['run','event'], keep='first')
sp = sp.with_columns((pl.col('polar')*180/np.pi).alias('th_sp'))
print(f"Spyral simple-estimator tracks: {len(sp)}")

at = pl.read_csv(f'{"/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad"}/attpc_pidpolar.csv')
at = at.sort('brho', descending=True).unique(subset=['run','event'], keep='first')
at = at.with_columns((pl.col('polar').abs()*180/np.pi).alias('th_at'))
at = at.with_columns(pl.when(pl.col('th_at')>90).then(180-pl.col('th_at')).otherwise(pl.col('th_at')).alias('th_at'))
print(f"ATTPCROOT simple-estimator tracks: {len(at)}")

m = at.join(sp.select(['run','event','th_sp','vertex_z']), on=['run','event'], how='inner')
print(f"MATCHED: {len(m)}")
d = (m['th_at']-m['th_sp']).to_numpy()
th = m['th_sp'].to_numpy()
ok = np.isfinite(d) & (np.abs(d)<10) & (th>20) & (th<100)
print(f"usable: {ok.sum()}")
if ok.sum()>500:
    print(f"\n  dTheta (simple estimators, ATTPCROOT - Spyral) = {d[ok].mean():+.4f} +- {d[ok].std()/np.sqrt(ok.sum()):.4f} deg   rms {d[ok].std():.3f}")
    print(f"  For comparison, the FITTED thetas differ by +0.45 deg")
    print(f"\n  {'theta':>10} {'n':>8} {'dTheta':>9}")
    for lo in range(30,85,10):
        s = ok & (th>=lo) & (th<lo+10)
        if s.sum()>100: print(f"  {lo:4d}-{lo+10:<5d} {s.sum():8d} {d[s].mean():+9.3f}")
