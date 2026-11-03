import glob, numpy as np, polars as pl
fs=sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet'))
sp=pl.concat([pl.read_parquet(f,columns=['orig_run','orig_event','ke','polar','azimuthal','vertex_z','redchisq']) for f in fs])
sp=sp.filter((pl.col('redchisq')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['orig_run','orig_event'],keep='first').rename({'orig_run':'run','orig_event':'event'})\
   .with_columns((pl.col('polar')*180/np.pi).alias('th_s'),(pl.col('azimuthal')*180/np.pi).alias('ph_s'),
                 pl.col('ke').alias('ke_s'),pl.col('vertex_z').alias('vz_s'))
at=pl.read_csv('/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_phi.csv')
at=at.filter((pl.col('chi2ndf')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['run','event'],keep='first').rename({'ke':'ke_a','theta':'th_a','phi':'ph_a','vz':'vz_a'})
m=at.join(sp.select(['run','event','ke_s','th_s','ph_s','vz_s']),on=['run','event'],how='inner')
d=(m['th_a']-m['th_s']).to_numpy(); ths=m['th_s'].to_numpy()
dp=((m['ph_a'].to_numpy()-m['ph_s'].to_numpy())%360+180)%360-180
kr=(m['ke_a']/m['ke_s']).to_numpy(); vz=m['vz_s'].to_numpy(); c2=m['chi2ndf'].to_numpy()
# clean core: tracks whose azimuth obeys the rigid rotation to +-3 deg
core=(np.abs(d)<10)&(np.abs(dp-9.55)<3)
print(f"core sample (|dphi-9.55|<3 deg): {core.sum()} of {(np.abs(d)<10).sum()}")
print(f"  <dtheta> = {d[core].mean():+.3f}   <KE ratio> = {kr[core].mean():.4f}\n")
print("=== dtheta vs theta  (a FIT bias should vary with theta; a convention offset should not) ===")
print(f"{'theta_S':>12} {'n':>7} {'dtheta':>9} {'err':>7} {'KEratio':>9}")
for lo in range(40,86,5):
    s=core&(ths>=lo)&(ths<lo+5)
    if s.sum()<200: continue
    print(f"{lo:5d}-{lo+5:<6d} {s.sum():7d} {d[s].mean():+9.3f} {d[s].std()/np.sqrt(s.sum()):7.3f} {kr[s].mean():9.4f}")
print("\n=== dtheta vs chi2/ndf of the ATTPCROOT fit (poorly sampled tracks fit worse) ===")
for lo,hi in ((0,0.2),(0.2,0.5),(0.5,1),(1,2),(2,3),(3,5)):
    s=core&(c2>=lo)&(c2<hi)
    if s.sum()<200: continue
    print(f"  chi2/ndf {lo:4.1f}-{hi:<4.1f}: n={s.sum():6d}  <dtheta>={d[s].mean():+6.3f} +- {d[s].std()/np.sqrt(s.sum()):.3f}   <KEratio>={kr[s].mean():.4f}")
print("\n=== dtheta vs vertex z ===")
for lo in range(0,1000,125):
    s=core&(vz>=lo)&(vz<lo+125)
    if s.sum()<200: continue
    print(f"  vz {lo:4d}-{lo+125:<4d}: n={s.sum():6d}  <dtheta>={d[s].mean():+6.3f}   <KEratio>={kr[s].mean():.4f}")
