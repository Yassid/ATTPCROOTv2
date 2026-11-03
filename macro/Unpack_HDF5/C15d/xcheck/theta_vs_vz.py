# Re-derive dTheta(vz) with the unit trap handled: Spyral's PARQUET stores the vertex in METRES.
import glob, numpy as np, polars as pl
fs=sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet'))
sp=pl.concat([pl.read_parquet(f,columns=['orig_run','orig_event','ke','polar','azimuthal','vertex_z','redchisq']) for f in fs])
sp=sp.filter((pl.col('redchisq')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['orig_run','orig_event'],keep='first').rename({'orig_run':'run','orig_event':'event'})\
   .with_columns((pl.col('polar')*180/np.pi).alias('th_s'),(pl.col('azimuthal')*180/np.pi).alias('ph_s'),
                 (pl.col('vertex_z')*1000.0).alias('vz_mm'))     # metres -> mm
at=pl.read_csv('/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_phi.csv')
at=at.filter((pl.col('chi2ndf')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['run','event'],keep='first').rename({'theta':'th_a','phi':'ph_a','vz':'vz_a'})
m=at.join(sp.select(['run','event','th_s','ph_s','vz_mm']),on=['run','event'],how='inner')
d=(m['th_a']-m['th_s']).to_numpy()
dp=((m['ph_a'].to_numpy()-m['ph_s'].to_numpy())%360+180)%360-180
vzs=m['vz_mm'].to_numpy(); vza=m['vz_a'].to_numpy()
core=(np.abs(d)<10)&(np.abs(dp-9.55)<3)
print(f"core {core.sum()} tracks;  Spyral vz range {vzs[core].min():.0f}-{vzs[core].max():.0f} mm,"
      f"  ATTPCROOT vz range {vza[core].min():.0f}-{vza[core].max():.0f} mm")
print(f"\n{'vz [mm]':>12} {'n':>7} {'dTheta':>9} {'err':>7} {'dVz':>8}")
cen,val,err=[],[],[]
for lo in range(0,1000,125):
    s=core&(vzs>=lo)&(vzs<lo+125)
    if s.sum()<100: continue
    mu=d[s].mean(); e=d[s].std()/np.sqrt(s.sum())
    print(f"{lo:5d}-{lo+125:<6d} {s.sum():7d} {mu:+9.3f} {e:7.3f} {(vza[s]-vzs[s]).mean():+8.1f}")
    cen.append(lo+62.5); val.append(mu); err.append(e)
cen,val,err=map(np.array,(cen,val,err))
if len(cen)>2:
    W=1/err**2
    A=np.vstack([np.ones_like(cen),cen]).T
    c=np.linalg.lstsq(A*np.sqrt(W)[:,None],val*np.sqrt(W),rcond=None)[0]
    cov=np.linalg.inv((A*W[:,None]).T@A)
    print(f"\n  slope = {c[1]*1000:+.3f} +- {np.sqrt(cov[1,1])*1000:.3f} deg per 1000 mm")
    print(f"  Garfield straight-line estimate predicts about +0.27 deg per 1000 mm (ATTPCROOT minus Spyral, RISING)")
