# The dtheta(phi) structure is not a 1-fold sinusoid (no tilt).  The AT-TPC pad plane is a
# hexagonal arrangement of triangular pads, so test harmonics: n=6 would be pad geometry,
# n=1 a tilt, n=2 a quadrupole field/scale effect.
import glob, numpy as np, polars as pl
fs = sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet'))
sp = pl.concat([pl.read_parquet(f, columns=['orig_run','orig_event','ke','polar','azimuthal','vertex_z','redchisq']) for f in fs])
sp = sp.filter((pl.col('redchisq')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
       .unique(subset=['orig_run','orig_event'],keep='first')\
       .rename({'orig_run':'run','orig_event':'event'})\
       .with_columns((pl.col('polar')*180/np.pi).alias('th_s'),
                     (pl.col('azimuthal')*180/np.pi).alias('ph_s'),
                     pl.col('ke').alias('ke_s'), pl.col('vertex_z').alias('vz_s'))
at = pl.read_csv(f'{"/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad"}/attpc_phi.csv')
at = at.filter((pl.col('chi2ndf')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
       .unique(subset=['run','event'],keep='first')\
       .rename({'ke':'ke_a','theta':'th_a','phi':'ph_a','vz':'vz_a'})
m = at.join(sp.select(['run','event','ke_s','th_s','ph_s','vz_s']),on=['run','event'],how='inner')
dth=(m['th_a']-m['th_s']).to_numpy(); ps=m['ph_s'].to_numpy()%360
pa=m['ph_a'].to_numpy()%360; ker=(m['ke_a']/m['ke_s']).to_numpy()
ok=np.abs(dth)<10
d,p,kr,paa = dth[ok],ps[ok],ker[ok],pa[ok]

NB=36
edges=np.arange(0,361,360/NB); idx=np.digitize(p,edges)-1
cen=np.array([(edges[i]+edges[i+1])/2 for i in range(NB)])
val=np.array([d[idx==i].mean() for i in range(NB)])
err=np.array([d[idx==i].std()/np.sqrt(max((idx==i).sum(),1)) for i in range(NB)])
cnt=np.array([(idx==i).sum() for i in range(NB)])
print(f"36 bins of 10 deg, {cnt.min()}-{cnt.max()} tracks each")
print("\n=== harmonic decomposition of dtheta(phi) ===")
print(f"{'n':>3} {'period':>8} {'amplitude':>11} {'err':>7} {'phase':>7}")
W=1/err**2; mu=np.average(val,weights=W)
tot=np.sum(W*(val-mu)**2)
for n in range(1,9):
    X=np.column_stack([np.ones(NB),np.cos(n*np.deg2rad(cen)),np.sin(n*np.deg2rad(cen))])
    c,*_=np.linalg.lstsq(X*np.sqrt(W)[:,None],val*np.sqrt(W),rcond=None)
    cov=np.linalg.inv((X*W[:,None]).T@X); s=np.sqrt(np.diag(cov))
    A=np.hypot(c[1],c[2]); sA=np.hypot(c[1]*s[1],c[2]*s[2])/max(A,1e-9)
    ph=np.degrees(np.arctan2(c[2],c[1]))%(360/n)
    mark=" <<<" if A>3*sA and A>0.05 else ""
    print(f"{n:3d} {360/n:7.1f}deg {A:11.4f} {sA:7.4f} {ph:7.1f}{mark}")
print(f"\n  weighted rms of dtheta(phi) about its mean: {np.sqrt(tot/NB):.3f} deg")
print(f"  mean offset: {mu:+.3f} deg")

print("\n=== is the 9.55 deg azimuth rotation uniform in phi? ===")
dp=(paa-p+180)%360-180
for b in range(6):
    s=(p>=b*60)&(p<b*60+60)
    print(f"  phi_S {b*60:3d}-{b*60+60:3d}: <dphi> = {dp[s].mean():+6.2f} +- {dp[s].std()/np.sqrt(s.sum()):.2f}  rms {dp[s].std():5.2f}")
print("\n=== and does dtheta track the per-event dphi? (a shared mis-rotation would correlate) ===")
for lo,hi in ((-30,-10),(-10,-5),(-5,0),(0,5),(5,10),(10,15),(15,20),(20,40)):
    s=(dp>=lo)&(dp<hi)
    if s.sum()<200: continue
    print(f"  dphi {lo:+4d}..{hi:+4d}: n={s.sum():6d}  <dtheta>={d[s].mean():+6.3f}  <KEratio>={kr[s].mean():.3f}")
