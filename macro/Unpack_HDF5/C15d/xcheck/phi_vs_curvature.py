# If the 9.55 deg azimuth offset is the two codes quoting the direction at different points along
# the trajectory, it must scale with CURVATURE: dphi = (arc length offset) / radius, and radius
# grows with momentum.  So dphi should fall as KE rises, and show no theta dependence beyond
# what p_T implies.  A frame/convention rotation would instead be constant in everything.
import glob, numpy as np, polars as pl
fs=sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet'))
sp=pl.concat([pl.read_parquet(f,columns=['orig_run','orig_event','ke','brho','polar','azimuthal','redchisq']) for f in fs])
sp=sp.filter((pl.col('redchisq')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['orig_run','orig_event'],keep='first').rename({'orig_run':'run','orig_event':'event'})\
   .with_columns((pl.col('polar')*180/np.pi).alias('th_s'),(pl.col('azimuthal')*180/np.pi).alias('ph_s'),
                 pl.col('ke').alias('ke_s'))
at=pl.read_csv('/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_phi.csv')
at=at.filter((pl.col('chi2ndf')<5)&(pl.col('ke')>6)).sort('ke',descending=True)\
   .unique(subset=['run','event'],keep='first').rename({'ke':'ke_a','theta':'th_a','phi':'ph_a'})
m=at.join(sp.select(['run','event','ke_s','brho','th_s','ph_s']),on=['run','event'],how='inner')
d=(m['th_a']-m['th_s']).to_numpy()
dp=((m['ph_a'].to_numpy()-m['ph_s'].to_numpy())%360+180)%360-180
ke=m['ke_s'].to_numpy(); br=m['brho'].to_numpy(); th=m['th_s'].to_numpy()
ok=(np.abs(d)<10)&(np.abs(dp)<40)
print("=== dphi vs KE  (curvature test) ===")
print(f"{'KE [MeV]':>12} {'n':>7} {'<dphi>':>9} {'err':>7} {'<Brho>':>8} {'pT-pred':>9}")
# pT = Brho * q * c  ~  Brho[T m]*300 MeV/c per unit charge; radius_T = Brho/B ; dphi ~ s/R
rows=[]
for lo,hi in ((6,12),(12,18),(18,24),(24,30),(30,36),(36,45),(45,60)):
    s=ok&(ke>=lo)&(ke<hi)
    if s.sum()<300: continue
    pt=br[s].mean()*np.sin(np.deg2rad(th[s].mean()))
    print(f"{lo:5d}-{hi:<6d} {s.sum():7d} {dp[s].mean():+9.3f} {dp[s].std()/np.sqrt(s.sum()):7.3f} {br[s].mean():8.4f} {1/max(pt,1e-9):9.3f}")
    rows.append((dp[s].mean(), br[s].mean(), pt))
r=np.array(rows)
if len(r)>2:
    cc=np.corrcoef(r[:,0], 1/r[:,2])[0,1]
    print(f"\n  correlation( <dphi> , 1/pT ) over the KE bins = {cc:+.3f}")
    print("  (near +1 => the offset is an arc-length/reference-point effect;")
    print("   near  0 => a rigid frame rotation, independent of momentum)")
print("\n=== and dtheta over the same KE bins (should NOT follow curvature) ===")
for lo,hi in ((6,12),(12,18),(18,24),(24,30),(30,36),(36,45),(45,60)):
    s=ok&(ke>=lo)&(ke<hi)
    if s.sum()<300: continue
    print(f"  KE {lo:3d}-{hi:<3d}: n={s.sum():6d}  <dtheta>={d[s].mean():+6.3f} +- {d[s].std()/np.sqrt(s.sum()):.3f}")
