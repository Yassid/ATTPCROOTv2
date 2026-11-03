# Is the ATTPCROOT-Spyral theta offset a DETECTOR TILT or a flat convention/scale offset?
# A tilt of the assumed beam axis by angle B toward azimuth phi0 makes the apparent polar angle
# differ by  dtheta = B*cos(phi - phi0)  -- a sinusoid in phi.  A drift-velocity or vertex-scale
# error instead shifts every track the same way, giving a FLAT dtheta(phi).
import glob, numpy as np, polars as pl

# ---- Spyral side
fs = sorted(glob.glob('/home/yassid/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet'))
sp = pl.concat([pl.read_parquet(f, columns=['orig_run','orig_event','ke','polar','azimuthal',
                                            'vertex_x','vertex_y','vertex_z','redchisq'])
                for f in fs])
sp = sp.filter((pl.col('redchisq') < 5) & (pl.col('ke') > 6))
# keep the most energetic track per (run,event), as every earlier comparison did
sp = sp.sort('ke', descending=True).unique(subset=['orig_run','orig_event'], keep='first')
sp = sp.rename({'orig_run':'run','orig_event':'event'}).with_columns(
        (pl.col('polar')*180/np.pi).alias('th_s'),
        (pl.col('azimuthal')*180/np.pi).alias('ph_s'),
        pl.col('ke').alias('ke_s'), pl.col('vertex_z').alias('vz_s'))
print(f"Spyral tracks: {len(sp)}")

# ---- ATTPCROOT side
at = pl.read_csv(f'{S}/attpc_phi.csv'.replace('$S','')) if False else pl.read_csv(
        '/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_phi.csv')
at = at.filter((pl.col('chi2ndf') < 5) & (pl.col('ke') > 6))
at = at.sort('ke', descending=True).unique(subset=['run','event'], keep='first')
at = at.rename({'ke':'ke_a','theta':'th_a','phi':'ph_a','vz':'vz_a'})
print(f"ATTPCROOT tracks: {len(at)}")

m = at.join(sp.select(['run','event','ke_s','th_s','ph_s','vz_s']), on=['run','event'], how='inner')
print(f"MATCHED on (run,event): {len(m)}")

dth = (m['th_a'] - m['th_s']).to_numpy()
pa  = m['ph_a'].to_numpy() % 360.0
ps  = m['ph_s'].to_numpy() % 360.0
ok  = np.abs(dth) < 10                      # drop mismatched-track outliers
print(f"  |dtheta|<10 deg: {ok.sum()}  ({100*ok.mean():.1f} %)   mean dtheta {dth[ok].mean():+.3f}  rms {dth[ok].std():.3f}")

# ---- which azimuth convention do the two codes share?
def circ(d):
    d = (d + 180) % 360 - 180
    r = np.abs(np.mean(np.exp(1j*np.deg2rad(d))))
    return d.mean(), d.std(), r
for name, d in (('phiA - phiS', pa-ps), ('phiA + phiS', pa+ps), ('phiA - phiS + 180', pa-ps+180)):
    mu, sd, r = circ(d[ok])
    print(f"  {name:20s}  mean {mu:+8.2f}  rms {sd:6.2f}  concentration R={r:.3f}")

# ---- THE TEST: dtheta vs azimuth.  Tilt -> sinusoid; scale/convention -> flat line.
print("\n=== dtheta (ATTPCROOT - Spyral) vs azimuth, 12 bins of 30 deg ===")
print(f"{'phi_S [deg]':>14} {'n':>7} {'dtheta':>9} {'err':>7} {'<dphi>':>8}")
cen, val, err = [], [], []
dph = (pa - ps + 180) % 360 - 180
for b in range(12):
    lo, hi = b*30, b*30+30
    s = ok & (ps >= lo) & (ps < hi)
    if s.sum() < 100: continue
    mu = dth[s].mean(); e = dth[s].std()/np.sqrt(s.sum())
    print(f"{lo:6.0f}-{hi:<7.0f} {s.sum():7d} {mu:+9.3f} {e:7.3f} {dph[s].mean():+8.2f}")
    cen.append(lo+15); val.append(mu); err.append(e)
cen, val, err = map(np.array, (cen, val, err))

# least squares for  A + B*cos(phi - phi0)   ==  A + C*cos(phi) + Sx*sin(phi)
X = np.column_stack([np.ones_like(cen), np.cos(np.deg2rad(cen)), np.sin(np.deg2rad(cen))])
W = 1/err**2
coef, *_ = np.linalg.lstsq(X*np.sqrt(W)[:,None], val*np.sqrt(W), rcond=None)
A, C, Sx = coef
B = np.hypot(C, Sx); phi0 = np.degrees(np.arctan2(Sx, C)) % 360
cov = np.linalg.inv((X*W[:,None]).T @ X)
sA, sC, sS = np.sqrt(np.diag(cov))
sB = np.hypot(C*sC, Sx*sS)/max(B, 1e-9)
chi2_fit  = np.sum(W*(val - X@coef)**2)
chi2_flat = np.sum(W*(val - np.average(val, weights=W))**2)
print(f"\n  flat offset      A    = {A:+.3f} +- {sA:.3f} deg")
print(f"  sinusoid ampl.   B    = {B:.3f} +- {sB:.3f} deg   toward phi0 = {phi0:.0f} deg")
print(f"  chi2/ndf  flat only : {chi2_flat/(len(cen)-1):8.2f}")
print(f"  chi2/ndf  A+Bcos    : {chi2_fit/(len(cen)-3):8.2f}")
print(f"  -> the modulation is {B/max(abs(A),1e-9)*100:.0f} % of the mean offset")
