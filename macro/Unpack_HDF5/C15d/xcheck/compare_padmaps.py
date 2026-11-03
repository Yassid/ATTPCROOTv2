# Compare the two codes' pad -> (x,y) tables.  If they differ by a rigid rotation, the per-event
# azimuth offset is explained by geometry rather than by anything physical.
import numpy as np, csv
sp=np.genfromtxt('/home/yassid/Data/spyral/venv_1.1.2/lib/python3.11/site-packages/spyral/data/padxy.csv',
                 delimiter=',',names=True,usecols=(0,1))
sx,sy=sp['x'],sp['y']
at=np.genfromtxt(f'{"/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad"}/attpc_padxy.csv',
                 delimiter=',',names=True)
ax=np.full(len(sx),np.nan); ay=np.full(len(sx),np.nan)
pid=at['pad'].astype(int); ok=pid<len(sx)
ax[pid[ok]]=at['x'][ok]; ay[pid[ok]]=at['y'][ok]
g=np.isfinite(ax)&np.isfinite(sx)&((sx**2+sy**2)>400)   # skip the centre, where angle is ill-defined
print(f"pads compared: {g.sum()}")
ra,rs=np.hypot(ax[g],ay[g]),np.hypot(sx[g],sy[g])
print(f"radius:  ATTPCROOT mean {ra.mean():.2f} mm   Spyral mean {rs.mean():.2f} mm   ratio {ra.mean()/rs.mean():.5f}")
print(f"         per-pad |r_A - r_S|  mean {np.abs(ra-rs).mean():.3f} mm  max {np.abs(ra-rs).max():.3f} mm")
da=(np.degrees(np.arctan2(ay[g],ax[g])-np.arctan2(sy[g],sx[g]))+180)%360-180
R=np.abs(np.mean(np.exp(1j*np.deg2rad(da))))
print(f"\nper-pad azimuth difference: mean {da.mean():+.4f} deg  median {np.median(da):+.4f}  rms {da.std():.4f}  R={R:.5f}")
# best rigid rotation by least squares (Procrustes on the 2-D point sets)
num=np.sum(ax[g]*sy[g]-ay[g]*sx[g]); den=np.sum(ax[g]*sx[g]+ay[g]*sy[g])
rot=np.degrees(np.arctan2(num,den))
print(f"best-fit rigid rotation (Spyral -> ATTPCROOT): {-rot:+.4f} deg")
c,s=np.cos(np.deg2rad(-rot)),np.sin(np.deg2rad(-rot))
rx,ry=c*sx[g]-s*sy[g], s*sx[g]+c*sy[g]
res=np.hypot(ax[g]-rx,ay[g]-ry)
print(f"residual after that rotation: mean {res.mean():.3f} mm  median {np.median(res):.3f}  max {res.max():.3f}")
res0=np.hypot(ax[g]-sx[g],ay[g]-sy[g])
print(f"residual with NO rotation:    mean {res0.mean():.3f} mm  median {np.median(res0):.3f}  max {res0.max():.3f}")
