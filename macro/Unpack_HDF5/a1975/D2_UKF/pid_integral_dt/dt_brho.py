import numpy as np
d=np.load('/mnt/f/a1975/sim_rebuild/trig/dt/joined_dt.npz'); T=d['T']; s=d['sel'].astype(bool)
print(T[:2]); ke=T[:,3]; m=2808.921; br=np.sqrt(ke**2+2*m*ke)/299.792458
for lab,w in [('adopted (sel)',s),('adopted, Ex<5 (bound)',s&(d['ex']<5))]:
    b=br[w]; print(lab,w.sum(),'Brho pct 5/50/95',np.percentile(b,[5,50,95]).round(3))
    for lo,hi in [(0,0.3),(0.3,0.6),(0.6,0.9),(0.9,1.2),(1.2,1.5),(1.5,9)]: print(f'   {lo}-{hi}: {np.mean((b>=lo)&(b<hi))*100:5.1f}%')
