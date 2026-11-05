# Does the mesh fire early because the heavy residual leaves the beam hole? r_res at the pad plane from two-body momentum balance of each track, against the time of the first mesh crossing.
import numpy as np, glob, os
D = '/mnt/f/a1975/sim_rebuild/trig/'; rows = []; names = None
for f in sorted(glob.glob(D + 'mesh/run_*.npz')):
    z = np.load(f); names = list(z['names']); rows.append(np.c_[np.full(len(z['d']), int(os.path.basename(f)[4:8])), z['d']])
A = np.vstack(rows); col = {n: i + 1 for i, n in enumerate(names)}; g = lambda n: A[:, col[n]]; print('runs', len(rows), 'events', len(A))
key = {(int(r), int(e)): i for i, (r, e) in enumerate(zip(A[:, 0], g('ev')))}; T = np.loadtxt(D + 'tracks.txt'); idx = np.array([key.get((int(r), int(e) + (1 if int(r) >= 31 else 0)), -1) for r, e in zip(T[:, 0], T[:, 1])]); T = T[idx >= 0]; idx = idx[idx >= 0]
c = dict(run=0, ev=1, gate=2, th=3, vz=4, ke=5, ex=6, tc=7, zmin=8, zmax=9, nh=10, rzmin=11, rzmax=12, dirB=13)
s = (T[:, 2] == 1) & (T[:, 13] == 1) & (T[:, 11] < 30); t = T[s]; i = idx[s]
th = np.radians(t[:, 3]); ke = t[:, 5]; vz = t[:, 4]; p = np.sqrt(ke**2 + 2 * ke * 938.272); pT = p * np.sin(th); pz = p * np.cos(th)
Eb = 183.873 - 0.011099 * vz; M = 16.0147013 * 931.494; pb = np.sqrt(Eb**2 + 2 * Eb * M); tres = np.arctan2(pT, pb - pz); L = 1000 - vz; rres = L * np.tan(tres)
m = A[i, col['m60']]; ic = A[i, col['icbin']]; rel = np.where(m >= 0, m - ic, 1e6); early = rel < 100; no = m < 0; pk = A[i, col['mpk']]
print('gated backward tracks', len(t))
print('\nfraction of events with mesh firing within 100 bins of the IC, by predicted residual radius at the pad plane')
for a, b in [(0,5),(5,10),(10,14),(14,17),(17,20),(20,23),(23,26),(26,30),(30,35),(35,45),(45,60),(60,100)]:
    w = (rres >= a) & (rres < b)
    if w.sum() > 5: print(f'  r_res {a:3d}-{b:3d}: n {w.sum():5d} early {early[w].mean():.3f}  no crossing {no[w].mean():.3f}  median first crossing {np.median(rel[w & ~no]) if (w & ~no).sum() else -1:6.0f}  mesh peak {np.median(pk[w]):6.0f}')
TB = [(10,20),(20,25),(25,30),(30,35),(35,40)]; ZB = [200,300,400,500,600,700,800,900]
for nm, ww in [('bound Ex -1.5..0.73', (t[:,6]>-1.5)&(t[:,6]<0.729)), ('Ex 3-5', (t[:,6]>3)&(t[:,6]<5))]:
    print('\n', nm, ': n and early fraction; then median r_res; rows vz, columns theta_cm 10-20 20-25 25-30 30-35 35-40')
    for z0, z1 in zip(ZB[:-1], ZB[1:]):
        l = f'  {z0}-{z1}'
        for a, b in TB:
            w = ww & (vz >= z0) & (vz < z1) & (t[:,7] >= a) & (t[:,7] < b); l += f'  {w.sum():3d} {early[w].mean() if w.sum() else -1:5.2f} {np.median(rres[w]) if w.sum() else -1:4.0f}'
        print(l)
# all triggered events: late-crossing class only, what is the mesh height lower edge?
print('\nmesh peak by first-crossing class, all events: percentiles 1 5 16 50')
x = g('m60'); r = np.where(x >= 0, x - g('icbin'), 1e6)
for lab, w in [('crossing < IC-100', r < -100), ('within -100..+100 of IC', (r >= -100) & (r < 100)), ('late +100..+900', (r >= 100) & (r < 1e5)), ('no crossing', r > 1e5)]:
    print(f'  {lab:28s} n {w.sum():7d} frac {w.mean():.3f}', np.round(np.percentile(g('mpk')[w], [1, 5, 16, 50])))
