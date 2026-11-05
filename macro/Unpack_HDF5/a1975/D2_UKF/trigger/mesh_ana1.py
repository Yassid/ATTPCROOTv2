# Event-level mesh study, all runs with features: scalers, IC timing, mesh crossing time and amplitude, join to tracks, GET vs mesh timing.
import numpy as np, glob, os
D = '/mnt/f/a1975/sim_rebuild/trig/'
fs = sorted(glob.glob(D + 'mesh/run_*.npz')); names = None; rows = []; sc = []
for f in fs:
    z = np.load(f); names = list(z['names']); d = z['d']; run = int(os.path.basename(f)[4:8]); rows.append(np.c_[np.full(len(d), run), d]); s = f.replace('run_', 'scaler_').replace('.npz', '.npy')
    if os.path.exists(s): sc.append(np.load(s))
A = np.vstack(rows); col = {n: i + 1 for i, n in enumerate(names)}; col['run'] = 0; g = lambda n: A[:, col[n]]
print('runs', len(fs), 'events', len(A)); S = np.sum(sc, axis=0)
print('scalers: ic %.4g mesh %.4g trig_free %.4g trig_live %.4g | mesh/ic %.4f trig/mesh %.3f trig/ic %.5f live %.3f' % (S[4], S[5], S[2], S[3], S[5]/S[4], S[2]/S[5], S[2]/S[4], S[3]/S[2]))
ib = g('icbin'); print('IC peak bin p1 p16 p50 p84 p99', np.percentile(ib, [1, 16, 50, 84, 99]), 'IC amp p16 p50 p84', np.percentile(g('icamp'), [16, 50, 84]), 'IC pulses per trace', np.percentile(g('icn'), [16, 50, 84]))
pk = g('mpk'); print('mesh peak ADC percentiles 1 5 10 20 30 50 70 90 99', np.round(np.percentile(pk, [1, 5, 10, 20, 30, 50, 70, 90, 99])))
h, e = np.histogram(pk, bins=[0, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 600, 800, 1000, 1500, 2000, 3000]); print('mesh peak hist'); print(' '.join(f'{int(a)}-{int(b)}:{c}' for a, b, c in zip(e[:-1], e[1:], h)))
for t in ['m60', 'm100', 'm200']:
    x = g(t); ok = x >= 0; rel = x[ok] - ib[ok]; print(t, 'present', round(ok.mean(), 3), 'first crossing minus IC bin, percentiles 1 5 16 50 84 95 99', np.round(np.percentile(rel, [1, 5, 16, 50, 84, 95, 99])))
x = g('m60'); ok = x >= 0; rel = x[ok] - ib[ok]; h, e = np.histogram(rel, bins=np.arange(-1200, 1000, 100)); print('m60 - IC per 100 bins:'); print(' '.join(f'{int(a)}:{c}' for a, c in zip(e[:-1], h)))
print('mesh pulse already above 60 before bin 900:', round(np.mean(g('m0900') > 60), 3))
# tracks
T = np.loadtxt(D + 'tracks.txt'); key = {(int(r), int(ev)): i for i, (r, ev) in enumerate(zip(A[:, 0], g('ev')))}; idx = np.array([key.get((int(r), int(e) + (1 if int(r) >= 31 else 0)), -1) for r, e in zip(T[:, 0], T[:, 1])]); ok = idx >= 0; print('tracks', len(T), 'joined', ok.sum(), 'with direct ev=evt')
T = T[ok]; idx = idx[ok]; c = dict(run=0, ev=1, gate=2, th=3, vz=4, ke=5, ex=6, tc=7, zmin=8, zmax=9, nh=10, rzmin=11, rzmax=12, dirB=13)
sel = (T[:, c['gate']] == 1) & (T[:, c['dirB']] == 1) & (T[:, c['rzmin']] < 30); t = T[sel]; i = idx[sel]
m = A[i, col['m60']]; ic = A[i, col['icbin']]; good = m >= 0; rel = (m - ic)[good]; zm = t[good, c['zmin']]
# expected first arrival: nearest end of the track to the pad plane, z_digi / v, in 80 ns bins, v = 1.10424 cm/us
pred = zm / (11.0424 * 0.08); res = rel - pred; print('\ngated backward tracks', len(t), 'with mesh crossing', good.sum(), 'of which no crossing', (~good).sum())
print('mesh first crossing - predicted from GET track start, bins: p5 p16 p50 p84 p95', np.round(np.percentile(res, [5, 16, 50, 84, 95])), ' (1 bin = 0.88 mm)')
print('by vz:'); 
for a, b in [(100,200),(200,300),(300,400),(400,500),(500,600),(600,700),(700,800),(800,900),(900,1000)]:
    w = (t[good, c['vz']] >= a) & (t[good, c['vz']] < b)
    if w.sum() > 10: print(f'  vz {a}-{b}: n {w.sum():5d} residual p16 p50 p84', np.round(np.percentile(res[w], [16, 50, 84])), ' mesh peak p16 p50 p84', np.round(np.percentile(A[i[good], col["mpk"]][w], [16, 50, 84])))
np.savez(D + 'joined.npz', t=t, i=i)
