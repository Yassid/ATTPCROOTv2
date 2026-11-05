# Is event N of the reco the same event as evt N in the FRIB waveforms? Run by run: correlation of the measured mesh peak with the pad-charge pulse, at offsets of the event index. IC too.
import numpy as np, glob, os
D = '/mnt/f/a1975/sim_rebuild/trig/'; res = []
for f in sorted(glob.glob(D + 'mesh/run_*.npz')):
    z = np.load(f); names = list(z['names']); run = int(os.path.basename(f)[4:8]); d = z['d']; q = D + f'qhits/q_{run}.txt'
    if not os.path.exists(q) or os.path.getsize(q) == 0: continue
    Q = np.loadtxt(q, ndmin=2); pk = d[:, names.index('mpk')]; idx = {int(e): i for i, e in enumerate(d[:, 0])}; out = []
    for off in (-3, -2, -1, 0, 1, 2, 3):
        j = np.array([idx.get(int(e) + off, -1) for e in Q[:, 1]]); w = j >= 0; x = np.log(np.maximum(Q[w, 3], 1)); y = pk[j[w]]; s = y < 1700; out.append(np.corrcoef(x[s], np.log(np.maximum(y[s], 1)))[0, 1] if s.sum() > 5 else 0)
    j = np.array([idx.get(int(e), -1) for e in Q[:, 1]]); w = j >= 0; nocross = (d[j[w], names.index('m60')] < 0).mean(); res.append((run, len(Q), len(d), int(d[:, 0].max()), nocross, *out))
print('run   nev_q  nFRIB maxev  noCross | corr at offset -3 -2 -1 0 +1 +2 +3')
for r in res: print(f'{r[0]:4d} {r[1]:6d} {r[2]:6d} {r[3]:6d}  {r[4]:.3f} | ' + ' '.join(f'{v:5.2f}' for v in r[5:]) + ('   <- off' if abs(r[8]) < .4 else ''))
