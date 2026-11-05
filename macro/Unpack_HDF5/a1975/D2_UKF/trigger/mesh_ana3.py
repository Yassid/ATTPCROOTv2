# Calibration: measured mesh peak against the pad-charge pulse model, data events. Which width describes the mesh, and what does the hardware edge correspond to in P units?
import numpy as np, glob, os
D = '/mnt/f/a1975/sim_rebuild/trig/'; A = {}; names = None; SG = [10, 20, 30, 40, 60, 80]; out = []
for f in sorted(glob.glob(D + 'mesh/run_*.npz')):
    z = np.load(f); names = list(z['names']); run = int(os.path.basename(f)[4:8]); q = D + f'qhits/q_{run}.txt'
    if f'done qhits/q_{run}.txt' not in open(D + 'run_hits.log').read(): continue
    d = z['d']; idx = {int(e): i for i, e in enumerate(d[:, 0])}
    for r in np.loadtxt(q, ndmin=2):
        i = idx.get(int(r[1]) + (1 if run >= 31 else 0));
        if i is not None: out.append(np.r_[run, d[i], r[2:]])
X = np.array(out); c = {n: i + 1 for i, n in enumerate(names)}; n0 = len(names) + 1; P = X[:, n0:n0 + 6]; print('events with both', len(X))
pk = X[:, c['mpk']]; rel = np.where(X[:, c['m60']] >= 0, X[:, c['m60']] - X[:, c['icbin']], 1e6); sat = pk > 1700
for lab, w in [('late fired, unsaturated', (rel >= 100) & (rel < 1e5) & ~sat), ('early fired, unsaturated', (rel >= -100) & (rel < 100) & ~sat)]:
    print(lab, 'n', w.sum())
    for k, s in enumerate(SG):
        l = np.log(P[w, k]); y = np.log(pk[w]); co = np.corrcoef(l, y)[0, 1]; r = np.median(pk[w] / P[w, k]); print(f'  sigma {s:3d} corr(log) {co:.3f}  median mesh/P {r:.5f}  scatter of ln ratio (p84-p16)/2 {(np.percentile(y - l, 84) - np.percentile(y - l, 16)) / 2:.3f}  mesh edge 320 ADC = P {320 / r:9.0f}')
print('\nP (sigma 10) percentiles of triggered data events: 1 5 16 50 84', np.round(np.percentile(P[:, 0], [1, 5, 16, 50, 84])))
print('no mesh crossing events among these', (rel > 1e5).sum(), 'median P sigma10', np.median(P[rel > 1e5, 0]) if (rel > 1e5).any() else -1)
S = np.loadtxt(D + 'qhits/q_s0.txt'); print('sim g.s. P (sigma 10) percentiles 1 5 16 50 84', np.round(np.percentile(S[:, 2], [1, 5, 16, 50, 84])))
