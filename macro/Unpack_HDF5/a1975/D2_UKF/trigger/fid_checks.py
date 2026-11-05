# Stress tests of the r_res fiducial: (a) Spyral chain, independent of genfit: data/sim per theta_cm vs cut, first and second half of the runs; (b) genfit with the pivot Ebeam frame; (c) tighter chi2, wider Ex window.
import numpy as np, runpy, io, contextlib
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'
with contextlib.redirect_stdout(io.StringIO()): G_ = runpy.run_path(D + 'dsigma_fid.py')
rres, kine = G_['rres'], G_['kine']
T = np.loadtxt(D + 'tracks.txt'); T = T[T[:, 2] == 1]; S = np.loadtxt(D + 'sim_tracks.txt'); S = S[S[:, 0] == 0]
drr = rres(T[:, 5], T[:, 3], T[:, 4]); srr = rres(S[:, 4], S[:, 2], S[:, 3]); TB = np.array([20, 25, 30, 35, 40.]); runs = np.unique(T[:, 0]); h1 = np.isin(T[:, 0], runs[:len(runs) // 2])
def ratio(mask, R):
    w = mask & (T[:, 6] > -1.5) & (T[:, 6] < 0.729) & (T[:, 4] >= 200) & (T[:, 4] < 900) & (drr >= R); Y = np.histogram(T[w, 7], TB)[0].astype(float); w2 = (S[:, 3] >= 200) & (S[:, 3] < 900) & (srr >= R); N = np.histogram(S[w2, 6], TB)[0].astype(float); return Y, N
print('(a) Spyral chain: data/sim normalised to 35-40, with stat error; DWBA l=2 sum ratios to 35-40: 20-25 4.91  25-30 3.02  30-35 1.72')
for lab, m in [('all', np.ones(len(T), bool)), ('first half runs', h1), ('second half', ~h1)]:
    for R in (0, 20, 24, 26, 28, 30):
        Y, N = ratio(m, R); r = np.where(N > 0, Y / np.maximum(N, 1), 0); nr = r[-1] if r[-1] > 0 else 1; o = ' '.join((f'{r[i] / nr:5.2f}({r[i] / nr * np.sqrt(1 / max(Y[i], 1) + 1 / max(Y[-1], 1) + 1 / max(N[i], 1) + 1 / max(N[-1], 1)):4.2f})' if Y[i] >= 3 else '     -     ') for i in range(4)); print(f'  {lab:16s} R>={R:2d} {o}  Y {Y.astype(int)}')
# (b) pivot frame, by redefining the beam energy in the genfit fiducial
g, d = G_['g'], G_['d']; EL, EH, SM = G_['EL'], G_['EH'], G_['SM']; K = np.loadtxt(D + 'keys.txt'); TB2 = np.array([25, 30, 35, 40.]); dO = 2 * np.pi * (np.cos(np.radians(TB2[:-1])) - np.cos(np.radians(TB2[1:]))); L2 = G_['L2'][3:6]
def calc(Eb, R, c2=13.5, top=EH, smear=SM):
    gz = K[:, 1] + 19.1; gx, gc = kine(Eb(gz), np.radians(K[:, 2]), K[:, 3]); gr = rres(K[:, 3], K[:, 2], gz); sx, sc = kine(Eb(g[:, 2]), np.radians(g[:, 1]), g[:, 0]); sr = rres(g[:, 0], g[:, 1], g[:, 2]); dx, dc = kine(Eb(d[:, 2]), np.radians(d[:, 1]), d[:, 0]); dr = rres(d[:, 0], d[:, 1], d[:, 2])
    ds = (d[:, 3] <= c2) & (d[:, 1] >= 89) & (d[:, 1] <= 161) & (d[:, 0] >= 1.2) & (d[:, 2] >= 200) & (d[:, 2] < 900) & (dr >= R); ss = (g[:, 3] <= c2) & (g[:, 1] >= 89) & (g[:, 1] <= 161) & (g[:, 0] >= 1.2) & (g[:, 2] >= 200) & (g[:, 2] < 900) & (sr >= R); pw = ndtr((top - sx) / smear) - ndtr((EL - sx) / smear)
    Y = np.histogram(dc[ds & (dx >= EL) & (dx <= top)], TB2)[0].astype(float); B = np.histogram(dc[ds & (dx >= -4) & (dx < EL)], TB2)[0].astype(float) * (top - EL) / 2.5; gs = (gz >= 200) & (gz < 900) & np.isfinite(gx); Na = np.histogram(gc[gs], TB2)[0].astype(float); Ng = np.histogram(gc[gs & (gr >= R)], TB2)[0].astype(float); Nr = np.histogram(sc[ss], TB2, weights=pw[ss])[0]; e = Nr / np.maximum(Ng, 1); Le = G_['LPM'] * 700 * Ng / np.maximum(Na, 1); x = (Y - B) / np.maximum(e, 1e-6) / dO / np.maximum(Le, 1e-6); er = np.sqrt(np.maximum(Y, 1) + B) / np.maximum(e, 1e-6) / dO / np.maximum(Le, 1e-6); u = Y >= 3; s = np.sum((x * L2 / er**2)[u]) / np.sum((L2 * L2 / er**2)[u]); c = np.sum((((x - s * L2) / er)**2)[u]); return x, er, s, c, Y
EB = G_['Eb']; PV = lambda z: 184.25 - 0.011030 * (z - 500)
print('\n(b,c) genfit, dsdo 25-30 30-35 35-40, l=2 S and chi2/2')
for lab, E, R, c2, top, sm in [('nominal R28', EB, 28, 13.5, EH, SM), ('nominal R26', EB, 26, 13.5, EH, SM), ('pivot frame R28', PV, 28, 13.5, EH, SM), ('pivot frame R26', PV, 26, 13.5, EH, SM), ('chi2<5 R28', EB, 28, 5, EH, SM), ('chi2<1 R28', EB, 28, 1, EH, SM), ('Ex top 1.0 R28', EB, 28, 13.5, 1.0, SM), ('smear .15 R28', EB, 28, 13.5, EH, .15), ('smear .30 R28', EB, 28, 13.5, EH, .30)]:
    x, er, s, c, Y = calc(E, R, c2, top, sm); print(f'  {lab:18s}', ' '.join(f'{a:5.2f}({b:4.2f})' for a, b in zip(x, er)), f' S {s:.3f} chi2 {c:5.2f}  Y', Y.astype(int))
