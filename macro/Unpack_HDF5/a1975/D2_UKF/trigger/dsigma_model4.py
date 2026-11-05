# Corrected (d,p) bound angular distribution with the two-path trigger model; per vz range, consistency between ranges, FRESCO fits, parameter variants.
import numpy as np, runpy, io, contextlib, sys
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'
with contextlib.redirect_stdout(io.StringIO()): G_ = runpy.run_path(D + 'dsigma_fid.py')
kine, Eb, rres = G_['kine'], G_['Eb'], G_['rres']; g, sc, sx, sr, sel, pwin, d, dx, dc, ds, gz, gx, gc = [G_[k] for k in 'g sc sx sr sel pwin d dx dc ds gz gx gc'.split()]; EL, EH = G_['EL'], G_['EH']; LPM = G_['LPM']; L2, L0, TH = G_['L2'], G_['L0'], G_['TH']
keys = np.loadtxt(D + 'keys.txt'); kmap = {(round(a, 3), round(b, 3)): int(e) for e, a, b in keys[:, :3]}; ev = np.array([kmap.get((round(a, 3), round(b, 3)), -1) for a, b in g[:, [6, 5]]]); qm = {int(r[1]): r[2:] for r in np.loadtxt(D + 'qhits/qp_all.txt')}; Pq = np.array([qm.get(e, np.full(9, 1.)) for e in ev]); bad = ev < 0
TB = np.array([10, 15, 20, 25, 30, 35, 40.]); nb = 6; dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:]))); L2 = L2[:nb]; L0 = L0[:nb]; TH = TH[:nb]
def dsig(z0, z1, p):
    w = ds & (d[:, 2] >= z0) & (d[:, 2] < z1); Y = np.histogram(dc[w & (dx >= EL) & (dx <= EH)], TB)[0].astype(float); B = np.histogram(dc[w & (dx >= -4) & (dx < EL)], TB)[0].astype(float) * (EH - EL) / 2.5
    Ng = np.histogram(gc[(gz >= z0) & (gz < z1) & np.isfinite(gx)], TB)[0].astype(float); s = sel & ~bad & (g[:, 2] >= z0) & (g[:, 2] < z1); Nr = np.histogram(sc[s], TB, weights=(pwin * p)[s])[0]; e = Nr / np.maximum(Ng, 1); L = LPM * (z1 - z0); x = (Y - B) / np.maximum(e, 1e-9) / dO / L; er = np.sqrt(np.maximum(Y, 1) + B) / np.maximum(e, 1e-9) / dO / L; return x, er, Y, e
def fit1(x, er, m, u): S = np.sum((x * m / er**2)[u]) / np.sum((m * m / er**2)[u]); return S, np.sqrt(1 / np.sum((m * m / er**2)[u])), np.sum((((x - S * m) / er)**2)[u])
PAR = eval(sys.argv[1]) if len(sys.argv) > 1 else []
for lab, k, T, w, R0, s_r, e0 in PAR:
    p = 1 - (1 - e0) * (1 - ndtr(np.log(Pq[:, k] / T) / w)) * (1 - ndtr((sr - R0) / s_r)); print('\n===', lab); A = []; E = []
    for z0, z1 in [(200, 400), (400, 600), (600, 800), (800, 900), (200, 900)]:
        x, er, Y, e = dsig(z0, z1, p); u = np.ones(nb, bool); u20 = TB[:-1] >= 20; a = fit1(x, er, L2, u); b = fit1(x, er, L2, u20); c = fit1(x, er, L0, u)
        print(f' vz {z0}-{z1}', ' '.join(f'{v:6.2f}({q:4.2f})' for v, q in zip(x, er)), f'| l=2 10-40 S {a[0]:.3f} chi2 {a[2]:5.1f}/5; 20-40 S {b[0]:.3f} chi2 {b[2]:5.1f}/3; l=0 chi2 {c[2]:5.1f}/5'); print('        eps', ' '.join(f'{v:6.3f}      ' for v in e), 'Y', Y.astype(int))
        if z1 - z0 < 700: A.append(x); E.append(er)
    A = np.array(A); E = np.array(E); wm = np.sum(A / E**2, axis=0) / np.sum(1 / E**2, axis=0); cs = np.sum(((A - wm) / E)**2, axis=0); print(' vz ranges consistency per angle (3 dof each)', np.round(cs, 1), 'total', round(cs.sum(), 1), '/ 18'); print(' weighted mean of the four ranges', np.round(wm, 2), ' +-', np.round(1 / np.sqrt(np.sum(1 / E**2, axis=0)), 2))
