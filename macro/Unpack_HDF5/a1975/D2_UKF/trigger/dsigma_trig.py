# Corrected (d,p) bound-state angular distribution: genfit r280 data and sim, sim-anchored Ebeam, data Ex smear, trigger model on each sim event. dsigma per vz range, whole range, against the FRESCO curves.
import numpy as np, sys
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'; U = 931.49401; M1, M2, M3, M4 = 16.0147013 * U, 2.0135532 * U, 1.00727646688 * U, 17.0225864 * U
def om(x, y, z): return np.sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z)
def kine(Kp, th, Ke):
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th) * om(s, M1*M1, M2*M2) * om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2) * (M2*M2 + M3*M3 - uu)) / (2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4)) / (om(s, M1*M1, M2*M2) * om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
Eb = lambda z: 183.871 - 0.011071 * z; ZOFF = 19.1; EL, EH = -1.5, 0.729; SM = 0.238
keys = np.loadtxt(D + 'keys.txt'); kmap = {(round(a, 3), round(b, 3)): int(e) for e, a, b in keys[:, :3]}; gz = keys[:, 1] + ZOFF; gx, gc = kine(Eb(gz), np.radians(keys[:, 2]), keys[:, 3])
Q = np.loadtxt(D + 'qhits/q_all.txt'); qm = {int(r[1]): r[2:] for r in Q}
g = np.loadtxt(D + 'gf_sim.txt'); ev = np.array([kmap.get((round(a, 3), round(b, 3)), -1) for a, b in g[:, [6, 5]]]); print('sim rows', len(g), 'unmatched', (ev < 0).sum())
ok = (ev >= 0) & np.array([e in qm for e in ev]); g = g[ok]; ev = ev[ok]; Pq = np.array([qm[e] for e in ev]); sx, sc = kine(Eb(g[:, 2]), np.radians(g[:, 1]), g[:, 0]); sel = (g[:, 3] <= 13.5) & (g[:, 1] >= 89) & (g[:, 1] <= 161) & (g[:, 0] >= 1.2) & np.isfinite(sx); pwin = ndtr((EH - sx) / SM) - ndtr((EL - sx) / SM)
d = np.loadtxt(D + 'gf_data.txt'); dx, dc = kine(Eb(d[:, 2]), np.radians(d[:, 1]), d[:, 0]); ds = (d[:, 3] <= 13.5) & (d[:, 1] >= 89) & (d[:, 1] <= 161) & (d[:, 0] >= 1.2) & np.isfinite(dx)
TB = np.array([10, 15, 20, 25, 30, 35, 40.]); dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))
def rd(f): return np.loadtxt(f, comments=['#', '&', '@'])
base = '/home/yassid/fair_install/16C_dp/fresco/outputs/dp17C_bound_adwa_state%d.dat'; G = [rd(base % i) for i in (1, 2, 3)]
def avg(g, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w)
TH = np.array([[avg(g_, a, b) for g_ in G] for a, b in zip(TB[:-1], TB[1:])]); L2 = TH[:, 0] + TH[:, 2]; L0 = TH[:, 1]
def dsig(z0, z1, ptrig):
    # background from the sideband, same vz range
    w = ds & (d[:, 2] >= z0) & (d[:, 2] < z1); Y = np.histogram(dc[w & (dx >= EL) & (dx <= EH)], TB)[0].astype(float); B = np.histogram(dc[w & (dx >= -4) & (dx < EL)], TB)[0].astype(float) * (EH - EL) / 2.5
    Ng = np.histogram(gc[(gz >= z0) & (gz < z1) & np.isfinite(gx)], TB)[0].astype(float); w = sel & (g[:, 2] >= z0) & (g[:, 2] < z1); Nr = np.histogram(sc[w], TB, weights=(pwin * ptrig)[w])[0]; e = Nr / np.maximum(Ng, 1); L = 316.98 * (z1 - z0) / 971.73
    x = (Y - B) / np.maximum(e, 1e-9) / dO / L; er = np.sqrt(np.maximum(Y, 1) + B) / np.maximum(e, 1e-9) / dO / L; return x, er, Y, B, e
def fitS(x, er, m, use): S = np.sum((x * m / er**2)[use]) / np.sum((m * m / er**2)[use]); c = np.sum((((x - S * m) / er)**2)[use]); return S, c
print('theta_cm                    10-15  15-20  20-25  25-30  30-35  35-40; DWBA l=2 unit S', np.round(L2, 1), ' l=0', np.round(L0, 2))
PARS = [('none', None), ('best s10 T31544 w.3 f.1', (0, 31544, .3, .1)), ('s20 T56252 w.3 f.1', (1, 56252, .3, .1)), ('bound alone s10 T25588 w.1 f.1', (0, 25588, .1, .1)), ('Ex3-5 alone s10 T43173 w.8 f0', (0, 43173, .8, 0.)), ('s10 T31544 w.3 f.05', (0, 31544, .3, .05)), ('s10 T31544 w.3 f.2', (0, 31544, .3, .2))]
for lab, pr in PARS:
    p = np.ones(len(g)) if pr is None else pr[3] + (1 - pr[3]) * ndtr(np.log(Pq[:, pr[0]] / pr[1]) / pr[2]); print('\n===', lab)
    Xs = []
    for z0, z1 in [(200, 400), (400, 600), (600, 800), (800, 900), (200, 900)]:
        x, er, Y, B, e = dsig(z0, z1, p); use = np.ones(6, bool); Sl2, cl2 = fitS(x, er, L2, use); Sl0, cl0 = fitS(x, er, L0, use); Xs.append((x, er))
        print(f' vz {z0}-{z1} dsdo', ' '.join(f'{a:6.2f}({b:4.2f})' for a, b in zip(x, er)), f'| l=2 S {Sl2:.3f} chi2/5 {cl2:6.1f}  l=0 S {Sl0:5.2f} chi2/5 {cl0:6.1f}'); print('            eps ', ' '.join(f'{a:6.3f}      ' for a in e), ' Y', Y.astype(int))
    # do the vz ranges agree? chi2 of the four ranges around their weighted mean, per angle
    A = np.array([q[0] for q in Xs[:4]]); E = np.array([q[1] for q in Xs[:4]]); wm = np.sum(A / E**2, axis=0) / np.sum(1 / E**2, axis=0); print(' vz-range consistency chi2 per angle (3 dof each):', np.round(np.sum(((A - wm) / E)**2, axis=0), 1), 'total', round(np.sum(((A - wm) / E)**2), 1), '/ 18')
