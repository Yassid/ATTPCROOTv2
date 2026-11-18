# 2026-11-18 copy of dt_trig4.py for the INTEGRAL-GATE production (trig/dt_hdbint, g.s. shift +0.2897)
# (d,t): rerun the vz-flatness test with Ebeam(vz) (the (d,p) beam-loss line, pivoted so vz 312 = 184.0) to separate kinematic migration from a loss.
import numpy as np, importlib.util, sys
from scipy.stats import chi2 as C2
spec = importlib.util.spec_from_file_location('t1', '/dev/null'); 
U = 931.49401; M1, M2, M3, M4 = 16.0147013*U, 2.0135532*U, 3.01550072*U, 15.0105993*U
def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))
def kine(Kp, th, Ke):
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
D = '/mnt/f/a1975/sim_rebuild/trig/dt_hdbint/'; z = np.load(D + 'joined_dt.npz'); T = z['T']; base = z['base']; early = z['early']; ke, th, vz = T[:, 3], T[:, 4], T[:, 5]
ex, tcm = kine(184.0 - 0.011071*(vz - 312.5), np.radians(th - 0.046*ke), ke); ex = ex + 0.2897
sel = base & np.isfinite(ex)
ZB = np.array([25, 125, 225, 325, 425, 525, 600.]); dz = np.diff(ZB)
print('median Ex by vz slice, Ebeam(vz):')
for nm, (a, b) in {'0.740': (0.3, 1.3), '3.103': (2.6, 3.7), '4.657': (4.1, 5.3), '6.700': (6.0, 7.5)}.items():
    print(f'  {nm:6s}', ' '.join(f'{np.median(ex[sel & (ex > a) & (ex < b) & (vz >= z0) & (vz < z1)]):6.3f}' for z0, z1 in zip(ZB[:-1], ZB[1:])))
def test(w, lab):
    n = np.histogram(vz[w], ZB)[0].astype(float); lam = n.sum()*dz/dz.sum(); c = np.sum((n - lam)**2/np.maximum(lam, 1))
    x = .5*(ZB[1:] + ZB[:-1]); y = n/dz; s = np.sqrt(np.maximum(n, 1))/dz; W = 1/s**2; xm = np.sum(W*x)/W.sum(); ym = np.sum(W*y)/W.sum(); b = np.sum(W*(x - xm)*(y - ym))/np.sum(W*(x - xm)**2); be = 1/np.sqrt(np.sum(W*(x - xm)**2))
    print(f'  {lab:14s} N {int(n.sum()):5d} per-100mm', ' '.join(f'{v*100:6.1f}' for v in n/dz), f'| chi2 {c:5.1f}/5 p {C2.sf(c, 5):.3f} | slope over 575 mm {b*575/ym:+.2f}({be*575/ym:.2f})')
L = {'g.s.+0.740': (-1, 1.3), '3.103': (2.6, 3.7), '4.657': (4.1, 5.3), '6.700': (6.0, 7.5), 'cont 8-14': (8, 14)}
TB = [5, 10, 15, 20, 25, 30, 35, 40, 45]; allc = []; alln = 0
for nm, (a, b) in L.items():
    wl = (ex > a) & (ex < b); print(f'== {nm}'); test(sel & wl & (tcm < 40), 'theta_cm<40')
    for p, q in zip(TB[:-1], TB[1:]): test(sel & wl & (tcm >= p) & (tcm < q), f'{p}-{q}')
