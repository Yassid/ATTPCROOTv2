# Trigger-safe fiducial: keep tracks whose heavy residual is predicted outside the beam hole at the pad plane, where the mesh fires on the residual whatever the proton does. Same cut on generated, reconstructed and data events. No trigger model. Stability against the cut is the test.
import numpy as np
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'; U = 931.49401; M1, M2, M3, M4 = 16.0147013 * U, 2.0135532 * U, 1.00727646688 * U, 17.0225864 * U
def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))
def kine(Kp, th, Ke):
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th) * om(s, M1*M1, M2*M2) * om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2) * (M2*M2 + M3*M3 - uu)) / (2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4)) / (om(s, M1*M1, M2*M2) * om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
Eb = lambda z: 183.871 - 0.011071 * z; ZOFF = 19.1; EL, EH = -1.5, 0.729; SM = 0.238
def rres(ke, th, vz):
    t = np.radians(th); p = np.sqrt(ke**2 + 2 * ke * 938.272); E = Eb(vz); pb = np.sqrt(E**2 + 2 * E * M1); return (1000 - vz) * np.tan(np.arctan2(p * np.sin(t), pb - p * np.cos(t)))
K = np.loadtxt(D + 'keys.txt'); gz = K[:, 1] + ZOFF; gx, gc = kine(Eb(gz), np.radians(K[:, 2]), K[:, 3]); gr = rres(K[:, 3], K[:, 2], gz)
g = np.loadtxt(D + 'gf_sim.txt'); sx, sc = kine(Eb(g[:, 2]), np.radians(g[:, 1]), g[:, 0]); sr = rres(g[:, 0], g[:, 1], g[:, 2]); sel = (g[:, 3] <= 13.5) & (g[:, 1] >= 89) & (g[:, 1] <= 161) & (g[:, 0] >= 1.2) & np.isfinite(sx); pwin = ndtr((EH - sx) / SM) - ndtr((EL - sx) / SM)
d = np.loadtxt(D + 'gf_data.txt'); dx, dc = kine(Eb(d[:, 2]), np.radians(d[:, 1]), d[:, 0]); dr = rres(d[:, 0], d[:, 1], d[:, 2]); ds = (d[:, 3] <= 13.5) & (d[:, 1] >= 89) & (d[:, 1] <= 161) & (d[:, 0] >= 1.2) & np.isfinite(dx)
TB = np.array([10, 15, 20, 25, 30, 35, 40, 45.]); dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:]))); LPM = 316.98 / 971.73
base = '/home/yassid/fair_install/16C_dp/fresco/outputs/dp17C_bound_adwa_state%d.dat'; G = [np.loadtxt(base % i, comments=['#', '&', '@']) for i in (1, 2, 3)]
def avg(g_, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g_[:, 0], g_[:, 1]) * w) / np.sum(w)
TH = np.array([[avg(g_, a, b) for g_ in G] for a, b in zip(TB[:-1], TB[1:])]); L2 = TH[:, 0] + TH[:, 2]; L0 = TH[:, 1]
ZS = np.arange(150, 915, 10.)
def run(R, z0=200, z1=900):
    w = ds & (d[:, 2] >= z0) & (d[:, 2] < z1) & (dr >= R); Y = np.histogram(dc[w & (dx >= EL) & (dx <= EH)], TB)[0].astype(float); B = np.histogram(dc[w & (dx >= -4) & (dx < EL)], TB)[0].astype(float) * (EH - EL) / 2.5
    gsel = (gz >= z0) & (gz < z1) & np.isfinite(gx); Nall = np.histogram(gc[gsel], TB)[0].astype(float); Ng = np.histogram(gc[gsel & (gr >= R)], TB)[0].astype(float); s = sel & (g[:, 2] >= z0) & (g[:, 2] < z1) & (sr >= R); Nr = np.histogram(sc[s], TB, weights=pwin[s])[0]; e = Nr / np.maximum(Ng, 1)
    Leff = LPM * (z1 - z0) * Ng / np.maximum(Nall, 1); x = (Y - B) / np.maximum(e, 1e-6) / dO / np.maximum(Leff, 1e-6); er = np.sqrt(np.maximum(Y, 1) + B) / np.maximum(e, 1e-6) / dO / np.maximum(Leff, 1e-6); return x, er, Y, B, e, Leff
print('theta_cm 10-15 15-20 20-25 25-30 30-35 35-40 40-45. DWBA unit S: l=2', np.round(L2, 1), 'l=0', np.round(L0, 2))
for R in (0, 16, 20, 22, 24, 26, 28, 30, 34):
    x, er, Y, B, e, Le = run(R); use = (Y >= 3) & (e > .3); 
    S = np.sum((x * L2 / er**2)[use]) / np.sum((L2 * L2 / er**2)[use]) if use.any() else 0; c = np.sum((((x - S * L2) / er)**2)[use]); S0 = np.sum((x * L0 / er**2)[use]) / np.sum((L0 * L0 / er**2)[use]) if use.any() else 0; c0 = np.sum((((x - S0 * L0) / er)**2)[use])
    print(f'r_res>={R:2d} dsdo', ' '.join(f'{a:5.2f}({b:4.2f})' if Y[i] >= 3 else '    -      ' for i, (a, b) in enumerate(zip(x, er))), f'| l=2 S {S:.3f} chi2 {c:5.1f}/{use.sum()-1} l=0 {S0:5.2f} chi2 {c0:5.1f}'); print('          Y ', ' '.join(f'{int(a):5d}      ' for a in Y), ' eps', np.round(e, 2)); print('          Leff', ' '.join(f'{a:5.1f}      ' for a in Le))
print('\nprevious numbers at 12.5..37.5: notebook far-slice 23.1 16.8 3.7 3.7 4.2 3.4; adopted -, 1.56, 3.68 5.85 6.82 3.49')
