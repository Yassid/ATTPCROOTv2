# 16C(d,p)17C unbound states: angular distributions on the ADOPTED 2026-11-16 products (HDBSCAN + float tau, data calibrated).
# NO TRIGGER CORRECTION: the fiducial cut that makes the bound states trigger-safe empties Ex 2-8 (C17_BOUND_STATES.md), so these
# are (yield / genfit+Spyral efficiency / luminosity) only. The vz-slice comparison is the test of whether the trigger distorts them.
#   ~/Spyral/venv/bin/python -I unb_ad.py [vzLo vzHi tag]
# Spectrum model per theta_cm bin (fit_unbound_r280.C model, Poisson likelihood, 100 keV, Ex 1-5): 1n phase space
# a*sqrt(Ex-0.729)*(1+t(Ex-3)) + 3 gaussians; positions, common sigma and tilt t FIXED from the fit of the whole vz window, a and the
# three amplitudes free per bin. Efficiency per reco theta_cm bin from the Ex 3.9 sim (r280_hdbtau), applied to all three states.
import sys, os, numpy as np, uproot
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("K = np.loadtxt")[0])   # kine, Eb, ZOFF, rres, masses
T = '/mnt/f/a1975/sim_rebuild/trig_hdbtau/'; U = '/mnt/f/a1975/sim_rebuild/unb_hdbtau/'; R5 = '/mnt/f/a1975/sim_rebuild/r280_hdbtau/'
CAL = (184.25, 11.0, 0.18, -3, -0.05, -0.6947); LPM = 316.98 / 971.73; ENV = 0.84
TLO, THI = 89, 113                     # theta_lab (corrected) window: every track can reach Ex 5.36 (kinematic ceiling)
TB = np.array([20, 25, 30, 35, 40.]); dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))
EXLO, EXHI, BW = 1.0, 5.0, 0.1; edges = np.arange(EXLO, EXHI + 1e-9, BW); cen = .5 * (edges[1:] + edges[:-1])
vz0, vz1, tag = (float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]) if len(sys.argv) > 3 else (100., 900., 'main')

d = np.loadtxt(T + 'gf_data.txt'); ke, th, vz, c2 = d.T.copy()
ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); dx, dc = kine(CAL[0] - CAL[1] / 1000 * (vz - 500), np.radians(th), ke); dx = dx + CAL[5]
dsel = (c2 <= 13.5) & (th >= TLO) & (th < THI) & (ke >= 1.2) & np.isfinite(dx)

def model(p, n):   # p: a, t, sigma, then (amp, mean) per peak; returns counts per bin (bin-centre approximation, 100 keV)
   x = cen; m = p[0] * np.sqrt(np.maximum(x - 0.729, 0)) * (1 + p[1] * (x - 3))
   for k in range(n): m = m + p[3 + 2 * k] * np.exp(-.5 * ((x - p[4 + 2 * k]) / p[2]) ** 2)
   return np.maximum(m, 1e-9)
def nll(p, h, n): m = model(p, n); return 2 * np.sum(m - h + np.where(h > 0, h * np.log(np.maximum(h, 1e-12) / m), 0))

# 1) global fit, this vz window, all theta: positions, sigma, tilt
w = dsel & (vz >= vz0) & (vz < vz1); hg = np.histogram(dx[w], edges)[0].astype(float)
p0 = [hg.mean() / 1.5, .2, .22, 20, 2.6, 30, 3.45, 25, 4.0]
bnd = [(0, None), (-.6, 2), (.08, .9), (0, None), (2.2, 3.0), (0, None), (3.0, 3.8), (0, None), (3.7, 4.6)]
G = minimize(nll, p0, args=(hg, 3), bounds=bnd, method='L-BFGS-B'); pg = G.x
print(f'=== {tag}: vz {vz0:.0f}-{vz1:.0f}, theta_lab {TLO}-{THI}, {int(hg.sum())} counts in Ex {EXLO}-{EXHI}; global -2lnL {G.fun:.1f}/{len(cen) - 9}')
print(f'    Ex {pg[4]:.3f} {pg[6]:.3f} {pg[8]:.3f}  sigma {pg[2]:.3f}  tilt {pg[1]:+.3f}  (fit_unbound_r280.C calibrated, vz 100-900: 2.592 3.442 4.006 / .220 / +.222)')

# 2) per theta_cm bin: a and amplitudes free; errors from the curvature of -2lnL (profile not needed at this level)
def binfit(h):
   q0 = [max(h.mean() / 1.5, .1), max(h.max() / 3, .5), max(h.max() / 3, .5), max(h.max() / 3, .5)]
   f = lambda q: nll([q[0], pg[1], pg[2], q[1], pg[4], q[2], pg[6], q[3], pg[8]], h, 3)
   r = minimize(f, q0, bounds=[(0, None)] * 4, method='L-BFGS-B'); q = r.x
   H = np.zeros((4, 4)); e = 1e-3 * np.maximum(q, 1)   # numerical Hessian of -2lnL -> cov = 2 H^-1
   for i in range(4):
      for j in range(4):
         qq = [q.copy() for _ in range(4)]
         qq[0][i] += e[i]; qq[0][j] += e[j]; qq[1][i] += e[i]; qq[1][j] -= e[j]; qq[2][i] -= e[i]; qq[2][j] += e[j]; qq[3][i] -= e[i]; qq[3][j] -= e[j]
         H[i, j] = (f(qq[0]) - f(qq[1]) - f(qq[2]) + f(qq[3])) / (4 * e[i] * e[j])
   try: cov = 2 * np.linalg.inv(H)
   except np.linalg.LinAlgError: cov = np.diag([np.inf] * 4)
   Y = q[1:] * pg[2] * np.sqrt(2 * np.pi) / BW; eY = np.sqrt(np.maximum(np.diag(cov)[1:], 0)) * pg[2] * np.sqrt(2 * np.pi) / BW
   eY = np.maximum(eY, np.sqrt(np.maximum(Y, 1)) * 0 + 1.0)
   return Y, eY, r.fun, int(h.sum())

# 3) efficiency from the Ex 3.9 sim: reco in (theta_cm bin, Ex 1-5, same cuts) / generated in (true theta_cm bin, vz window); Spyral validity
K = np.loadtxt(U + 'keys_ex39.txt'); gz = K[:, 1] + ZOFF; gx, gc = kine(Eb(gz), np.radians(K[:, 2]), K[:, 3])
gsel = (gz >= vz0) & (gz < vz1) & np.isfinite(gx); Ng = np.histogram(gc[gsel], TB)[0].astype(float)
s = uproot.open(R5 + 'p_ex39_dump.root')['r'].arrays(['ke', 'th', 'vz', 'c2'], library='np')
sx, sc = kine(Eb(s['vz']), np.radians(s['th']), s['ke'])
ss = (s['c2'] <= 13.5) & (s['th'] >= TLO) & (s['th'] < THI) & (s['ke'] >= 1.2) & np.isfinite(sx) & (s['vz'] >= vz0) & (s['vz'] < vz1) & (sx >= EXLO) & (sx < EXHI)
Nr = np.histogram(sc[ss], TB)[0].astype(float); eps = Nr / np.maximum(Ng, 1)
X = uproot.open(R5 + 'p_ex39_census.root')['c'].arrays(['thl', 'vz', 'ke', 'ex', 'tcm', 'sv'], library='np')
xb = (X['vz'] >= vz0) & (X['vz'] < vz1) & (X['thl'] >= TLO) & (X['thl'] < THI) & (X['ke'] >= 1.2) & (X['ex'] > 2.5) & (X['ex'] < 6)
spy = np.array([np.mean(X['sv'][xb & (X['tcm'] >= a) & (X['tcm'] < b)] > 0) if (xb & (X['tcm'] >= a) & (X['tcm'] < b)).sum() > 20 else np.nan for a, b in zip(TB[:-1], TB[1:])])
Leff = LPM * (vz1 - vz0)

names = ['2.6', '3.4', '4.0']; rows = []
print(f'\n theta_cm   N   -2lnL/36 | eps(gen+fit)  spy  | yields {names}  | dsdo [mb/sr] (env {ENV}, NO trigger term)')
for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
   h = np.histogram(dx[w & (dc >= a) & (dc < b)], edges)[0].astype(float)
   Y, eY, L, n = binfit(h); f = 1 / (eps[i] * spy[i] * dO[i] * Leff * ENV)
   rows.append((a, b, n, eps[i], spy[i], *Y, *eY, *(Y * f), *(eY * f)))
   print(f'  {a:.0f}-{b:.0f}  {n:4d}  {L:6.1f}   |   {eps[i]:.3f}      {spy[i]:.3f} | ' + '  '.join(f'{y:6.1f}({e:4.1f})' for y, e in zip(Y, eY)) + ' | ' + '  '.join(f'{y * f:5.2f}({e * f:4.2f})' for y, e in zip(Y, eY)))
hdr = 'thcm_lo,thcm_hi,N,eps,spy,' + ','.join(f'Y{x}' for x in names) + ',' + ','.join(f'eY{x}' for x in names) + ',' + ','.join(f'dsdo{x}' for x in names) + ',' + ','.join(f'edsdo{x}' for x in names)
with open(U + f'unb_ad_{tag}.csv', 'w') as fo:
   fo.write(f'# NOT ADOPTED, NO TRIGGER CORRECTION. unb_ad.py vz {vz0:.0f}-{vz1:.0f} theta_lab {TLO}-{THI}; global Ex {pg[4]:.3f} {pg[6]:.3f} {pg[8]:.3f} sigma {pg[2]:.3f} tilt {pg[1]:+.3f}; DPCAL {CAL}\n' + hdr + '\n')
   for r in rows: fo.write(','.join(f'{v:.5g}' for v in r) + '\n')
np.save(U + f'unb_hist_{tag}.npy', np.array([np.histogram(dx[w & (dc >= a) & (dc < b)], edges)[0] for a, b in zip(TB[:-1], TB[1:])]))
np.save(U + f'unb_pars_{tag}.npy', pg)
