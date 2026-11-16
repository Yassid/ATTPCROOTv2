# Per-bin spectrum fits behind unb_vzscan.csv: same data, model and fixed parameters as unb_ad2.py scan.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("K = np.loadtxt")[0])
T = '/mnt/f/a1975/sim_rebuild/trig_hdbtau/'; CAL = (184.25, 11.0, 0.18, -3, -0.05, -0.6947); TLO, THI = 89, 113
TB = np.array([22, 26, 30, 34, 38.]); edges = np.arange(1.0, 5.0001, .1); cen = .5 * (edges[1:] + edges[:-1]); BW = .1
d = np.loadtxt(T + 'gf_data.txt'); ke, th, vz, c2 = d.T.copy()
ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); dx, dc = kine(CAL[0] - CAL[1] / 1000 * (vz - 500), np.radians(th), ke); dx = dx + CAL[5]
dsel = (c2 <= 13.5) & (th >= TLO) & (th < THI) & (ke >= 1.2) & np.isfinite(dx)
def model(p, n, x=cen):
   m = p[0] * np.sqrt(np.maximum(x - 0.729, 0)) * (1 + p[1] * (x - 3))
   for k in range(n): m = m + p[3 + 2 * k] * np.exp(-.5 * ((x - p[4 + 2 * k]) / p[2]) ** 2)
   return np.maximum(m, 1e-9)
def nll(p, h, n): m = model(p, n); return 2 * np.sum(m - h + np.where(h > 0, h * np.log(np.maximum(h, 1e-12) / m), 0))
w = dsel & (vz >= 100) & (vz < 900); hg = np.histogram(dx[w], edges)[0].astype(float)
pg = minimize(nll, [hg.mean() / 1.5, .2, .22, 20, 2.6, 30, 3.45, 25, 4.0], args=(hg, 3), method='L-BFGS-B',
              bounds=[(0, None), (-.6, 2), (.08, .9), (0, None), (2.2, 3.0), (0, None), (3.0, 3.8), (0, None), (3.7, 4.6)]).x
xf = np.linspace(1, 5, 400); cols = ['tab:green', 'tab:orange', 'tab:purple']; names = ['2.59', '3.45', '4.01']
fig = plt.figure(figsize=(20, 11)); gs = fig.add_gridspec(4, 4, height_ratios=[3, 1, 3, 1], hspace=.08, wspace=.18)
for r, vz0 in enumerate((100, 500)):
   for j, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
      h = np.histogram(dx[dsel & (vz >= vz0) & (vz < 900) & (dc >= a) & (dc < b)], edges)[0].astype(float)
      f = lambda q: nll([q[0], pg[1], pg[2], q[1], pg[4], q[2], pg[6], q[3], pg[8]], h, 3)
      res = minimize(f, [max(h.mean() / 1.5, .1)] + [max(h.max() / 3, .5)] * 3, bounds=[(0, None)] * 4, method='L-BFGS-B'); q = res.x
      p = [q[0], pg[1], pg[2], q[1], pg[4], q[2], pg[6], q[3], pg[8]]
      ax = fig.add_subplot(gs[2 * r, j]); rx = fig.add_subplot(gs[2 * r + 1, j], sharex=ax)
      ax.errorbar(cen, h, np.sqrt(np.maximum(h, 1)), fmt='o', ms=3, color='k', lw=1)
      ax.plot(xf, model(p, 3, xf), 'r', lw=1.5, label=f'total, -2lnL {res.fun:.1f}/36')
      ax.plot(xf, model([q[0], pg[1], pg[2]], 0, xf), 'b--', lw=1, label='1n phase space x tilt')
      for k in range(3):
         Y = q[1 + k] * pg[2] * np.sqrt(2 * np.pi) / BW
         ax.plot(xf, q[1 + k] * np.exp(-.5 * ((xf - pg[4 + 2 * k]) / pg[2]) ** 2), color=cols[k], lw=1.3, label=f'{names[k]}: {Y:.0f} counts')
      ax.set_title(f'vz $\\geq$ {vz0} mm, $\\theta_{{cm}}$ {a:.0f}-{b:.0f}  (N {int(h.sum())})' + ('   [NOT USED: eps 3.4/4.0 .29/.11]' if a == 34 else ''), fontsize=9)
      ax.legend(fontsize=7, loc='upper left'); ax.grid(alpha=.3); ax.set_ylim(0, None); plt.setp(ax.get_xticklabels(), visible=False)
      m = model(p, 3); pull = (h - m) / np.sqrt(np.maximum(m, 1)); rx.bar(cen, pull, width=.09, color='gray'); rx.axhline(0, color='k', lw=.6)
      rx.set_ylim(-3.5, 3.5); rx.grid(alpha=.3); rx.set_xlabel(r'E$_x$ [MeV]' if r == 1 else ''); rx.set_ylabel('pull', fontsize=8)
      if j == 0: ax.set_ylabel('counts / 100 keV')
fig.suptitle(f'16C(d,p) unbound, per-bin fits: positions {pg[4]:.3f} / {pg[6]:.3f} / {pg[8]:.3f} MeV, $\\sigma$ {pg[2]:.3f}, tilt {pg[1]:+.3f} fixed from vz 100-900; '
             'free: phase-space scale + 3 amplitudes. $\\theta_{lab}$ 89-113, calibrated data.', fontsize=11)
fig.savefig('unb_fits_2026-11-16.png', dpi=105, bbox_inches='tight'); print('saved')
