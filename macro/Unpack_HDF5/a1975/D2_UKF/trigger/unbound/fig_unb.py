# Figure: unbound ADs per state for vz >= 100/300/500 (no trigger term) + L templates normalised to vz >= 500; spectrum fit inset.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
_L = [l.strip().split(',') for l in open('unb_vzscan.csv') if not l.startswith('#')]; D = {k: np.array([float(r[i]) for r in _L[1:]]) for i, k in enumerate(_L[0])}
ST = [('2.6', 1, '2.59 (lit. 2.763 1/2$^-$)'), ('3.4', 3, '3.44 (lit. 3.661)'), ('4.0', 4, '4.01 (lit. 4.231 3/2$^+$)')]
def avg(g, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w)
fig, ax = plt.subplots(1, 4, figsize=(19, 4.8))
for j, (s, k, lab) in enumerate(ST):
   a = ax[j]
   for vz0, c, dx in ((100, 'C0', -.5), (300, 'C2', 0), (500, 'C3', .5)):
      u = (D['vzLo'] == vz0) & (D['thcm_hi'] <= 34) & (D[f'eps{s}'] >= .5); mid = .5 * (D['thcm_lo'][u] + D['thcm_hi'][u])
      a.errorbar(mid + dx, D[f'dsdo{s}'][u], D[f'edsdo{s}'][u], fmt='o', color=c, capsize=3, label=f'vz $\\geq$ {vz0} mm')
      if vz0 == 500: x, e, lo, hi = D[f'dsdo{s}'][u], D[f'edsdo{s}'][u], D['thcm_lo'][u], D['thcm_hi'][u]
   tt = np.linspace(15, 40, 200)
   for L, ls in ((0, ':'), (1, '-'), (2, '--'), (3, '-.')):
      g = np.loadtxt(F % (L, k)); m = np.array([avg(g, p, q) for p, q in zip(lo, hi)]); S = np.sum(x * m / e ** 2) / np.sum(m * m / e ** 2)
      c2 = np.sum(((x - S * m) / e) ** 2); a.plot(tt, S * np.interp(tt, g[:, 0], g[:, 1]), color='gray', ls=ls, lw=1.2, label=f'L={L} $\\chi^2$ {c2:.1f}/2 (vz$\\geq$500)')
   a.set_title(f'E$_x$ {lab}', fontsize=10); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr], NO trigger term'); a.set_xlim(18, 38); a.set_ylim(0, 1.8 * np.max(D[f'dsdo{s}'][(D['thcm_hi'] <= 34)] + D[f'edsdo{s}'][(D['thcm_hi'] <= 34)])); a.grid(alpha=.3); a.legend(fontsize=7)
h = np.load('unb_hist_main.npy').sum(0); p = np.load('unb_pars_main.npy'); e = np.arange(1, 5.0001, .1); c = .5 * (e[1:] + e[:-1])
a = ax[3]; a.step(c, h, where='mid', color='k', label='data, $\\theta_{lab}$ 89-113, vz 100-900')
bg = p[0] * np.sqrt(np.maximum(c - .729, 0)) * (1 + p[1] * (c - 3)); a.plot(c, bg, 'b--', label='1n phase space x tilt')
tot = bg.copy()
for k in range(3): gk = p[3 + 2 * k] * np.exp(-.5 * ((c - p[4 + 2 * k]) / p[2]) ** 2); tot += gk; a.plot(c, gk, color='g', lw=1)
a.plot(c, tot, 'r', label=f'fit: {p[4]:.2f} / {p[6]:.2f} / {p[8]:.2f}, $\\sigma$ {p[2]:.2f}')
a.set_xlabel(r'E$_x$ [MeV] (calibrated; bound blend at 0)'); a.set_ylabel('counts / 100 keV'); a.legend(fontsize=7); a.grid(alpha=.3); a.set_title('global spectrum fit (positions fixed per bin)', fontsize=10)
fig.tight_layout(); fig.savefig('unb_ad_2026-11-16.png', dpi=115); print('saved')
