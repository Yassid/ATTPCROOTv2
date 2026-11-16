# Unbound ADs at BW offset -0.32, all-free fits, vz >= 100 and >= 500, with L=0..3 templates normalised to vz >= 500. NO trigger term.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
R = [l.strip().split(',') for l in open('unb_ad32.csv') if not l.startswith('#')]; hdr = R[0]; R = R[1:]
mid = np.array([24, 28, 32, 36.]); TB = np.array([22, 26, 30, 34, 38.])
def avg(gt, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, gt[:, 0], gt[:, 1]) * w) / np.sum(w)
ST = [('2.763 1/2-', 1, 'shell model 1/2-: L=1'), ('2.980 3/2+', 2, 'shell model 3/2+: L=2'), ('3.661 (3/2+)', 3, 'not in shell model'), ('4.231 3/2+', 4, 'shell model 3/2+: L=2')]
fig, ax = plt.subplots(1, 4, figsize=(20, 5))
for j, (nm, k, note) in enumerate(ST):
   a = ax[j]; ref = None
   for vz0, c, dx in (('100', 'C0', -.4), ('500', 'C3', .4)):
      r = [x for x in R if x[0] == vz0 and x[1] == 'FREE' and x[2] == nm][0]
      x = np.array([float(v) for v in r[3:7]]); e = np.array([float(v) for v in r[7:11]]); chi = [float(v) for v in r[11:15]]
      u = ~((nm.startswith('4.231')) & (np.arange(4) == 3))
      a.errorbar(mid[u] + dx, x[u], e[u], fmt='o', color=c, capsize=3, label=f'vz $\\geq$ {vz0}: $\\chi^2$ L0-3 ' + '/'.join(f'{v:.1f}' for v in chi))
      if vz0 == '500': ref = (x[u], e[u], u)
   tt = np.linspace(20, 40, 200)
   for Lt, ls in ((0, ':'), (1, '-'), (2, '--'), (3, '-.')):
      gt = np.loadtxt(F % (Lt, k)); m = np.array([avg(gt, p, q) for p, q in zip(TB[:-1], TB[1:])])[ref[2]]
      S = np.sum(ref[0] * m / ref[1] ** 2) / np.sum(m * m / ref[1] ** 2); a.plot(tt, S * np.interp(tt, gt[:, 0], gt[:, 1]), color='gray', ls=ls, lw=1.2, label=f'L={Lt} (norm. to vz$\\geq$500)')
   a.set_title(f'{nm}  ({note})', fontsize=10); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr], NO trigger term'); a.set_xlim(20, 39)
   top = max(np.max(np.array([float(v) for v in xx[3:7]]) + np.array([float(v) for v in xx[7:11]])) for xx in R if xx[2] == nm and xx[1] == 'FREE'); a.set_ylim(0, 1.6 * top)
   a.grid(alpha=.3); a.legend(fontsize=7)
fig.suptitle('16C(d,p)17C unbound, BW fits with offset -0.32 MeV (bound-blend centroid), all states free per bin; acceptance in the model; env 0.84; NOT ADOPTED', fontsize=11)
fig.tight_layout(); fig.savefig('unb_ad32_2026-11-17.png', dpi=110); print('saved')
