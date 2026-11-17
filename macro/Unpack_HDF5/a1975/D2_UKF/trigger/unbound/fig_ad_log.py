# Log-scale unbound ADs (working result: THI=137 KEMIN=0.8 TBLO=10 TBSTEP=4 SUF=_137ke08, offset -0.32, all states free). NO TRIGGER TERM.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
src = open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/fig_ad32b.py').read(); exec(src.split("mix = {'YSOX'")[0])   # keep, ad(), tpl(), lfit(), mid
PAN = [('2.763 1/2-', 0, 1, 'SM: L=1'), ('2.980 3/2+', 1, 2, 'SM: L=2'), ('low peak 2.763 + 2.980', None, None, 'sum, covariance'), ('3.661', 2, 3, 'no SM state'), ('4.231 3/2+', 3, 4, 'SM: L=2')]
D = {vz0: ad(vz0) for vz0 in (100, 500)}
fig, ax = plt.subplots(1, 5, figsize=(25, 5.4)); tt = np.linspace(TB[0], TB[-1], 300); styles = {0: ':', 1: '-', 2: '--', 3: '-.'}
for j, (nm, s, k, note) in enumerate(PAN):
   a = ax[j]
   for vz0, col, dx in ((100, 'C0', -.3), (500, 'C3', .3)):
      X, E, low, elow, acc, L, g = D[vz0]
      if s is None: x, e, u = low, elow, acc[:, 1] > 0.3
      else: x, e, u = X[:, s], E[:, s], acc[:, s] > 0.3
      pos = u & (x > 0); lo = np.minimum(e, 0.999 * x)          # log axis: asymmetric lower bar, zero points drawn as upper limits
      a.errorbar(mid[pos] + dx, x[pos], [lo[pos], e[pos]], fmt='o', color=col, capsize=2, ms=4, label=f'vz $\\geq$ {vz0}')
      z = u & (x <= 0)
      if z.any(): a.errorbar(mid[z] + dx, e[z], yerr=0.3 * e[z], uplims=True, fmt='none', color=col)
   X, E, low, elow, acc, L, g = D[100]
   x, e, u = (low, elow, acc[:, 1] > 0.3) if s is None else (X[:, s], E[:, s], acc[:, s] > 0.3)
   if s is None:
      cand = [('L=1 (2.763 shape)', [(1, 1, 1.)]), ('L=2 (2.980 shape)', [(2, 2, 1.)]), ('YSOX mix', [(1, 1, .0684), (2, 2, .2043)]), ('WBP mix', [(1, 1, .0372), (2, 2, .2029)])]
      lss = [':', '-', '--', '-.']
   else:
      cand = [(f'L={Lt}', [(Lt, k, 1.)]) for Lt in range(4)]; lss = [styles[Lt] for Lt in range(4)]
   for (n, kk), ls in zip(cand, lss):
      m = sum(w * tpl(Lt, kq)[0] for Lt, kq, w in kk); Sv, eS, ch = lfit(x[u], np.maximum(e[u], 1e-3), m[u])
      curve = sum(w * np.interp(tt, tpl(Lt, kq)[1][:, 0], tpl(Lt, kq)[1][:, 1]) for Lt, kq, w in kk)
      a.plot(tt, Sv * curve, color='gray', ls=ls, lw=1.3, label=f'{n}: $\\chi^2$ {ch:.1f}/{u.sum() - 1} (vz$\\geq$100)')
   a.set_yscale('log'); a.set_xlim(TB[0], TB[-1]); a.set_title(f'{nm}  ({note})', fontsize=10); a.set_xlabel(r'$\theta_{cm}$ [deg]')
   a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr], NO trigger term'); a.grid(alpha=.3, which='both'); a.legend(fontsize=7, loc='lower left')
   vals = np.concatenate([D[v][0][:, s] if s is not None else D[v][2] for v in (100, 500)]); vals = vals[vals > 0]
   a.set_ylim(max(vals.min() / 5, 0.03), vals.max() * 4)
fig.suptitle('16C(d,p)17C unbound, BW fits (offset -0.32, all states free), theta_lab 89-137, KE >= 0.8 MeV, 4-deg bins; curves normalised to vz >= 100; NOT ADOPTED', fontsize=11)
fig.tight_layout(); fig.savefig(U + f'unb_ad_log{SUF}_2026-11-17.png', dpi=110); print('saved')
