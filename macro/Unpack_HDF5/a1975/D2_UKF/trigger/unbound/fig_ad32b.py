# Unbound ADs at BW offset -0.32, any bin width (TBSTEP/SUF as unb_ad3.py): all-free fits; the 2.763 + 2.980 LOW PEAK as a sum with the Fisher
# covariance (the split is unstable at fine binning), L fits incl. the shell-model mixture for the low peak; figure. NO TRIGGER TERM.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw_force.py').read().split("\nif __name__ == '__main__':")[0])
keep = np.load(U + f'unb_bw_force32{SUF}_keep.npy', allow_pickle=True).item(); LAB = 'all 6 free, offset fixed -0.32'
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
def avg(gt, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, gt[:, 0], gt[:, 1]) * w) / np.sum(w)
def tpl(L, k): gt = np.loadtxt(F % (L, k)); return np.array([avg(gt, a, b) for a, b in zip(TB[:-1], TB[1:])]), gt
mid = .5 * (TB[:-1] + TB[1:])
def ad(vz0):
   L, g, Y, H, A, qs, mode, fo, c = keep[(vz0, LAB)]; S = shapes(g[0], g[1], 1.0); fct = 1 / (dO * LPM * (900 - vz0) * ENV)
   X = Y * fct[:, None]; E = np.zeros_like(X); low = np.zeros(NB); elow = np.zeros(NB); acc = np.zeros((NB, NS))
   for i in range(NB):
      B = np.vstack([np.sqrt(np.maximum(cen - SN, 0)) * (1 + g[2] * (cen - 3)), *S]); q = qs[i]; Bm = A[i][None, :] * B; m = np.maximum(q @ Bm, 1e-9)
      I = (Bm / m) @ Bm.T; act = np.r_[True, q[1:] > 1e-6]; C = np.zeros_like(I); C[np.ix_(act, act)] = np.linalg.pinv(I[np.ix_(act, act)])
      for k in range(1, len(q)): C[k, k] = C[k, k] if act[k] else 1 / max(I[k, k], 1e-12)
      E[i] = np.sqrt(np.diag(C)[1:]) * fct[i]; low[i] = (q[1] + q[2]) * fct[i]; elow[i] = np.sqrt(C[1, 1] + C[2, 2] + 2 * C[1, 2]) * fct[i]
      acc[i] = [np.interp(ST[s][0] + g[0], cen, A[i]) for s in range(NS)]
   return X, np.maximum(E, 1e-3), low, np.maximum(elow, 1e-3), acc, L, g
def lfit(x, e, m):
   S = np.sum(x * m / e ** 2) / np.sum(m * m / e ** 2); return S, 1 / np.sqrt(np.sum(m * m / e ** 2)), np.sum(((x - S * m) / e) ** 2)
mix = {'YSOX': (.0684, .2043), 'WBP': (.0372, .2029)}
PAN = [('low peak 2.763 + 2.980', None), ('3.661', 2), ('4.231 3/2+', 3)]
fig, ax = plt.subplots(1, 3, figsize=(19, 5.4))
for vz0, col, dx in ((100, 'C0', -.25), (500, 'C3', .25)):
   X, E, low, elow, acc, L, g = ad(vz0)
   print(f'=== vz >= {vz0}, bins {TB[0]:.0f}-{TB[-1]:.0f} by {TB[1] - TB[0]:.0f} deg, fit -2lnL {L:.1f}, sigma {g[1]:.3f}')
   for j, (nm, s) in enumerate(PAN):
      if s is None:
         x, e = low, elow; u = acc[:, 1] > 0.3
         cands = {f'L1 (2.763 shape)': tpl(1, 1)[0], f'L2 (2.980 shape)': tpl(2, 2)[0]}
         for mn, (a1, a2) in mix.items(): cands[f'{mn} mix'] = a1 * tpl(1, 1)[0] + a2 * tpl(2, 2)[0]
      else:
         x, e = X[:, s], E[:, s]; u = acc[:, s] > 0.3; k = s + 1
         cands = {f'L{Lt}': tpl(Lt, k)[0] for Lt in range(4)}
      res = {n: lfit(x[u], e[u], m[u]) for n, m in cands.items()}
      print(f'   {nm:24s} ' + ' '.join(f'{v:5.2f}({w:4.2f})' if uu else '    -      ' for v, w, uu in zip(x, e, u)) + '\n      ' + ' | '.join(f'{n}: chi2 {r[2]:4.1f}/{u.sum() - 1}' + (f' scale {r[0]:.2f}({r[1]:.2f})' if 'mix' in n else '') for n, r in res.items()))
      a = ax[j]; a.errorbar(mid[u] + dx, x[u], e[u], fmt='o', color=col, capsize=2, ms=4, label=f'vz $\\geq$ {vz0}')
      if vz0 == 500:
         tt = np.linspace(TB[0], TB[-1], 200); X1, E1, low1, elow1, acc1, _, _ = ad(100)
         xr, er = (low1, elow1) if s is None else (X1[:, s], E1[:, s]); ur = acc1[:, 1] > 0.3 if s is None else acc1[:, s] > 0.3
         styles = [':', '-', '--', '-.', (0, (5, 1)), (0, (1, 1))]
         for (n, m), ls in zip(cands.items(), styles):
            if s is None: kk = {'L1 (2.763 shape)': [(1, 1, 1.)], 'L2 (2.980 shape)': [(2, 2, 1.)], 'YSOX mix': [(1, 1, .0684), (2, 2, .2043)], 'WBP mix': [(1, 1, .0372), (2, 2, .2029)]}[n]
            else: kk = [(int(n[1]), s + 1, 1.)]
            curve = sum(w * np.interp(tt, tpl(Lt, k)[1][:, 0], tpl(Lt, k)[1][:, 1]) for Lt, k, w in kk)
            Sv, eS, ch = lfit(xr[ur], er[ur], m[ur]); a.plot(tt, Sv * curve, color='gray', ls=ls, lw=1.2, label=f'{n} (fit vz$\\geq$100) $\\chi^2$ {ch:.1f}/{ur.sum() - 1}')
      a.set_title(nm, fontsize=10); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr], NO trigger term'); a.grid(alpha=.3); a.legend(fontsize=7)
for a in ax: a.set_ylim(0, None)
fig.suptitle(f'16C(d,p)17C unbound, BW offset -0.32, all states free, {TB[1] - TB[0]:.0f}-deg bins; low peak = 2.763 + 2.980 summed with covariance; NOT ADOPTED', fontsize=11)
fig.tight_layout(); fig.savefig(U + f'unb_ad32b{SUF}_2026-11-17.png', dpi=110); print('saved')
