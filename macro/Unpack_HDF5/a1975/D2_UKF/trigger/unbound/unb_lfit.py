# L assignment test for the three unbound states on unb_vzscan.csv bins 22-34 (eps >= 0.5): one normalisation per L template
# (C17_dp_fits/fresco/outputs/dp17C_L<L>_state<k>.dat, ADWA at E_d 23.55, post coupling fix 2026-10-03), sin-weighted bin average.
import numpy as np, sys
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
_L = [l.strip().split(',') for l in open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_vzscan.csv') if not l.startswith('#')]
D = {k: np.array([float(r[i]) for r in _L[1:]]) for i, k in enumerate(_L[0])}
ST = {'2.6': 1, '3.4': 3, '4.0': 4}   # template state (2.763, 3.661, 4.231)
def avg(g, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w)
for vz0 in (100, 300, 500):
   print(f'=== vz >= {vz0}')
   for s, k in ST.items():
      sel = (D['vzLo'] == vz0) & (D['thcm_hi'] <= 34); r = {kk: vv[sel] for kk, vv in D.items()}
      x, e, ep = r[f'dsdo{s}'], r[f'edsdo{s}'], r[f'eps{s}']; u = ep >= 0.5; x, e = x[u], np.maximum(e[u], 1e-3)
      line = f'  {s}: ' + ' '.join(f'{v:.2f}({w:.2f})' for v, w in zip(x, e)) + '  |'
      for L in range(4):
         g = np.loadtxt(F % (L, k)); m = np.array([avg(g, a, b) for a, b in zip(r['thcm_lo'][u], r['thcm_hi'][u])])
         S = np.sum(x * m / e ** 2) / np.sum(m * m / e ** 2); c = np.sum(((x - S * m) / e) ** 2)
         line += f' L{L} S {S:.3f} chi2 {c:4.1f}/{u.sum() - 1} |'
      print(line)
