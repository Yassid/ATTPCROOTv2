# Check (Yassid 2026-11-17: "do the check with the 2.980 forced in"): can the low unbound peak be the 2.980 3/2+ the shell models predict?
# Variants vs the free BW fit of unb_bw.py: A drop 2.763 (offset free); B drop 2.763, offset fixed -0.21; C/D tie Y(2.980) = R_bin Y(2.763) with
# R_bin = C2S(3/2+)/C2S(1/2-) x [L2 state2 / L1 state1 FRESCO dsigma, bin average] for YSOX (.2043/.0684) / WBP (.2029/.0372).
import os, numpy as np
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw.py').read().split("\nif __name__ == '__main__':")[0])
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
def avg(g, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w)
g1, g2 = np.loadtxt(F % (1, 1)), np.loadtxt(F % (2, 2))
XR = np.array([avg(g2, a, b) / avg(g1, a, b) for a, b in zip(TB[:-1], TB[1:])])   # sigma(2.980, L2) / sigma(2.763, L1) per bin, unit C2S

def run(vz0, mode, fixoff=None, c2s=None):
   ww = dsel & (vz >= vz0) & (vz < 900); H = np.array([np.histogram(dx[ww & (dc >= a) & (dc < b)], edges)[0] for a, b in zip(TB[:-1], TB[1:])]).astype(float)
   A = acceptance(vz0, 900)
   def basis_bin(g, i):
      S = shapes(g[0], g[1], 1.0); ps = np.sqrt(np.maximum(cen - SN, 0)) * (1 + g[2] * (cen - 3))
      if mode == 'free': rows = [ps, *S]
      elif mode == 'drop': rows = [ps, *S[1:]]                       # 2.980 3.661 4.231 4.841 5.91
      else: R = c2s * XR[i]; rows = [ps, S[0] + R * S[1], *S[2:]]    # combined low state, then 3.661 ...
      return np.vstack(rows)
   lo_b, hi_b = [-0.7, 0.05, -0.6], [0.3, 0.5, 2.0]; cache = {}
   def outer(gfree):
      g = np.clip(np.r_[fixoff, gfree] if fixoff is not None else gfree, lo_b, hi_b); tot = 0; qs = []
      for i in range(NB): q, L = _binfit(H[i], A[i], basis_bin(g, i)); tot += L; qs.append(q)
      cache[tuple(g)] = qs; return tot
   best = None
   starts = [(d0, t0) for d0 in ((-0.45, -0.3, -0.15, 0.0) if fixoff is None else (None,)) for t0 in (0.25, -0.3)]   # nonzero tilt starts: NM step on a 0 start is 0.00025
   for d0, t0 in starts:
      g0 = [0.18, t0] if fixoff is not None else [d0, 0.18, t0]
      r = minimize(outer, g0, method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 800})
      if best is None or r.fun < best.fun: best = r
   g = np.clip(np.r_[fixoff, best.x] if fixoff is not None else best.x, lo_b, hi_b); outer(best.x); qs = cache[tuple(g)]
   Y = []
   for i, q in enumerate(qs):
      if mode == 'free': Y.append(q[1:7])
      elif mode == 'drop': Y.append(np.r_[0, q[1:6]])
      else: R = c2s * XR[i]; Y.append(np.r_[q[1], R * q[1], q[2:6]])
   npar = len(best.x) + sum(len(q) for q in qs)
   return best.fun, H.size - npar, g, np.array(Y), H, A, qs

if __name__ == '__main__':
   print('sigma(2.980,L2)/sigma(2.763,L1) per bin, unit C2S:', np.round(XR, 2), '-> YSOX ratio x', round(.2043 / .0684, 2), ', WBP x', round(.2029 / .0372, 2))
   V = [('baseline: all 6 free, offset free', 'free', None, None), ('A: 2.763 dropped, offset free', 'drop', None, None),
        ('B: 2.763 dropped, offset fixed -0.21', 'drop', -0.21, None), ('C: 2.980 tied to 2.763 at YSOX ratio', 'tie', None, .2043 / .0684),
        ('D: 2.980 tied to 2.763 at WBP ratio', 'tie', None, .2029 / .0372)]
   out = {}
   for vz0 in (100, 500):
      print(f'\n=== vz >= {vz0}    yields summed over the 4 bins (acceptance-corrected counts):  2.763  2.980  3.661  4.231  4.841')
      for lab, mode, fo, c in V:
         L, ndf, g, Y, H, A, qs = run(vz0, mode, fo, c); Ys = Y.sum(0)
         base = out.get((vz0, 'base'))
         if mode == 'free' and fo is None: out[(vz0, 'base')] = L
         dl = '' if base is None else f'  d(-2lnL) vs baseline {L - base:+6.1f}'
         print(f'  {lab:40s} -2lnL {L:6.1f}/{ndf}{dl}  offset {g[0]:+.3f} sigma {g[1]:.3f} | ' + '  '.join(f'{v:6.0f}' for v in Ys[:5]))
         out[(vz0, lab)] = (L, g, Y, H, A, qs, mode, fo, c)
   np.save(U + 'unb_bw_force_keep.npy', {k: v for k, v in out.items() if k[1] != 'base'}, allow_pickle=True)
