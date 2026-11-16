# 16C(d,p)17C unbound states with BREIT-WIGNER lineshapes (2026-11-16, Yassid: "add BW and consider the states in the repository"). NOT ADOPTED,
# NO TRIGGER TERM. States = the C17_dp_fits table (Movilla & Ayyad 2026-02 total fit, the "shell-model priors" of C16_pd_AngBins.C):
#   E0 2.763 1/2- (L1) G 50 keV, 2.980 3/2+ (L2) 200, 3.661 3/2+(5/2+) (L2) 398, 4.231 3/2+ (L2) 500, 4.841 1/2+ (L0) 300, 5.91 3/2+ (L2) 1500.
# Lineshape: Gamma(E) = G0 P_L(E-Sn)/P_L(E0-Sn) (neutron hard-sphere penetrability, R = 1.25 (16^1/3 + 1) fm), BW ~ Gamma(E)/((E-E0)^2 + Gamma(E)^2/4),
# unit area, folded with a gaussian resolution sigma. Model per theta_cm bin = A_bin(Ex) x [1n phase space x tilt + sum_k Y_k shape_k]
# (A = acceptance as in unb_ad3.py), so Y_k are acceptance-corrected counts of state k (full lineshape area).
# Shared: offset delta (E0 -> E0 + delta; the calibration puts the bound BLEND at 0, so ~ -0.2 expected), sigma, tilt. Per bin: PS scale, 6 Y.
#   ~/Spyral/venv/bin/python -I unb_bw.py [vzLo ...]        (default 100 500)     options via env: GSCALE (all widths x), GFREE=1 (shared width scale free)
import os, sys, numpy as np
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_ad3.py').read().split("if __name__ == '__main__':")[0])
SN = 0.729; MN = 939.565; MC = 16.0147013 * 931.494 - 6 * 0.511; MU = MN * MC / (MN + MC); HBC = 197.327; RAD = 1.25 * (16 ** (1 / 3) + 1)
ST = [(2.763, 1, .050, '2.763 1/2-'), (2.980, 2, .200, '2.980 3/2+'), (3.661, 2, .398, '3.661 (3/2+)'), (4.231, 2, .500, '4.231 3/2+'), (4.841, 0, .300, '4.841 1/2+'), (5.910, 2, 1.500, '5.91 3/2+')]
NS = len(ST); GSCALE = float(os.environ.get('GSCALE', 1.0)); GFREE = os.environ.get('GFREE', '0') == '1'
def pen(L, e):
   k = np.sqrt(2 * MU * np.maximum(e, 1e-9)) / HBC; r = k * RAD; r2 = r * r
   v = [np.ones_like(r), r2 / (1 + r2), r2 ** 2 / (9 + 3 * r2 + r2 ** 2), r2 ** 3 / (225 + 45 * r2 + 6 * r2 ** 2 + r2 ** 3)][L]
   return r * v
EG = np.arange(SN + 1e-4, 10.0, 0.005)
def bw_shape(E0, L, G0):
   G = G0 * pen(L, EG - SN) / pen(L, E0 - SN); s = G / ((EG - E0) ** 2 + G * G / 4); return s / (s.sum() * 0.005)
KX = np.arange(-1.5, 1.5001, 0.005)
def shapes(delta, sig, gsc):
   out = []
   g = np.exp(-.5 * (KX / sig) ** 2); g /= g.sum()
   for E0, L, G0, _ in ST:
      s = np.convolve(bw_shape(E0 + delta, L, G0 * GSCALE * gsc), g, mode='same'); out.append(np.interp(cen, EG, s) * BW)   # counts per 100 keV bin per unit area
   return np.array(out)

def _binfit(h, Ai, basis, q0=None):
   """Poisson fit linear in q (basis rows: phase space, then the NS shapes), q >= 0, analytic gradient. Returns q, -2lnL."""
   Bm = Ai[None, :] * basis
   def f(q):
      m = np.maximum(q @ Bm, 1e-9); val = 2 * np.sum(m - h + np.where(h > 0, h * np.log(np.maximum(h, 1e-12) / m), 0))
      return val, 2 * Bm @ (1 - h / m)
   if q0 is None: q0 = np.r_[max(h.sum() / 40, .1), [max(h.sum() / 20, .5)] * (len(basis) - 1)]
   r = minimize(f, q0, jac=True, bounds=[(0, None)] * len(basis), method='L-BFGS-B', options={'maxiter': 5000})
   return r.x, r.fun

def _basis(x, S):
   return np.vstack([np.sqrt(np.maximum(cen - SN, 0)) * (1 + x[2] * (cen - 3)), S])

def fit_bw(H, A):
   """outer: shared (offset, sigma, tilt[, Gamma scale]) by Nelder-Mead; inner: per-bin linear Poisson fits. Returns an object with .x/.fun
   laid out as before: [shared..., per bin (PS, Y1..YNS)]."""
   ng = 4 if GFREE else 3
   lo_b = [-0.6, 0.05, -0.6] + ([0.3] if GFREE else []); hi_b = [0.3, 0.5, 2.0] + ([3.0] if GFREE else [])
   cache = {}
   def outer(g):
      g = np.clip(g, lo_b, hi_b); S = shapes(g[0], g[1], g[3] if GFREE else 1.0); B = _basis(g, S); tot = 0; qs = []
      for i in range(NB): q, L = _binfit(H[i], A[i], B); tot += L; qs.append(q)
      cache[tuple(g)] = qs; return tot
   best = None
   for d0 in (-0.35, -0.2, -0.05):
      g0 = [d0, 0.15, 0.2] + ([1.0] if GFREE else [])
      r = minimize(outer, g0, method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-3, 'maxiter': 600})
      if best is None or r.fun < best.fun: best = r
   g = np.clip(best.x, lo_b, hi_b); outer(g); qs = cache[tuple(g)]
   class R: pass
   R.x = np.r_[g, np.concatenate(qs)]; R.fun = best.fun
   return R, ng, lambda x: shapes(x[0], x[1], x[3] if GFREE else 1.0)

def yerr(H, A, x, ng, S):
   """per-bin amplitude errors from the Fisher matrix of the linear Poisson model at the best fit (shared parameters held)."""
   E = np.zeros((NB, NS)); B = _basis(x, S)
   for i in range(NB):
      o = ng + (NS + 1) * i; q = x[o:o + 1 + NS]; Bm = A[i][None, :] * B; m = np.maximum(q @ Bm, 1e-9)
      I = (Bm / m) @ Bm.T
      act = np.r_[True, q[1:] > 1e-6]                 # amplitudes at the bound: error from the 1-parameter curvature
      C = np.zeros_like(I); C[np.ix_(act, act)] = np.linalg.pinv(I[np.ix_(act, act)])
      for k in range(NS): E[i, k] = np.sqrt(C[k + 1, k + 1]) if act[k + 1] else 1 / np.sqrt(max(I[k + 1, k + 1], 1e-12))
   return E

if __name__ == '__main__':
   VZS = [int(v) for v in sys.argv[1:]] or [100, 500]
   rows = []; keep = {}
   tag = os.environ.get('TAG', 'bw')
   for vz0 in VZS:
      vz1 = 900; ww = dsel & (vz >= vz0) & (vz < vz1); Leff = LPM * (vz1 - vz0)
      H = np.array([np.histogram(dx[ww & (dc >= a) & (dc < b)], edges)[0] for a, b in zip(TB[:-1], TB[1:])]).astype(float)
      A = acceptance(vz0, vz1); r, ng, parts = fit_bw(H, A); x = r.x; S = parts(x); E = yerr(H, A, x, ng, S)
      npar = len(x); print(f'=== BW vz >= {vz0}: -2lnL {r.fun:.1f}/{H.size - npar}  offset {x[0]:+.3f}  sigma_res {x[1]:.3f}  tilt {x[2]:+.3f}' + (f'  Gamma scale {x[3]:.2f}' if GFREE else '') + f'  (GSCALE {GSCALE})')
      print('   E0 in this frame: ' + '  '.join(f'{e + x[0]:.3f}' for e, *_ in ST))
      print('   bin   ' + ' '.join(f'{s[3]:>17s}' for s in ST) + '   [acceptance-corrected yield (err) -> dsdo mb/sr, env .84, NO trigger]')
      for i in range(NB):
         o = ng + (NS + 1) * i; Y = x[o + 1:o + 1 + NS]; fct = 1 / (dO[i] * Leff * ENV)
         m = A[i] * (x[o] * np.sqrt(np.maximum(cen - SN, 0)) * (1 + x[2] * (cen - 3)) + Y @ S)
         print(f'   {TB[i]:.0f}-{TB[i+1]:.0f} ' + ' '.join(f'{y:5.0f}({e:3.0f}){y * fct:5.2f}' for y, e in zip(Y, E[i])) + f'   -2lnL {pnll(m, H[i]):.1f}/{len(cen)}')
         rows.append((vz0, TB[i], TB[i + 1], int(H[i].sum()), *Y, *E[i], *(Y * fct), *(E[i] * fct)))
      keep[vz0] = dict(H=H, A=A, x=x, ng=ng)
      print()
   names = [f'{s[0]:.3f}' for s in ST]
   hdr = 'vzLo,thcm_lo,thcm_hi,N,' + ','.join('Y' + n for n in names) + ',' + ','.join('eY' + n for n in names) + ',' + ','.join('dsdo' + n for n in names) + ',' + ','.join('edsdo' + n for n in names)
   with open(U + f'unb_{tag}.csv', 'w') as fo:
      fo.write(f'# NOT ADOPTED, NO TRIGGER CORRECTION. unb_bw.py BW states {[s[3] for s in ST]}, GSCALE {GSCALE} GFREE {GFREE}, Ex {EXLO}-{EXHI}, theta_lab {TLO}-{THI}, env {ENV}\n' + hdr + '\n')
      for rr in rows: fo.write(','.join(f'{v:.5g}' for v in rr) + '\n')
   np.save(U + f'unb_{tag}_keep.npy', keep, allow_pickle=True)
