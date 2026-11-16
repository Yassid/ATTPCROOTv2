# 16C(d,p)17C unbound states, step 3 (2026-11-16, after looking at the per-bin fits of unb_ad2.py). NOT ADOPTED, NO TRIGGER TERM.
# Fixes: (1) the model of each theta_cm bin is multiplied by that bin's acceptance A(Ex) (theta_lab window + KE >= 1.2 by forward kinematics,
# x reco efficiency(true theta_lab, vz) of the 3.9 sim, x Spyral validity), so the theta_lab-89 edge is in the model and the amplitudes are
# acceptance-corrected; (2) peak positions: diagnostic fit with free positions per bin, production fit with one common shift per bin;
# (3) fit range Ex 1.3-5.0. The bins are fitted SIMULTANEOUSLY: shared sigma, tilt, three positions; per bin phase-space scale, 3 amplitudes, shift.
#   ~/Spyral/venv/bin/python -I unb_ad3.py [vzLo ...]      (default: 100 300 500 for the table; figure at 100 and 500)
import sys, numpy as np, uproot
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_ad2.py').read().split("if MODE == 'val':")[0].replace("MODE = sys.argv[1] if len(sys.argv) > 1 else 'scan'", "MODE = 'lib'"))
U = '/mnt/f/a1975/sim_rebuild/unb_hdbtau/'
TB = np.array([22, 26, 30, 34, 38.]); NB = len(TB) - 1; dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))
EXLO, EXHI, BW = 1.3, 5.0, 0.1; edges = np.arange(EXLO, EXHI + 1e-9, BW); cen = .5 * (edges[1:] + edges[:-1])
d = np.loadtxt(T + 'gf_data.txt'); ke, th, vz, c2 = d.T.copy()
ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); dx, dc = kine(CAL[0] - CAL[1] / 1000 * (vz - 500), np.radians(th), ke); dx = dx + CAL[5]
dsel = (c2 <= 13.5) & (th >= TLO) & (th < THI) & (ke >= 1.2) & np.isfinite(dx)

def acceptance(vz0, vz1):   # A[bin, Ex-bin centre]
   return np.array([eps_state(x, TB, vz0, vz1) for x in cen]).T

def shape(p_bg, tilt, amps, pos, sig, x=cen):
   m = p_bg * np.sqrt(np.maximum(x - 0.729, 0)) * (1 + tilt * (x - 3))
   for a_, m_ in zip(amps, pos): m = m + a_ * np.exp(-.5 * ((x - m_) / sig) ** 2)
   return m
def pnll(m, h): m = np.maximum(m, 1e-9); return 2 * np.sum(m - h + np.where(h > 0, h * np.log(np.maximum(h, 1e-12) / m), 0))

def fit_all(H, A, mode):
   """mode 'shift': shared pos(3), sig, tilt + per bin [bg, a1, a2, a3, shift]; 'free': shared sig, tilt + per bin [bg, a1..3, p1..3]."""
   if mode == 'shift':
      x0 = [2.6, 3.45, 4.0, .22, .2] + sum([[5, 10, 15, 10, 0.] for _ in range(NB)], [])
      bnd = [(2.2, 3.0), (3.0, 3.8), (3.7, 4.6), (.08, .6), (-.6, 2)] + sum([[(0, None)] * 4 + [(-.3, .3)] for _ in range(NB)], [])
      def unpack(x, i): o = 5 + 5 * i; return x[o], x[4], x[o + 1:o + 4], np.array(x[0:3]) + x[o + 4], x[3]
   else:
      x0 = [.22, .2] + sum([[5, 10, 15, 10, 2.6, 3.45, 4.0] for _ in range(NB)], [])
      bnd = [(.08, .6), (-.6, 2)] + sum([[(0, None)] * 4 + [(2.2, 3.0), (3.0, 3.85), (3.7, 4.6)] for _ in range(NB)], [])
      def unpack(x, i): o = 2 + 7 * i; return x[o], x[1], x[o + 1:o + 4], x[o + 4:o + 7], x[0]
   f = lambda x: sum(pnll(A[i] * shape(*unpack(x, i)), H[i]) for i in range(NB))
   best = None
   for trial in range(3):   # a few starts
      xs = np.array(x0, float); xs[5 if mode == 'shift' else 2::5 if mode == 'shift' else 7] *= (1 + .5 * trial)
      r = minimize(f, xs, bounds=bnd, method='L-BFGS-B', options={'maxiter': 20000, 'maxfun': 200000})
      if best is None or r.fun < best.fun: best = r
   return best, unpack, f

def amp_errors(H, A, x, unpack, mode):
   """1-sigma per amplitude from the -2lnL profile of that bin, shared parameters held at the best fit."""
   E = np.zeros((NB, 3))
   for i in range(NB):
      bg, tilt, amps, pos, sig = unpack(x, i)
      def fb(q): return pnll(A[i] * shape(q[0], tilt, q[1:4], pos, sig), H[i])
      q0 = np.r_[bg, amps]; L0 = fb(q0)
      for k in range(3):
         def prof(v):
            g = lambda qq: fb(np.insert(qq, k + 1, v)); return minimize(g, np.delete(q0, k + 1), bounds=[(0, None)] * 3, method='L-BFGS-B').fun
         lo, hi = q0[k + 1], q0[k + 1] + 1.0
         while prof(hi) - L0 < 1 and hi < 1e5: hi *= 1.5
         for _ in range(30): mid = .5 * (lo + hi); lo, hi = (mid, hi) if prof(mid) - L0 < 1 else (lo, mid)
         E[i, k] = hi - q0[k + 1]
   return E

if __name__ == '__main__':
   VZS = [int(v) for v in sys.argv[1:]] or [100, 300, 500]
   rows = []; keep = {}
   for vz0 in VZS:
      vz1 = 900; ww = dsel & (vz >= vz0) & (vz < vz1); Leff = LPM * (vz1 - vz0)
      H = np.array([np.histogram(dx[ww & (dc >= a) & (dc < b)], edges)[0] for a, b in zip(TB[:-1], TB[1:])]).astype(float)
      A = acceptance(vz0, vz1)
      rf, uf, ff = fit_all(H, A, 'free')
      rs, us, fs = fit_all(H, A, 'shift')
      print(f'=== vz >= {vz0}: counts per bin {H.sum(1).astype(int)}  | FREE positions -2lnL {rf.fun:.1f}/{H.size - len(rf.x)}  SHIFT model {rs.fun:.1f}/{H.size - len(rs.x)}  (diff {rs.fun - rf.fun:.1f} for {len(rf.x) - len(rs.x)} fewer pars)')
      for i in range(NB):
         _, _, af, pf, sf = uf(rf.x, i); _, _, as_, ps, ss = us(rs.x, i)
         print(f'   {TB[i]:.0f}-{TB[i+1]:.0f}  free pos ' + ' '.join(f'{v:.3f}' for v in pf) + f'   | shift {ps[0] - rs.x[0]:+.3f} -> ' + ' '.join(f'{v:.3f}' for v in ps) + f'   -2lnL bin {pnll(A[i] * shape(*us(rs.x, i)), H[i]):.1f}/{len(cen)}')
      E = amp_errors(H, A, rs.x, us, 'shift'); sig = rs.x[3]; c = sig * np.sqrt(2 * np.pi) / BW
      print(f'   shared: positions {rs.x[0]:.3f} {rs.x[1]:.3f} {rs.x[2]:.3f}  sigma {sig:.3f}  tilt {rs.x[4]:+.3f}')
      print('   bin    accept-corrected yields 2.6 3.4 4.0        | dsdo [mb/sr] (env .84, NO trigger term)   | <A> at the peaks')
      for i in range(NB):
         _, _, amps, pos, _ = us(rs.x, i); Y = amps * c; eY = E[i] * c; fct = 1 / (dO[i] * Leff * ENV)
         Ap = [np.interp(p_, cen, A[i]) for p_ in pos]
         rows.append((vz0, TB[i], TB[i + 1], int(H[i].sum()), *Ap, *Y, *eY, *(Y * fct), *(eY * fct), pos[0] - rs.x[0]))
         print(f'   {TB[i]:.0f}-{TB[i+1]:.0f}  ' + '  '.join(f'{y:6.1f}({e:5.1f})' for y, e in zip(Y, eY)) + ' | ' + '  '.join(f'{y * fct:5.2f}({e * fct:4.2f})' for y, e in zip(Y, eY)) + ' | ' + ' '.join(f'{a_:.2f}' for a_ in Ap))
      keep[vz0] = dict(H=H, A=A, xs=rs.x, xf=rf.x)
      print()
   hdr = 'vzLo,thcm_lo,thcm_hi,N,A2.6,A3.4,A4.0,Y2.6,Y3.4,Y4.0,eY2.6,eY3.4,eY4.0,dsdo2.6,dsdo3.4,dsdo4.0,edsdo2.6,edsdo3.4,edsdo4.0,shift'
   with open(U + 'unb_ad3.csv', 'w') as fo:
      fo.write(f'# NOT ADOPTED, NO TRIGGER CORRECTION. unb_ad3.py: acceptance in the model, common shift per bin, Ex {EXLO}-{EXHI}, theta_lab {TLO}-{THI}, vz [vzLo,900), env {ENV}\n' + hdr + '\n')
      for r in rows: fo.write(','.join(f'{v:.5g}' for v in r) + '\n')
   np.save(U + 'unb_ad3_keep.npy', keep, allow_pickle=True)
