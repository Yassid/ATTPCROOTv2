# Unbound ADs from the BW fits at offset -0.32 (unb_bw_force32.py keep), L = 0..3 template fits and shell-model comparison. NOT ADOPTED, NO TRIGGER TERM.
# Errors: Fisher matrix of the per-bin linear Poisson model (shared offset/sigma/tilt held). dsdo = Y / (dOmega Leff env 0.84), Y acceptance-corrected.
import numpy as np
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw_force.py').read().split("\nif __name__ == '__main__':")[0])
keep = np.load(U + f'unb_bw_force32{SUF}_keep.npy', allow_pickle=True).item()
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
TPL = [(1, 1, '1/2-', .0684, .0372), (2, 2, '3/2+', .2043, .2029), (3, 2, 'new', None, None), (4, 2, '3/2+', .4407, .4909), (5, 0, '1/2+', .0178, .0198), (6, 2, '3/2+', .1401, .1227)]
VAR = [('all 6 free, offset fixed -0.32', 'FREE'), ('2.980 tied at YSOX ratio, offset -0.32', 'TIED-YSOX')]
def errs(H, A, g, qs, mode, c):
   S = shapes(g[0], g[1], 1.0); E = []
   for i in range(NB):
      ps = np.sqrt(np.maximum(cen - SN, 0)) * (1 + g[2] * (cen - 3))
      if mode == 'free': B = np.vstack([ps, *S])
      else: R = c * XR[i]; B = np.vstack([ps, S[0] + R * S[1], *S[2:]])
      q = qs[i]; Bm = A[i][None, :] * B; m = np.maximum(q @ Bm, 1e-9); I = (Bm / m) @ Bm.T
      act = np.r_[True, q[1:] > 1e-6]; C = np.zeros_like(I); C[np.ix_(act, act)] = np.linalg.pinv(I[np.ix_(act, act)])
      e = np.array([np.sqrt(C[k, k]) if act[k] else 1 / np.sqrt(max(I[k, k], 1e-12)) for k in range(1, len(q))])
      E.append(e if mode == 'free' else np.r_[e[0], R * e[0], e[1:]])
   return np.array(E)
def avg(gt, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, gt[:, 0], gt[:, 1]) * w) / np.sum(w)
rows = []
for vz0 in (100, 500):
   Leff = LPM * (900 - vz0)
   for lab, short in VAR:
      L, g, Y, H, A, qs, mode, fo, c = keep[(vz0, lab)]; E = errs(H, A, g, qs, mode, c)
      fct = 1 / (dO * Leff * ENV); X = Y * fct[:, None]; EX = np.maximum(E * fct[:, None], 1e-3)
      print(f'\n=== vz >= {vz0}, {short}  (fit -2lnL {L:.1f}; sigma {g[1]:.3f})   dsdo [mb/sr] per theta_cm bin 22-26 26-30 30-34 34-38, NO trigger term')
      print('   state        dsdo                                              |  L=0         L=1         L=2         L=3      | S at SM L    YSOX   WBP')
      Sres = {}
      for s, (k, Lsm, jp, ys, wb) in enumerate(TPL):
         x, e = X[:, s], EX[:, s]; u = np.isfinite(x)
         if True: u = u & (np.array([np.interp(p_, cen, A[i]) for i, p_ in enumerate([ST[s][0] + g[0]] * NB)]) > 0.3)   # drop bins where the state is outside the acceptance
         line = f'   {ST[s][3]:13s}' + ' '.join(f'{a_:5.2f}({b_:4.2f})' if uu else '    -      ' for a_, b_, uu in zip(x, e, u)) + ' |'
         fits = []
         for Lt in range(4):
            gt = np.loadtxt(F % (Lt, k)); m = np.array([avg(gt, a, b) for a, b in zip(TB[:-1], TB[1:])])[u]
            Sv = np.sum(x[u] * m / e[u] ** 2) / np.sum(m * m / e[u] ** 2); eS = 1 / np.sqrt(np.sum(m * m / e[u] ** 2)); ch = np.sum(((x[u] - Sv * m) / e[u]) ** 2)
            fits.append((Sv, eS, ch)); line += f' {ch:4.1f}/{u.sum() - 1}     '
         Ssm = fits[Lsm]; Sres[s] = Ssm
         line += f'| {Ssm[0]:.3f}({Ssm[1]:.3f}) ' + (f'{ys:.3f} {wb:.3f}' if ys else '  -     -  ')
         print(line)
         rows.append((vz0, short, ST[s][3], *x, *e, *[f[2] for f in fits], Ssm[0], Ssm[1]))
      S4 = Sres[3]
      print('   ratios S/S(4.231):  ' + '  '.join(f'{ST[s][3].split()[0]} {Sres[s][0] / S4[0]:.2f}({abs(Sres[s][0] / S4[0]) * np.hypot(Sres[s][1] / max(Sres[s][0], 1e-9), S4[1] / S4[0]):.2f})'
                                         + (f' [YSOX {TPL[s][3] / .4407:.2f} WBP {TPL[s][4] / .4909:.2f}]' if TPL[s][3] else '') for s in (0, 1, 2)))
hdr = 'vzLo,variant,state,' + ','.join(f'dsdo{b}' for b in range(NB)) + ',' + ','.join(f'err{b}' for b in range(NB)) + ',chi2L0,chi2L1,chi2L2,chi2L3,S_smL,eS_smL'
with open(U + f'unb_ad32{SUF}.csv', 'w') as fo:
   fo.write(f'# NOT ADOPTED, NO TRIGGER TERM. BW offset -0.32, theta_cm bin edges {list(TB)}, env 0.84; unb_ad32.py\n' + hdr + '\n')
   for r in rows: fo.write(','.join(str(v) if isinstance(v, str) else f'{v:.5g}' for v in r) + '\n')
