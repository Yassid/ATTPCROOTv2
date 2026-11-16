# Extracted strengths of the BW states (unb_bw.csv) vs shell model (YSOX, WBP; files from Yassid 2026-11-16 in shellmodel/).
# S_k = one normalisation of the state's FRESCO ADWA template (C17_dp_fits/fresco/outputs/dp17C_L<L>_state<k>.dat, unit C2S, L from J^pi)
# to its dsigma/dOmega over the theta_cm bins. NO TRIGGER TERM: absolute S carry the open trigger scale; RATIOS to 4.231 largely cancel it.
#   ~/Spyral/venv/bin/python -I unb_sm.py [csv]
import sys, numpy as np
CSV = sys.argv[1] if len(sys.argv) > 1 else '/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw.csv'
F = '/home/yassid/C17_dp_fits/fresco/outputs/dp17C_L%d_state%d.dat'
ST = [('2.763', 1, 1, '1/2-', .0684, .0372), ('2.980', 2, 2, '3/2+', .2043, .2029), ('3.661', 3, 2, 'new (L2 assumed)', None, None),
      ('4.231', 4, 2, '3/2+', .4407, .4909), ('4.841', 5, 0, '1/2+', .0178, .0198), ('5.910', 6, 2, '3/2+', .1401, .1227)]
_L = [l.strip().split(',') for l in open(CSV) if not l.startswith('#')]; D = {k: np.array([float(r[i]) for r in _L[1:]]) for i, k in enumerate(_L[0])}
def avg(g, a, b): t = np.arange(a, b, .05); w = np.sin(np.radians(t)); return np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w)
for vz0 in sorted(set(D['vzLo'])):
   sel = D['vzLo'] == vz0; lo, hi = D['thcm_lo'][sel], D['thcm_hi'][sel]
   print(f'=== vz >= {vz0:.0f}   S = fitted normalisation of the unit-C2S template (NO trigger term)')
   print('   state  J^pi              L   dsdo per bin [mb/sr]                       S(fit)        chi2/ndf |  YSOX    WBP  | S/S(4.231)  YSOX  WBP')
   res = {}
   for name, k, L, jp, ys, wb in ST:
      x, e = D['dsdo' + name][sel], np.maximum(D['edsdo' + name][sel], 1e-3)
      g = np.loadtxt(F % (L, k)); m = np.array([avg(g, a, b) for a, b in zip(lo, hi)])
      S = np.sum(x * m / e ** 2) / np.sum(m * m / e ** 2); eS = 1 / np.sqrt(np.sum(m * m / e ** 2)); c = np.sum(((x - S * m) / e) ** 2)
      res[name] = (S, eS, c, len(x) - 1, ys, wb, jp, L, x, e)
   S4, eS4 = res['4.231'][0], res['4.231'][1]
   for name, (S, eS, c, ndf, ys, wb, jp, L, x, e) in res.items():
      r = S / S4; er = abs(r) * np.sqrt((eS / max(S, 1e-9)) ** 2 + (eS4 / S4) ** 2)
      f = lambda v: f'{v:.3f}' if v is not None else '  -  '
      fr = lambda v: f'{v / .4407:.2f}' if v is not None else '  - '
      fw = lambda v: f'{v / .4909:.2f}' if v is not None else '  - '
      print(f'   {name}  {jp:16s}  {L}   ' + ' '.join(f'{a:5.2f}({b:4.2f})' for a, b in zip(x, e)) + f'   {S:.3f}({eS:.3f})  {c:5.1f}/{ndf} | {f(ys)}  {f(wb)} |  {r:.2f}({er:.2f})   {fr(ys)}  {fw(wb)}')
   print()
