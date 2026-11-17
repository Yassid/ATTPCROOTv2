# 16C(d,p)17C unbound states, step 2: STATE-SPECIFIC efficiency + vz-threshold scan (trigger test). NOT ADOPTED.
# Efficiency per state E and theta_cm bin: take the generated events of the Ex 3.9 sim (truth vz, truth theta_cm), move each to state E
# at the same vz and theta_cm (two-body kinematics), accept if theta_lab in [TLO,THI) and KE >= 1.2, weight by the reconstruction
# efficiency looked up in (true theta_lab, true vz) from the 3.9 sim (genfit + chi2 + reco theta_lab window), times Spyral validity.
# theta_cm migration ignored (checked below on the sim). VALIDATED on the g.s. sim against the directly counted efficiency (mode 'val').
#   ~/Spyral/venv/bin/python -I unb_ad2.py val | scan
import sys, os, numpy as np, uproot
from scipy.optimize import minimize
exec(open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("K = np.loadtxt")[0])
T = '/mnt/f/a1975/sim_rebuild/trig_hdbtau/'; U = '/mnt/f/a1975/sim_rebuild/unb_hdbtau/'; R5 = '/mnt/f/a1975/sim_rebuild/r280_hdbtau/'
CAL = (184.25, 11.0, 0.18, -3, -0.05, -0.6947); LPM = 316.98 / 971.73; ENV = 0.84; MP = 938.272
TLO, THI = 89, float(os.environ.get('THI', 113)); KEMIN = float(os.environ.get('KEMIN', 1.2))   # KEMIN: 2026-11-17, sim shows KE 0.8-1.2 reconstructs as well as > 1.2   # 2026-11-17: THI=137 extends the theta_lab window (acceptance is in the model)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'scan'

def forward(Kb, tcm_deg, ex):
   """beam KE, proton theta_cm (deg, the convention of kine(): 180 - CM angle of the proton), Ex -> proton theta_lab (deg), KE."""
   m1, m2, m3, m4 = M1, M2, M3, M4 + ex
   E1 = Kb + m1; s = m1 * m1 + m2 * m2 + 2 * m2 * E1; rs = np.sqrt(s)
   p1 = np.sqrt(E1 ** 2 - m1 ** 2); beta = p1 / (E1 + m2); gam = 1 / np.sqrt(1 - beta ** 2)
   pc = om(s, m3 * m3, m4 * m4) / (2 * rs); E3c = np.sqrt(pc ** 2 + m3 ** 2)
   a = np.radians(180 - tcm_deg)           # CM angle of the proton w.r.t. the beam
   pz = gam * (pc * np.cos(a) + beta * E3c); pt = pc * np.sin(a); E3 = gam * (E3c + beta * pc * np.cos(a))
   return np.degrees(np.arctan2(pt, pz)), E3 - m3

# check forward() inverts kine()
_t, _k = forward(np.array([180.]), np.array([30.]), np.array([3.4])); _x, _c = kine(np.array([180.]), np.radians(_t), _k)
assert abs(_x[0] - 3.4) < 1e-6 and abs(_c[0] - 30) < 1e-6, (_x, _c)

def reco_eff_table(keys, dump, vzb, thb):
   K = np.loadtxt(keys); gz = K[:, 1] + ZOFF
   s = uproot.open(dump)['r'].arrays(['ke', 'th', 'vz', 'c2', 'tth', 'tvz'], library='np')
   ok = (s['c2'] <= 13.5) & (s['th'] >= TLO) & (s['th'] < THI) & (s['ke'] >= KEMIN)
   Ng = np.histogram2d(gz, K[:, 2], [vzb, thb])[0]; Nr = np.histogram2d(s['tvz'][ok] + ZOFF, s['tth'][ok], [vzb, thb])[0]
   return Nr / np.maximum(Ng, 1), Ng

VZB = np.arange(0, 1001, 100.); THB = np.arange(85, THI + 5, 2.)
EFF39, NG39 = reco_eff_table(U + 'keys_ex39.txt', R5 + 'p_ex39_dump.root', VZB, THB)
# 2026-11-17: where the 3.9 sim has no generated protons (theta_lab >~ 121 cannot be reached by Ex 3.9 with KE >= 1.2) fall back to the g.s. sim
# table; reco eff is flat 0.85-0.94 above theta_lab 101 in both sims, whatever KE. Inside 89-113 nothing changes (3.9 table populated).
EFFGS, NGGS = reco_eff_table(T + 'keys.txt', R5 + 'p_gs_dump.root', VZB, THB); EFF39 = np.where(NG39 >= 10, EFF39, EFFGS)
K39 = np.loadtxt(U + 'keys_ex39.txt'); gz39 = K39[:, 1] + ZOFF; _, gc39 = kine(Eb(gz39), np.radians(K39[:, 2]), K39[:, 3])
X = uproot.open(R5 + 'p_ex39_census.root')['c'].arrays(['thl', 'vz', 'ke', 'ex', 'sv'], library='np')
SPB = np.arange(89, THI + 1.001, 4.); xs = (X['ke'] >= KEMIN) & (X['ex'] > 2.0) & (X['ex'] < 6) & (X['vz'] > 0) & (X['vz'] < 1000)
SPY = np.array([np.mean(X['sv'][xs & (X['thl'] >= a) & (X['thl'] < b)] > 0) for a, b in zip(SPB[:-1], SPB[1:])])
XG = uproot.open(R5 + 'p_gs_census.root')['c'].arrays(['thl', 'ke', 'sv'], library='np')   # fallback where the 3.9 census is empty
SPYGS = np.array([np.mean(XG['sv'][(XG['ke'] >= KEMIN) & (XG['thl'] >= a) & (XG['thl'] < b)] > 0) for a, b in zip(SPB[:-1], SPB[1:])]); SPY = np.where(np.isfinite(SPY), SPY, SPYGS)

def eps_state(ex, TB, vz0, vz1, gz=gz39, gc=gc39, eff=EFF39, spy=True):
   w = (gz >= vz0) & (gz < vz1) & np.isfinite(gc); thl, kp = forward(Eb(gz[w]), gc[w], ex)
   iv = np.clip(np.digitize(gz[w], VZB) - 1, 0, len(VZB) - 2); it = np.digitize(thl, THB) - 1
   inside = (thl >= TLO) & (thl < THI) & (kp >= KEMIN) & (it >= 0) & (it < len(THB) - 1)
   wt = np.zeros(w.sum()); wt[inside] = eff[iv[inside], it[inside]]
   if spy: wt[inside] *= SPY[np.clip(np.digitize(thl[inside], SPB) - 1, 0, len(SPY) - 1)]
   return np.histogram(gc[w], TB, weights=wt)[0] / np.maximum(np.histogram(gc[w], TB)[0], 1)

if MODE == 'val':
   # g.s. sim: direct efficiency (reco theta_cm bin, same cuts, no Ex window) vs the transported 3.9 table applied at Ex 0 (no Spyral)
   TB = np.arange(float(os.environ.get('VALLO', 24)), 47, 2.)
   for vz0, vz1 in ((100, 900), (100, 400), (400, 900)):
      K = np.loadtxt(T + 'keys.txt'); gz = K[:, 1] + ZOFF; gx, gc = kine(Eb(gz), np.radians(K[:, 2]), K[:, 3])
      g = np.loadtxt(T + 'gf_sim.txt'); sx, sc = kine(Eb(g[:, 2]), np.radians(g[:, 1]), g[:, 0])
      sel = (g[:, 3] <= 13.5) & (g[:, 1] >= TLO) & (g[:, 1] < THI) & (g[:, 0] >= KEMIN) & np.isfinite(sx) & (g[:, 2] >= vz0) & (g[:, 2] < vz1)
      direct = np.histogram(sc[sel], TB)[0] / np.maximum(np.histogram(gc[(gz >= vz0) & (gz < vz1) & np.isfinite(gc)], TB)[0], 1)
      transp = eps_state(0.0, TB, vz0, vz1, spy=False)
      own = eps_state(0.0, TB, vz0, vz1, gz=gz, gc=gc, eff=reco_eff_table(T + 'keys.txt', '/mnt/f/a1975/sim_rebuild/r280_hdbtau/p_gs_dump.root', VZB, THB)[0], spy=False)
      print(f'g.s. vz {vz0}-{vz1}  theta_cm bins {TB[0]:.0f}..{TB[-1]:.0f} by 2')
      print('   direct      ', np.round(direct, 3)); print('   own table   ', np.round(own, 3)); print('   3.9 table   ', np.round(transp, 3))
   # theta_cm resolution of the 3.9 sim
   s = uproot.open(R5 + 'p_ex39_dump.root')['r'].arrays(['ke', 'th', 'vz', 'c2', 'tke', 'tth', 'tvz'], library='np')
   _, rc = kine(Eb(s['vz']), np.radians(s['th']), s['ke']); _, tc = kine(Eb(s['tvz'] + ZOFF), np.radians(s['tth']), s['tke'])
   ok = (s['c2'] <= 13.5) & (s['th'] >= TLO) & (s['th'] < THI) & (s['ke'] >= KEMIN) & np.isfinite(rc) & np.isfinite(tc)
   print('3.9 sim theta_cm reco - true: p16/p50/p84', np.round(np.percentile((rc - tc)[ok], [16, 50, 84]), 2))
   print('Spyral validity by theta_lab', SPB, np.round(SPY, 3))
   sys.exit()

# --- scan: per-bin fits as unb_ad.py, efficiency per state, vz threshold scan
TB = np.array([22, 26, 30, 34, 38.]); dO = 2 * np.pi * (np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))
EXLO, EXHI, BW = 1.0, 5.0, 0.1; edges = np.arange(EXLO, EXHI + 1e-9, BW); cen = .5 * (edges[1:] + edges[:-1])
d = np.loadtxt(T + 'gf_data.txt'); ke, th, vz, c2 = d.T.copy()
ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); dx, dc = kine(CAL[0] - CAL[1] / 1000 * (vz - 500), np.radians(th), ke); dx = dx + CAL[5]
dsel = (c2 <= 13.5) & (th >= TLO) & (th < THI) & (ke >= KEMIN) & np.isfinite(dx)
def model(p, n):
   m = p[0] * np.sqrt(np.maximum(cen - 0.729, 0)) * (1 + p[1] * (cen - 3))
   for k in range(n): m = m + p[3 + 2 * k] * np.exp(-.5 * ((cen - p[4 + 2 * k]) / p[2]) ** 2)
   return np.maximum(m, 1e-9)
def nll(p, h, n): m = model(p, n); return 2 * np.sum(m - h + np.where(h > 0, h * np.log(np.maximum(h, 1e-12) / m), 0))
w = dsel & (vz >= 100) & (vz < 900); hg = np.histogram(dx[w], edges)[0].astype(float)
G = minimize(nll, [hg.mean() / 1.5, .2, .22, 20, 2.6, 30, 3.45, 25, 4.0], args=(hg, 3), method='L-BFGS-B',
             bounds=[(0, None), (-.6, 2), (.08, .9), (0, None), (2.2, 3.0), (0, None), (3.0, 3.8), (0, None), (3.7, 4.6)]); pg = G.x
EXS = [pg[4], pg[6], pg[8]]
print(f'global (vz 100-900): Ex {pg[4]:.3f} {pg[6]:.3f} {pg[8]:.3f} sigma {pg[2]:.3f} tilt {pg[1]:+.3f}')
def binfit(h):
   f = lambda q: nll([q[0], pg[1], pg[2], q[1], pg[4], q[2], pg[6], q[3], pg[8]], h, 3)
   r = minimize(f, [max(h.mean() / 1.5, .1)] + [max(h.max() / 3, .5)] * 3, bounds=[(0, None)] * 4, method='L-BFGS-B'); q = r.x
   e = np.zeros(3)
   for k in range(1, 4):   # 1-sigma from the -2lnL profile along the amplitude (others re-minimised)
      def prof(v):
         g = lambda qq: f(np.insert(qq, k, v)); rr = minimize(g, np.delete(q, k), bounds=[(0, None)] * 3, method='L-BFGS-B'); return rr.fun
      lo, hi = q[k], q[k] + 1.0
      while prof(hi) - r.fun < 1 and hi < 1e4: hi *= 1.5
      for _ in range(30):
         mid = .5 * (lo + hi); lo, hi = (mid, hi) if prof(mid) - r.fun < 1 else (lo, mid)
      e[k - 1] = hi - q[k]
   c = pg[2] * np.sqrt(2 * np.pi) / BW; return q[1:] * c, e * c, r.fun, int(h.sum())
out = []
print('\nvz>=   bin    N  -2lnL/36 | eps 2.6 3.4 4.0       | dsdo [mb/sr]  2.6            3.4            4.0       (NO trigger term)')
for vz0 in (100, 200, 300, 400, 500, 600, 700):
   vz1 = 900; ww = dsel & (vz >= vz0) & (vz < vz1); Leff = LPM * (vz1 - vz0)
   E = np.array([eps_state(x, TB, vz0, vz1) for x in EXS])
   for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
      h = np.histogram(dx[ww & (dc >= a) & (dc < b)], edges)[0].astype(float); Y, eY, L, n = binfit(h)
      f = 1 / (E[:, i] * dO[i] * Leff * ENV); x, ex_ = Y * f, eY * f
      out.append((vz0, a, b, n, *E[:, i], *Y, *eY, *x, *ex_))
      print(f'{vz0:4d}  {a:.0f}-{b:.0f}  {n:4d}  {L:5.1f}   | {E[0, i]:.2f} {E[1, i]:.2f} {E[2, i]:.2f}   | ' + '  '.join(f'{v:5.2f}({e:4.2f})' for v, e in zip(x, ex_)))
   print()
hdr = 'vzLo,thcm_lo,thcm_hi,N,eps2.6,eps3.4,eps4.0,Y2.6,Y3.4,Y4.0,eY2.6,eY3.4,eY4.0,dsdo2.6,dsdo3.4,dsdo4.0,edsdo2.6,edsdo3.4,edsdo4.0'
with open(U + 'unb_vzscan.csv', 'w') as fo:
   fo.write(f'# NOT ADOPTED, NO TRIGGER CORRECTION. unb_ad2.py: vz [vzLo,900), theta_lab {TLO}-{THI}, state-specific eps, env {ENV}; global Ex {EXS} sigma {pg[2]:.3f}\n' + hdr + '\n')
   for r in out: fo.write(','.join(f'{v:.5g}' for v in r) + '\n')
