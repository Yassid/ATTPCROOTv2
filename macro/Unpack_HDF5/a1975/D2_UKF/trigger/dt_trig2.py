# (d,t): is the yield per mm of target flat in vz, bin by bin in theta_cm? Acceptance (no trigger term) is 0.94-1.00 and flat, so any vz dependence is a loss the acceptance does not carry.
# Also: the late-fire (triton-path) fraction per cell vs the yield in that cell, the (d,p) signature of a trigger loss.
import numpy as np
from scipy.stats import chi2 as C2
D = '/mnt/f/a1975/sim_rebuild/trig/dt/'; z = np.load(D + 'joined_dt.npz'); T = z['T']; ex = z['ex']; tcm = z['tcm']; early = z['early']; base = z['base']; vz = T[:, 5]
ZB = np.array([25, 125, 225, 325, 425, 525, 600.]); dz = np.diff(ZB); TB = [5, 10, 15, 20, 25, 30, 35, 40, 45]
L = {'g.s.+0.740': (ex > -1) & (ex < 1.3), '3.103': (ex > 2.7) & (ex < 3.6), '4.657': (ex > 4.2) & (ex < 5.2), '6.700': (ex > 6.2) & (ex < 7.3), 'continuum 8-14': (ex > 8) & (ex < 14)}
def test(w, lab):
    n = np.histogram(vz[w], ZB)[0].astype(float); e = early[w]; lam = n.sum() * dz / dz.sum(); c = np.sum((n - lam)**2 / np.maximum(lam, 1)); ndf = len(n) - 1
    # slope of density vs vz, weighted LSQ: fractional change over 25-600
    x = 0.5 * (ZB[1:] + ZB[:-1]); y = n / dz; s = np.sqrt(np.maximum(n, 1)) / dz; W = 1 / s**2; xm = np.sum(W * x) / W.sum(); ym = np.sum(W * y) / W.sum(); b = np.sum(W * (x - xm) * (y - ym)) / np.sum(W * (x - xm)**2); be = 1 / np.sqrt(np.sum(W * (x - xm)**2))
    fr = b * 575 / ym; fe = be * 575 / ym
    ef = [e[(vz[w] >= a) & (vz[w] < b_)].mean() if ((vz[w] >= a) & (vz[w] < b_)).sum() >= 5 else np.nan for a, b_ in zip(ZB[:-1], ZB[1:])]
    print(f'{lab:28s} N {int(n.sum()):5d} per-100mm', ' '.join(f'{v*100:6.1f}' for v in n / dz), f'| chi2 {c:5.1f}/{ndf} p {C2.sf(c, ndf):.3f} | far/near over 575 mm {fr:+.2f}({fe:.2f}) | early', ' '.join(f'{v:.2f}' for v in ef))
sel = base
print('vz slices', ZB, ' counts per 100 mm; flat expected')
for nm, wl in L.items():
    print(f'\n== {nm}')
    test(sel & wl & (tcm < 40), 'theta_cm < 40')
    for a, b in zip(TB[:-1], TB[1:]): test(sel & wl & (tcm >= a) & (tcm < b), f'theta_cm {a}-{b}')
# reference: same test on the (d,p)-like setup is not available here; control = all tritons, Ex 8-14 is continuum
