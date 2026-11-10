# Spyral validity on the 16C(d,d) elastic locus, data vs sim (2026-11-10).
# The r280 sim acceptance in dd_scale.py includes the Spyral validity of SIMULATED tracks; if data tracks fail Spyral
# more often, the scale S absorbs it. Here every Spyral result matched to a genfit fit is used (dump_dd_valid.C),
# the locus is defined from the FIT kinematics only (Ex within the window around the slice median, theta_lab 60-89.5,
# KE 0.3-10, carbon IC for data), and the valid fraction is compared per theta_cm and per vz.
# Usage: python -I dd_valid_frac.py VALID_DIR OUT_PREFIX
import sys, glob
import numpy as np, uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VD, OUT = sys.argv[1], sys.argv[2]
U = 931.49401; M1 = 16.0147013 * U; Md = 2.0135532 * U
Eb = lambda z: 183.871 - 0.011071 * z
out = []
P = lambda s='': (print(s), out.append(s))


def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))


def kine(Kp, th, Ke):
    M2, M3, M4 = Md, Md, M1
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    A = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    A = np.where(A > 0, A, np.nan); m4 = np.sqrt(A); t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return m4 - M4, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))


def load(files):
    A = [uproot.open(f)['v'].arrays(library='np') for f in files]
    return {k: np.concatenate([x[k] for x in A]) for k in A[0]}


def locus(d, isdata, nsig=2.5):
    ex, cm = kine(Eb(d['vz']), np.radians(d['th']), d['ke'])
    b = (d['th'] > 60) & (d['th'] < 89.5) & (d['ke'] > 0.3) & (d['ke'] < 10) & np.isfinite(ex) & (d['vz'] > 100) & (d['vz'] < 900)
    if isdata:
        b &= (d['ic'] >= 1000) & (d['ic'] <= 1350)
    s = b & (np.abs(ex) < 1.5)
    mu = np.median(ex[s]); q = np.percentile(ex[s], [16, 84]); sg = (q[1] - q[0]) / 2
    pk = b & (np.abs(ex - mu) < nsig * sg)
    sb = b & (ex > -4.5) & (ex < -1.5)   # sideband below the peak (data background under the locus)
    return ex, cm, pk, sb, mu, sg


D = load(sorted(glob.glob(VD + '/run_*.root'))); S = load(sorted(glob.glob(VD + '/sim_*.root')))
exD, cmD, pkD, sbD, muD, sgD = locus(D, True); exS, cmS, pkS, _, muS, sgS = locus(S, False)
P(f'Spyral validity on the elastic locus from the genfit kinematics (dd_valid_frac.py)')
P(f'data: {len(glob.glob(VD + "/run_*.root"))} runs, fitted tracks {len(D["ke"])}; locus Ex {muD:+.3f} +- 2.5x{sgD:.3f}: {pkD.sum()}')
P(f'sim : {len(glob.glob(VD + "/sim_*.root"))} seeds, fitted tracks {len(S["ke"])}; locus Ex {muS:+.3f} +- 2.5x{sgS:.3f}: {pkS.sum()}')
wpk = 2 * 2.5 * sgD; wsb = 3.0
P('\nValid fraction (data sideband-subtracted, sideband Ex -4.5..-1.5 scaled by width) | sim | data/sim')
TB = np.arange(12, 33, 2.)
rows = []
for lo, hi in zip(TB[:-1], TB[1:]):
    md, ms = pkD & (cmD >= lo) & (cmD < hi), pkS & (cmS >= lo) & (cmS < hi)
    bd = sbD & (cmD >= lo) & (cmD < hi)
    nd = md.sum() - bd.sum() * wpk / wsb; vd = (D['valid'][md] > 0).sum() - (D['valid'][bd] > 0).sum() * wpk / wsb
    fd = vd / nd if nd > 20 else np.nan; fs = np.mean(S['valid'][ms] > 0) if ms.sum() > 20 else np.nan
    rows.append((lo, hi, fd, fs))
    P(f'  theta_cm {lo:2.0f}-{hi:<2.0f} data {fd:.3f} (n {md.sum():5d}, bg {bd.sum()*wpk/wsb:6.0f}) | sim {fs:.3f} (n {ms.sum():4d}) | {fd/fs:.3f}')
P('\nBy vz, theta_cm 16-22:')
for z0 in range(100, 900, 100):
    md = pkD & (cmD >= 16) & (cmD < 22) & (D['vz'] >= z0) & (D['vz'] < z0 + 100)
    ms = pkS & (cmS >= 16) & (cmS < 22) & (S['vz'] >= z0) & (S['vz'] < z0 + 100)
    fd, fs = np.mean(D['valid'][md] > 0), np.mean(S['valid'][ms] > 0)
    P(f'  vz {z0}-{z0+100}: data {fd:.3f} (n {md.sum()}) | sim {fs:.3f} (n {ms.sum()}) | {fd/fs:.3f}')
open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
xm = [.5 * (r[0] + r[1]) for r in rows]
ax[0].plot(xm, [r[2] for r in rows], 'o-', label='data (sideband-subtracted)'); ax[0].plot(xm, [r[3] for r in rows], 's-', label='sim r280')
ax[0].set_xlabel(r'$\theta_{cm}$ [deg]'); ax[0].set_ylabel('Spyral valid fraction of fitted locus tracks'); ax[0].legend(); ax[0].grid(alpha=.3)
ax[0].set_ylim(0, 1.05)
ax[1].hist(exD[(cmD > 16) & (cmD < 22) & (D['ic'] >= 1000) & (D['ic'] <= 1350) & (D['th'] > 60)], np.linspace(-5, 3, 161), histtype='step', density=True, label='data, all fits')
ax[1].hist(exD[(cmD > 16) & (cmD < 22) & (D['ic'] >= 1000) & (D['ic'] <= 1350) & (D['th'] > 60) & (D['valid'] > 0)], np.linspace(-5, 3, 161), histtype='step', density=True, label='data, Spyral valid')
ax[1].hist(exS[(cmS > 16) & (cmS < 22) & (S['th'] > 60)], np.linspace(-5, 3, 161), histtype='step', density=True, label='sim, all fits')
ax[1].set_xlabel('Ex from the genfit kinematics [MeV], theta_cm 16-22'); ax[1].legend(fontsize=8); ax[1].set_yscale('log')
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=100)
