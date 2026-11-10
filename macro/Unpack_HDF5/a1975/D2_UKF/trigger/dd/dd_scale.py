# 16C(d,d) elastic scale test with the r280 acceptance (2026-11-10).
# Data: dump_dd.C caches (carbon IC gate, runs given), luminosity L over 971.73 mm.
# Sim:  dump_dd_sim.C outputs of the elastic r280 simulation (run_dd_el.sh): "g" generated, "d" reconstructed.
# Same selection on both: theta_lab 60-89.5, KE 0.3-10 MeV, Ex within +-2.5 sigma of the slice median (each sample
# its own median and sigma), lower sideband subtraction for data. Acceptance per (vz, theta_cm) cell =
# selected reco (reco vz, reco theta_cm) / generated (true z + 19.1 mm, true theta_cm = 180 - 2 theta_lab).
# S = Y / (L dz/971.73 * eps * sigma_OM * dOmega): the data/theory scale per cell and per potential.
# Usage: python -I dd_scale.py DATA_CACHE_DIR L_MB SIM_DUMP_GLOB OUT_PREFIX EXCLUDE_RUNS(comma or -) [NSIG=2.5]
# Data Ex peak degrades beyond theta_cm 22 (shifts down, dissolves at 26+; sim does not): quote only theta_cm <= 22.
import sys, glob
import numpy as np, uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CDIR, LTOT, SGLOB, OUT = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
EXCL = [int(x) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 and sys.argv[5] != '-' else []
NSIG = float(sys.argv[6]) if len(sys.argv) > 6 else 2.5   # Ex window half-width in units of each sample's slice sigma
LEN = 971.73; ZOFF = 19.1; ICLO, ICHI = 1000, 1350
U = 931.49401; M1 = 16.0147013 * U; Md = 2.0135532 * U
Eb = lambda z: 183.871 - 0.011071 * z
OMD = '/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/trigger/dd_om/'
POT = {'A': 'An-Cai', 'H': 'Han-Shi-Shen', 'Z': 'Zhang-Pang-Lou', 'D': 'Daehnick'}
EGRID = [21.9, 22.4, 23.1]
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


def load(files, tree):
    A = [uproot.open(f)[tree].arrays(library='np') for f in files]
    return {k: np.concatenate([x[k] for x in A]) for k in A[0]}


ZB = np.arange(100, 901, 100.); TB = np.arange(12, 31, 2.)
dO = 2*np.pi*(np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))


def select(d, sideband=False):
    ex, cm = kine(Eb(d['vz']), np.radians(d['th']), d['ke'])
    base = (d['ic'] >= ICLO) & (d['ic'] <= ICHI) & (d['th'] > 60) & (d['th'] < 89.5) & (d['ke'] > 0.3) & (d['ke'] < 10) & np.isfinite(ex)
    iz = np.clip(np.digitize(d['vz'], ZB) - 1, 0, len(ZB) - 2)
    mu = np.zeros(len(ZB) - 1); sg = np.zeros(len(ZB) - 1)
    sall = base & (np.abs(ex) < 1.5)
    q = np.percentile(ex[sall], [16, 50, 84]) if sall.sum() else np.array([-.25, 0, .25])
    for j in range(len(ZB) - 1):
        s = base & (iz == j) & (np.abs(ex) < 1.5) & (d['vz'] >= ZB[j]) & (d['vz'] < ZB[j + 1])
        if s.sum() < 30:  # sparse slice: global peak
            mu[j], sg[j] = q[1], (q[2] - q[0]) / 2; continue
        mu[j] = np.median(ex[s]); qq = np.percentile(ex[s], [16, 84]); sg[j] = (qq[1] - qq[0]) / 2
    pk = base & (np.abs(ex - mu[iz]) < NSIG * sg[iz])
    sb = base & (ex > -4.5) & (ex < -1.5)
    return ex, cm, pk, sb, sg, mu


files = [f for f in sorted(glob.glob(CDIR + '/run_*.root')) if int(f.split('/')[-1][4:8]) not in EXCL]
D = load(files, 'd'); exD, cmD, pkD, sbD, sgD, muD = select(D)
sims = sorted(glob.glob(SGLOB))
G = load(sims, 'g'); S = load(sims, 'd'); exS, cmS, pkS, _, sgS, muS = select(S)
gok = G['tke'] > 0
gvz = G['tvz'] + ZOFF; gcm = 180 - 2 * G['tth']
P(f'16C(d,d) scale with r280 acceptance. Data: {len(files)} runs ({",".join(f.split("/")[-1][4:8] for f in files)}), L = {LTOT} mb^-1, IC {ICLO}-{ICHI}.')
P(f'Sim: {len(sims)} files, {gok.sum()} generated elastic events, {len(S["ke"])} reco rows, {pkS.sum()} selected.')
P('Ex peak sigma per vz slice (100-900): data ' + ' '.join(f'{x:.3f}' for x in sgD) + ' | sim ' + ' '.join(f'{x:.3f}' for x in sgS))
P('Ex peak median: data ' + ' '.join(f'{x:+.3f}' for x in muD) + ' | sim ' + ' '.join(f'{x:+.3f}' for x in muS))

T = {k: [np.loadtxt(OMD + f'dd_{k}_{str(e).replace(".", "p")}_dsdo.dat', comments=['#', '&', '@']) for e in EGRID] for k in POT}


def sig_bin(k, Ed, a, b):
    t = np.arange(a, b, .05); w = np.sin(np.radians(t))
    v = [np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w) for g in T[k]]
    return np.interp(Ed, EGRID, v)


eps = np.zeros((len(ZB) - 1, len(TB) - 1)); epe = np.zeros_like(eps)
Sc = {k: np.full_like(eps, np.nan) for k in POT}; Se = np.full_like(eps, np.nan); Ycell = np.zeros_like(eps)
Den = {k: np.zeros_like(eps) for k in POT}; Num = np.zeros_like(eps); Var = np.zeros_like(eps)
P('\nAcceptance x efficiency (sim, r280), rows vz, columns theta_cm (12-30 by 2):')
for j, (z0, z1) in enumerate(zip(ZB[:-1], ZB[1:])):
    for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
        ng = np.sum(gok & (gvz >= z0) & (gvz < z1) & (gcm >= a) & (gcm < b))
        nr = np.sum(pkS & (S['vz'] >= z0) & (S['vz'] < z1) & (cmS >= a) & (cmS < b))
        eps[j, i] = nr / ng if ng else np.nan; epe[j, i] = np.sqrt(max(eps[j, i] * (1 - eps[j, i]), 1 / ng) / ng) if ng else np.nan
    P(f'  {z0:3.0f}-{z1:<4.0f} ' + ' '.join(f'{x:4.2f}' for x in eps[j]))

P('\nScale S = Y/(L eps sigma_ZPL dOmega), carbon; "-" = Y < 10 or eps < 0.3. Rows vz, columns theta_cm:')
for j, (z0, z1) in enumerate(zip(ZB[:-1], ZB[1:])):
    w = (D['vz'] >= z0) & (D['vz'] < z1)
    Y = np.histogram(cmD[w & pkD], TB)[0].astype(float)
    B = np.histogram(cmD[w & sbD], TB)[0].astype(float) * (2 * NSIG * sgD[j]) / 3.0
    L = LTOT / LEN * (z1 - z0); Ed = Eb(.5 * (z0 + z1)) * Md / M1
    cells = []
    for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
        if Y[i] < 10 or not eps[j, i] >= 0.3:
            cells.append('    -     '); continue
        th = {k: sig_bin(k, Ed, a, b) for k in POT}
        for k in POT:
            Sc[k][j, i] = (Y[i] - B[i]) / (L * eps[j, i] * th[k] * dO[i])
            Den[k][j, i] = L * eps[j, i] * th[k] * dO[i]
        Num[j, i] = Y[i] - B[i]; Var[j, i] = Y[i] + B[i]
        Se[j, i] = Sc['Z'][j, i] * np.hypot(np.sqrt(Y[i] + B[i]) / max(Y[i] - B[i], 1), epe[j, i] / eps[j, i])
        Ycell[j, i] = Y[i]
        cells.append(f'{Sc["Z"][j, i]:5.3f}({int(round(1000*Se[j, i])):3d})')
    P(f'  {z0:3.0f}-{z1:<4.0f} ' + ' '.join(cells))

P('\nScale per window and potential (stat incl. acceptance stat):')
for lab, zlo, zhi, tlo, thi in [('plateau vz 250-900, theta_cm 16-22', 250, 900, 16, 22), ('vz 200-900, theta_cm 14-24', 200, 900, 14, 24),
                                ('vz 300-900, theta_cm 18-22', 300, 900, 18, 22), ('vz 200-900, theta_cm 22-28', 200, 900, 22, 28)]:
    m = np.zeros_like(eps, bool)
    for j, (z0, z1) in enumerate(zip(ZB[:-1], ZB[1:])):
        for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
            m[j, i] = z0 >= zlo and z1 <= zhi and a >= tlo and b <= thi and np.isfinite(Se[j, i])
    s = []; r = []
    for k, n in POT.items():
        v = Sc[k][m]; e = Se[m] * Sc[k][m] / Sc['Z'][m]
        s.append(f'{n} {np.sum(v/e**2)/np.sum(1/e**2):.3f}')
        # ratio of sums: linear in eps, so the per-cell acceptance noise does not bias it (the weighted mean is
        # biased LOW: a cell with eps fluctuating up gets a low S AND a small error)
        den = Den[k][m]; num = Num[m]
        Sr = num.sum() / den.sum()
        er = Sr * np.hypot(np.sqrt(Var[m].sum()) / num.sum(), np.sqrt(np.sum((den * epe[m] / eps[m])**2)) / den.sum())
        r.append(f'{n} {Sr:.3f}({er:.3f})')
    v = Sc['Z'][m]; e = Se[m]; mean = np.sum(v/e**2)/np.sum(1/e**2); chi = np.sum(((v-mean)/e)**2)
    P(f'  {lab} ({m.sum()} cells):\n     ratio of sums: ' + '  '.join(r) + f'\n     weighted mean (biased low): ' + '  '.join(s) + f'   [ZPL chi2 {chi:.1f}/{m.sum()-1}]')
open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
cmap = plt.get_cmap('viridis')
for j, (z0, z1) in enumerate(zip(ZB[:-1], ZB[1:])):
    xm = .5 * (TB[:-1] + TB[1:]) + (j - 4) * .12
    ok = np.isfinite(Se[j])
    ax[0].errorbar(xm[ok], Sc['Z'][j][ok], Se[j][ok], fmt='o', ms=4, capsize=2, color=cmap(j / 8), label=f'vz {z0:.0f}-{z1:.0f}')
    ax[1].plot(.5 * (TB[:-1] + TB[1:]), eps[j], 'o-', ms=3, color=cmap(j / 8), label=f'vz {z0:.0f}-{z1:.0f}')
ax[0].axhline(1, color='k', ls=':', lw=.8); ax[0].set_xlabel(r'$\theta_{cm}$ [deg]'); ax[0].set_ylabel('data / Zhang-Pang-Lou')
ax[0].set_title('16C(d,d) scale, r280 acceptance applied', fontsize=10); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=.3); ax[0].set_ylim(0, 2)
ax[1].set_xlabel(r'$\theta_{cm}$ [deg]'); ax[1].set_ylabel('acceptance x efficiency'); ax[1].set_title('elastic sim r280', fontsize=10)
ax[1].legend(fontsize=7, ncol=2); ax[1].grid(alpha=.3); ax[1].set_ylim(0, 1.1)
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=100)
