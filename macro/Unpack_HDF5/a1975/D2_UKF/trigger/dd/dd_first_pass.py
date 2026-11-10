# 16C(d,d) elastic scale test, first pass on the runs dumped by dump_dd.C (2026-11-10).
# Yields per (vz slice, theta_cm bin), with and without the trigger-safe cut r_res >= 28 mm, and the ratio
#     R = Y / (L(vz slice) * sigma_OM(theta_cm, E_d(vz)) * dOmega)
# R is (acceptance x efficiency) x (luminosity scale): NO acceptance is applied, because the elastic simulation
# has not yet been fitted with the r280 geometry. Where the acceptance is ~1 and flat, R is the scale.
# The beam is ~62 % 16C and ~31 % 21O at the same velocity; d+21O elastic recoils sit on the same (theta, KE) locus,
# so the 16C is selected with the IC amplitude (carbon window ICLO-ICHI). L must be the CARBON luminosity of the runs
# kept. Cells the r_res cut only partly keeps are corrected by the kept fraction of the cell (uniform vz,
# sigma*sin(theta) weighting), as the (d,p) fiducial analysis does for its luminosity.
# Usage: python -I dd_first_pass.py CACHE_DIR L_MB OUT_PREFIX [ICLO ICHI] [EXCLUDE_RUNS,comma]
import sys, glob
import numpy as np, uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CDIR, LTOT, OUT = sys.argv[1], float(sys.argv[2]), sys.argv[3]
ICLO, ICHI = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (-1e9, 1e9)
EXCL = [int(x) for x in sys.argv[6].split(',')] if len(sys.argv) > 6 else []
LEN = 971.73                      # mm, the length LTOT refers to (as in the (d,p) luminosity)
U = 931.49401; M1 = 16.0147013 * U; Md = 2.0135532 * U
Eb = lambda z: 183.871 - 0.011071 * z         # 16C energy at the vertex (data calibration, vz_fit frame)
OMD = '/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/trigger/dd_om/'
POT = {'A': 'An-Cai', 'H': 'Han-Shi-Shen', 'Z': 'Zhang-Pang-Lou', 'D': 'Daehnick'}
EGRID = [21.9, 22.4, 23.1]


def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))


def kine(Kp, th, Ke):
    """Ex of the 16C partner and theta_cm (d scattering angle), elastic masses."""
    M2, M3, M4 = Md, Md, M1
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    A = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    A = np.where(A > 0, A, np.nan); m4 = np.sqrt(A); t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return m4 - M4, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))


def rres(ke, th, vz):
    """Radius of the 16C residual at the pad plane, from the measured deuteron (as dsigma_fid.py)."""
    t = np.radians(th); p = np.sqrt(ke**2 + 2*ke*Md); E = Eb(vz); pb = np.sqrt(E**2 + 2*E*M1)
    return (1000 - vz) * np.tan(np.arctan2(p*np.sin(t), pb - p*np.cos(t)))


files = [f for f in sorted(glob.glob(CDIR + '/run_*.root')) if int(f.split('/')[-1][4:8]) not in EXCL]
A_ = [uproot.open(f)['d'].arrays(library='np') for f in files]
d = {k: np.concatenate([x[k] for x in A_]) for k in A_[0]}
ex, cm = kine(Eb(d['vz']), np.radians(d['th']), d['ke'])
rr = rres(d['ke'], d['th'], d['vz'])
base = (d['ic'] >= ICLO) & (d['ic'] <= ICHI) & (d['th'] > 60) & (d['th'] < 89.5) & (d['ke'] > 0.3) & (d['ke'] < 10) & np.isfinite(ex)

out = []
P = lambda s='': (print(s), out.append(s))
P(f'16C(d,d) first pass: {len(files)} runs ({", ".join(f.split("/")[-1][4:8] for f in files)}), {len(d["ke"])} tracks, L = {LTOT} mb^-1 over {LEN} mm')
P(f'IC window {ICLO:g}-{ICHI:g}; excluded runs {EXCL}.')
P('No acceptance applied. R = Y / (L_slice * sigma_OM * dOmega) = (acceptance x efficiency) x (luminosity scale).')

# Ex peak: centroid and width per vz slice, Gaussian fit on the core
ZB = np.arange(0, 901, 100.)
P('\nEx peak per vz slice (theta_lab 60-89.5, all theta_cm): median and robust sigma of |Ex|<1.5')
mu, sg = [], []
for z0, z1 in zip(ZB[:-1], ZB[1:]):
    s = base & (d['vz'] >= z0) & (d['vz'] < z1) & (np.abs(ex) < 1.5)
    m = np.median(ex[s]); q = np.percentile(ex[s], [16, 84]); mu.append(m); sg.append((q[1] - q[0]) / 2)
    P(f'  vz {z0:4.0f}-{z1:4.0f}: n {s.sum():6d}  median {m:+.3f}  sigma {sg[-1]:.3f} MeV')
mu, sg = np.array(mu), np.array(sg)
iz = np.clip(np.digitize(d['vz'], ZB) - 1, 0, len(ZB) - 2)
pull = (ex - mu[iz]) / sg[iz]
WIN = 2.5                          # +-2.5 sigma peak window
SB = (ex > -4.5) & (ex < -1.5)     # lower sideband only: the 1.766 MeV 2+ of 16C sits on the upper side
inpk = np.abs(pull) < WIN
P(f'Peak window +-{WIN} sigma around the slice median; background from the Ex -4.5..-1.5 sideband scaled by width.')


def theory():
    """sigma[pot][E index] as callable of theta_cm (deg)."""
    T = {}
    for k in POT:
        T[k] = [np.loadtxt(OMD + f'dd_{k}_{str(e).replace(".", "p")}_dsdo.dat', comments=['#', '&', '@']) for e in EGRID]
    return T


T = theory()


def sig_bin(k, Ed, a, b):
    """Bin average (sin-weighted) at Ed, linear in E between the three curves."""
    t = np.arange(a, b, .05); w = np.sin(np.radians(t))
    v = [np.sum(np.interp(t, g[:, 0], g[:, 1]) * w) / np.sum(w) for g in T[k]]
    return np.interp(Ed, EGRID, v)


TB = np.arange(12, 31, 2.)
dO = 2*np.pi*(np.cos(np.radians(TB[:-1])) - np.cos(np.radians(TB[1:])))
LPM = LTOT / LEN
KMAX = 4 * M1 * Md / (M1 + Md)**2


def kept_frac(R, z0, z1, a, b, k='Z'):
    """Fraction of a (vz, theta_cm) cell with r_res >= R, elastic kinematics, sigma*sin weighting."""
    if R <= 0:
        return 1.0
    z = np.linspace(z0, z1, 41)[:, None]; t = np.linspace(a, b, 41)[None, :]
    Ed = Eb(z) * Md / M1; w = np.sin(np.radians(t)) * np.interp(t, T[k][1][:, 0], T[k][1][:, 1]) * np.ones_like(z)
    thl = (180 - t) / 2; ke = KMAX * Eb(z) * np.cos(np.radians(thl))**2
    return float(np.sum(w * (rres(ke, thl, z) >= R)) / np.sum(w))


def table(R, label):
    P(f'\n=== {label} ===')
    P('R = Y/(L f sigma dOmega) with Zhang-Pang-Lou, f = kept fraction of the cell; "-" = fewer than 10 counts or f < 0.25. Rows vz slice, columns theta_cm.')
    P('vz       ' + ' '.join(f'{a:4.0f}-{b:<3.0f}' for a, b in zip(TB[:-1], TB[1:])))
    res = {}
    for j, (z0, z1) in enumerate(zip(ZB[:-1], ZB[1:])):
        w = base & (d['vz'] >= z0) & (d['vz'] < z1) & (rr >= R)
        Y = np.histogram(cm[w & inpk], TB)[0].astype(float)
        B = np.histogram(cm[w & SB], TB)[0].astype(float) * (2*WIN*sg[j]) / 3.0
        Ed = Eb(.5*(z0 + z1)) * Md / M1
        L = LPM * (z1 - z0)
        row = {}
        cells = []
        for i, (a, b) in enumerate(zip(TB[:-1], TB[1:])):
            th = {k: sig_bin(k, Ed, a, b) for k in POT}
            f = kept_frac(R, z0, z1, a, b)
            if Y[i] < 10 or f < 0.25:
                cells.append('   -    '); continue
            r = (Y[i] - B[i]) / (L * f * th['Z'] * dO[i]); e = np.sqrt(Y[i] + B[i]) / (L * f * th['Z'] * dO[i])
            row[(a, b)] = (Y[i], B[i], r, e, th)
            cells.append(f'{r:5.3f}({int(round(1000*e)):2d})' if e < .1 else f'{r:5.2f}({int(round(100*e)):2d})')
        res[(z0, z1)] = row
        P(f'{z0:3.0f}-{z1:<4.0f}' + ' '.join(cells) + f'   Ed {Ed:.2f}')
    return res


r0 = table(0, 'all events (no trigger cut)')
r28 = table(28, 'r_res >= 28 mm (trigger-safe; theta_cm and vz coupled, see dd_om/fiducial_window.txt)')

P('\nYields Y (background B) in the trigger-safe region, r_res >= 28:')
for zz, row in r28.items():
    if row:
        P(f'  vz {zz[0]:.0f}-{zz[1]:.0f}: ' + '  '.join(f'{a:.0f}-{b:.0f} {int(v[0])}({v[1]:.1f})' for (a, b), v in row.items()))

P('\nPotential dependence: R with each potential, r_res >= 28, cells with >= 10 counts, weighted mean over cells')
for k, n in POT.items():
    num = den = 0
    for row in r28.values():
        for (a, b), (Y, B, r, e, th) in row.items():
            rk = r * th['Z'] / th[k]; ek = e * th['Z'] / th[k]; num += rk / ek**2; den += 1 / ek**2
    if den:
        P(f'  {n:15s}: R = {num/den:.3f} +- {1/np.sqrt(den):.3f}')

P('\nSame, all events, theta_cm 18-24 only, vz 0-500:')
for k, n in POT.items():
    num = den = 0
    for (z0, z1), row in r0.items():
        if z1 > 500: continue
        for (a, b), (Y, B, r, e, th) in row.items():
            if a < 18 or b > 24: continue
            rk = r * th['Z'] / th[k]; ek = e * th['Z'] / th[k]; num += rk / ek**2; den += 1 / ek**2
    if den:
        P(f'  {n:15s}: R = {num/den:.3f} +- {1/np.sqrt(den):.3f}')

open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
cmap = plt.get_cmap('viridis')
for axi, (res, ttl) in zip(ax, [(r0, 'all events'), (r28, r'r$_{res}\geq$28 mm (trigger-safe)')]):
    for j, (zz, row) in enumerate(res.items()):
        if not row: continue
        x = [.5*(a + b) + (j - 4)*.12 for a, b in row]; y = [v[2] for v in row.values()]; e = [v[3] for v in row.values()]
        axi.errorbar(x, y, e, fmt='o', ms=4, capsize=2, color=cmap(j / 8), label=f'vz {zz[0]:.0f}-{zz[1]:.0f}')
    axi.axhline(1, color='k', lw=.8, ls=':')
    axi.set_xlabel(r'$\theta_{cm}$ [deg]'); axi.set_ylabel(r'R = Y / (L $\sigma_{ZPL}$ d$\Omega$)')
    axi.set_title(f'16C(d,d), {len(files)} runs, no acceptance: {ttl}', fontsize=10)
    axi.set_xlim(11, 31); axi.grid(alpha=.3); axi.legend(fontsize=7, ncol=2)
    yy = [v[2] for row in res.values() for v in row.values()]
    if yy: axi.set_ylim(0, min(max(yy) * 1.2, 3))
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=100)
