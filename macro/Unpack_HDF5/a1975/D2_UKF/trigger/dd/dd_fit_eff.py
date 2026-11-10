# Genfit fit efficiency on the 16C(d,d) elastic locus defined WITHOUT the fit (2026-11-10).
# Tree "s" of dump_dd_valid.C: every Spyral result with brho, polar, valid, hasfit. The locus is built from Spyral
# alone (deuteron KE from brho, theta_lab = 180 - polar: Spyral's polar is measured in the opposite z convention,
# Ebeam at mid-target 177.5 MeV); the fraction with a genfit fit is compared, data (carbon IC) vs sim.
# Usage: python -I dd_fit_eff.py VALID_DIR OUT.txt
import sys, glob
import numpy as np, uproot

VD, OUT = sys.argv[1], sys.argv[2]
U = 931.49401; M1 = 16.0147013 * U; Md = 2.0135532 * U
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


def lds(fs):
    A = [uproot.open(f)['s'].arrays(library='np') for f in fs]
    return {k: np.concatenate([a[k] for a in A]) for k in A[0]}


P('Genfit fit efficiency on the Spyral-defined elastic locus (dd_fit_eff.py); data sideband-subtracted')
for nm, pat, isd in (('data', '/run_*.root', True), ('sim', '/sim_*.root', False)):
    d = lds(sorted(glob.glob(VD + pat)))
    p = 299.792458 * d['brho']; ke = np.sqrt(p*p + Md*Md) - Md; th = 180 - d['polar']
    ex, cm = kine(np.full_like(ke, 177.5), np.radians(th), ke)
    b = (d['valid'] > 0) & (th > 60) & (th < 89.5) & (ke > 0.3) & (ke < 10) & np.isfinite(ex)
    if isd:
        b &= (d['ic'] >= 1000) & (d['ic'] <= 1350)
    s = b & (np.abs(ex) < 2); mu = np.median(ex[s]); q = np.percentile(ex[s], [16, 84]); sg = (q[1] - q[0]) / 2
    pk = b & (np.abs(ex - mu) < 2.5 * sg); sb = b & (ex > mu - 5) & (ex < mu - 2); f = 5 * sg / 3
    P(f'{nm}: locus Ex {mu:+.3f} sigma {sg:.3f}, n {pk.sum()}')
    for lo in range(12, 28, 2):
        m = pk & (cm >= lo) & (cm < lo + 2); mb = sb & (cm >= lo) & (cm < lo + 2)
        n = m.sum() - mb.sum() * f; nf = (d['hasfit'][m] > 0).sum() - (d['hasfit'][mb] > 0).sum() * f
        P(f'   theta_cm {lo}-{lo+2}: n {m.sum():5d} bg {mb.sum()*f:6.1f}  with a genfit fit {nf/n if n > 20 else np.nan:.3f}')
open(OUT, 'w').write('\n'.join(out) + '\n')
