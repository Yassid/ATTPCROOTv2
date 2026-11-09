# 16C(d,p)17C bound states: l=0 + l=2 mixture on the adopted fiducial points, plus stability against R and the l=2 curve.
# Usage: python -I fit_lmix_dp.py OUT_PREFIX    (writes OUT_PREFIX.txt and OUT_PREFIX.png)
import sys, io, contextlib, runpy
import numpy as np, uproot
from scipy.stats import chi2 as chi2dist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = sys.argv[1]
ENV = 0.84
with contextlib.redirect_stdout(io.StringIO()):
    G = runpy.run_path('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py')
run, TH, TB, rres = G['run'], G['TH'], G['TB'], G['rres']
L0, D32, D52 = TH[:, 1], TH[:, 0], TH[:, 2]
mid = .5 * (TB[:-1] + TB[1:])
X = uproot.open('/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root')['c'].arrays(['thl', 'vz', 'ke', 'ex', 'tcm', 'sv'], library='np')


def points(R):
    x, er, Y, B, e, Le = run(R)
    b = (X['ex'] > -1.5) & (X['ex'] < 1.0) & (X['vz'] > 200) & (X['vz'] < 900) & (X['thl'] >= 89) & (X['thl'] <= 161) & (X['ke'] >= 1.2) & (rres(X['ke'], X['thl'], X['vz']) >= R)
    sv = []
    for lo, hi in zip(TB[:-1], TB[1:]):
        c = b & (X['tcm'] >= lo) & (X['tcm'] < hi)
        sv.append(np.mean(X['sv'][c] > 0) if c.sum() > 20 else 1.0)
    sv = np.array(sv)
    use = (Y >= 3) & (TB[:-1] >= 25) & (TB[1:] <= 45)
    return x / sv / ENV, er / sv / ENV, use


def fit(y, e, curves):
    A = np.column_stack([c / e for c in curves])
    C = np.linalg.inv(A.T @ A)
    p = C @ A.T @ (y / e)
    return p, C, float(np.sum((y / e - A @ p)**2))


def fit_bounded(y, e, c0, c2):
    """Both coefficients >= 0: returns S0, S2, chi2."""
    p, _, q = fit(y, e, [c0, c2])
    if np.all(p >= 0):
        return p[0], p[1], q
    p2, _, q2 = fit(y, e, [c2])
    p0, _, q0 = fit(y, e, [c0])
    return (0., p2[0], q2) if q2 <= q0 else (p0[0], 0., q0)


def upper_s0(y, e, c0, c2, cl=0.90):
    """Upper limit on S0 >= 0 from the profile chi2, one-sided on the boundary."""
    s0best, _, qmin = fit_bounded(y, e, c0, c2)
    grid = np.arange(0, 60, .01)
    dq = chi2dist.ppf(2 * cl - 1, 1)
    for s0 in grid:
        ys = y - s0 * c0
        s2 = max(np.sum(ys * c2 / e**2) / np.sum(c2 * c2 / e**2), 0)
        if np.sum(((ys - s2 * c2) / e)**2) - qmin > dq and s0 > s0best:
            return s0
    return np.nan


out = []
P = lambda s='': (print(s), out.append(s))

y, e, use = points(28)
y, e = y[use], e[use]
c0, c32, c52 = L0[use], D32[use], D52[use]
c2 = c32 + c52
P('Adopted points r_res>=28, 25-45 deg, Spyral term, env 0.84')
P('bin      dsdo       err    | unit-S curves: l=0 (0.217 1/2+)  d3/2 (g.s.)  d5/2 (0.331)')
for i, m in enumerate(mid[use]):
    P(f'{m:5.1f} {y[i]:8.3f} {e[i]:8.3f}    | {c0[i]:10.3f} {c32[i]:10.3f} {c52[i]:10.3f}')
P()

P('Fits, ndf = 4 - npar')
for name, curves in [('l=2 sum', [c2]), ('l=0', [c0]), ('d3/2', [c32]), ('d5/2', [c52])]:
    p, C, q = fit(y, e, curves)
    P(f'{name:9s} S {p[0]:.4f}({np.sqrt(C[0,0]):.4f}) chi2 {q:.2f}/3  p={chi2dist.sf(q,3):.3f}')
p, C, q = fit(y, e, [c0, c2])
rho = C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])
P(f'l=0 + l=2 sum (free): S0 {p[0]:.4f}({np.sqrt(C[0,0]):.4f})  S2 {p[1]:.4f}({np.sqrt(C[1,1]):.4f})  rho {rho:.3f}  chi2 {q:.2f}/2  p={chi2dist.sf(q,2):.3f}')
b = fit_bounded(y, e, c0, c2)
P(f'l=0 + l=2 sum (S>=0): S0 {b[0]:.4f}  S2 {b[1]:.4f}  chi2 {b[2]:.2f}')
P(f'S0 upper limit, 90% CL (profile, S0>=0): {upper_s0(y, e, c0, c2):.3f}')
p3, C3, q3 = fit(y, e, [c0, c32, c52])
P(f'l=0 + d3/2 + d5/2 (free, ndf 1): {p3[0]:.3f}({np.sqrt(C3[0,0]):.3f}) {p3[1]:.3f}({np.sqrt(C3[1,1]):.3f}) {p3[2]:.3f}({np.sqrt(C3[2,2]):.3f}) chi2 {q3:.2f}/1')
P()

plb = np.array([.03, .64, .62])  # g.s. 3/2+, 1/2+, 5/2+
th = plb[0] * c32 + plb[1] * c0 + plb[2] * c52
f = np.sum(y * th / e**2) / np.sum(th * th / e**2)
P(f'PLB C2S mix (.03 .64 .62) with one free scale: scale {f:.3f}({1/np.sqrt(np.sum(th*th/e**2)):.3f}) chi2 {np.sum(((y-f*th)/e)**2):.2f}/3')
lp = np.array([.03, 0, .62])
th2 = lp[0] * c32 + lp[2] * c52
f2 = np.sum(y * th2 / e**2) / np.sum(th2 * th2 / e**2)
P(f'PLB without the l=0 state, one free scale: scale {f2:.3f} chi2 {np.sum(((y-f2*th2)/e)**2):.2f}/3')
P()

P('Stability against R, free l=0 + l=2 sum')
for R in (24, 26, 28, 30, 34):
    yr, er, ur = points(R)
    if ur.sum() < 3:
        P(f'R>={R}: {ur.sum()} bins, skipped')
        continue
    pr, Cr, qr = fit(yr[ur], er[ur], [L0[ur], (D32 + D52)[ur]])
    P(f'R>={R}: bins {int(ur.sum())}  S0 {pr[0]:.3f}({np.sqrt(Cr[0,0]):.3f})  S2 {pr[1]:.3f}({np.sqrt(Cr[1,1]):.3f})  chi2 {qr:.2f}/{int(ur.sum())-2}')

open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

m = mid[use]
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(m, y, e, fmt='D', color='k', capsize=2, label=r'adopted, r$_{res}\geq$28 mm')
lo = TB[:-1][use]
tt = np.r_[lo, lo[-1] + 5]
def st(v, **k): ax.step(tt, np.r_[v, v[-1]], where='post', **k)
s2 = fit(y, e, [c2])
st(p[0] * c0 + p[1] * c2, color='tab:red', label=f'l=0 + l=2: S0 {p[0]:.3f}({np.sqrt(C[0,0]):.3f}), S2 {p[1]:.3f}({np.sqrt(C[1,1]):.3f}), chi2 {q:.1f}/2')
st(s2[0][0] * c2, color='tab:orange', ls=':', label=f'l=2 only: {s2[0][0]:.3f}, chi2 {s2[2]:.1f}/3')
st(f * th, color='tab:blue', ls='--', label='PLB C2S mix, scale free')
ax.set_yscale('log')
ax.set_xlim(24, 46)
ax.set_ylim(np.min(y - e) * .6, np.max(y + e) * 1.8)
ax.grid(alpha=.3, which='both')
ax.legend(fontsize=8)
ax.set_xlabel(r'$\theta_{cm}$ [deg]')
ax.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]')
ax.set_title(r'$^{16}$C(d,p)$^{17}$C bound states: l=0 + l=2 mixture (stat. errors)', fontsize=10)
plt.tight_layout()
plt.savefig(OUT + '.png', dpi=100)
