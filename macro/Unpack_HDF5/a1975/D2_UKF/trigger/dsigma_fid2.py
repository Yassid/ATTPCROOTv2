# Fiducial result: fits on theta_cm 25-40 (and 20-40 where the bin exists) at each cut, each state's curve, and the Pereira-Lopez spectroscopic factors applied to our FRESCO curves with no fit. Figure.
import numpy as np, runpy, io, contextlib, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()): G_ = runpy.run_path('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py')
run, L2, L0, TH, TB = G_['run'], G_['L2'], G_['L0'], G_['TH'], G_['TB']; mid = .5 * (TB[:-1] + TB[1:])
PLB = {'gs 3/2+': 0.0, '1/2+': 0.0, '5/2+': 0.0}
def fit1(x, er, m, u): S = np.sum((x * m / er**2)[u]) / np.sum((m * m / er**2)[u]); return S, np.sqrt(1 / np.sum((m * m / er**2)[u])), np.sum((((x - S * m) / er)**2)[u])
print('fits on the bins with Y>=3 and eps>.3 inside 20-40 deg (no 40-45)')
for R in (16, 20, 22, 24, 26, 28, 30):
    x, er, Y, B, e, Le = run(R); u = (Y >= 3) & (e > .3) & (TB[:-1] >= 20) & (TB[1:] <= 40); n = u.sum()
    s2 = fit1(x, er, L2, u); s0 = fit1(x, er, L0, u); s5 = fit1(x, er, TH[:, 2], u); sg = fit1(x, er, TH[:, 0], u)
    print(f' R>={R:2d} bins {n}: l=2 sum S {s2[0]:.3f}({s2[1]:.3f}) chi2 {s2[2]:5.1f}/{n-1} | l=0 S {s0[0]:5.2f}({s0[1]:.2f}) chi2 {s0[2]:5.1f}/{n-1} | d5/2 alone S {s5[0]:.3f} chi2 {s5[2]:5.1f} | d3/2 alone S {sg[0]:.3f} chi2 {sg[2]:5.1f}')
fig, ax = plt.subplots(1, 2, figsize=(13, 5)); Rs = [0, 16, 20, 22, 24, 26, 28, 30]; col = {2: 'tab:red', 3: 'tab:orange', 4: 'tab:green', 5: 'tab:blue'}
for b, lab in [(2, '20-25'), (3, '25-30'), (4, '30-35'), (5, '35-40')]:
    xs, es, rs = [], [], []
    for R in Rs:
        x, er, Y, B, e, Le = run(R)
        if Y[b] >= 3 and e[b] > .3: xs.append(x[b]); es.append(er[b]); rs.append(R)
    ax[0].errorbar(np.array(rs) + (b - 3.5) * .25, xs, es, fmt='o-', color=col[b], label=lab + r'$^\circ$', capsize=2)
ax[0].set_xlabel('cut on predicted residual radius at the pad plane [mm]'); ax[0].set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); ax[0].set_ylim(0, 13); ax[0].grid(alpha=.3); ax[0].legend(title=r'$\theta_{cm}$'); ax[0].set_title('Cross section levels off once the residual is outside the hole')
a0 = np.array([3.68, 5.85, 6.82, 3.49]); e0 = np.array([.62, .64, .76, .49]); ax[1].errorbar(mid[2:6] - .4, a0, e0, fmt='s', color='gray', mfc='white', label='adopted 2026-09-24', capsize=2)
for R, c, m in [(26, 'tab:blue', 'o'), (28, 'k', 'D')]:
    x, er, Y, B, e, Le = run(R); u = (Y >= 3) & (e > .3) & (TB[:-1] >= 20) & (TB[1:] <= 40); ax[1].errorbar(mid[u] + (R - 27) * .3, x[u], er[u], fmt=m, color=c, label=f'r_res >= {R} mm', capsize=2)
    if R == 28: S = fit1(x, er, L2, u); ax[1].step(np.r_[TB[:-1], TB[-1]], np.r_[S[0] * L2, S[0] * L2[-1]], where='post', color='tab:red', label=f'l=2 sum x {S[0]:.2f}')
ax[1].set_yscale('log'); ax[1].set_ylim(.3, 40); ax[1].set_xlim(10, 45); ax[1].set_xlabel(r'$\theta_{cm}$ [deg]'); ax[1].set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); ax[1].legend(); ax[1].grid(alpha=.3, which='both'); ax[1].set_title(r'$^{16}$C(d,p)$^{17}$C bound states, statistical errors')
plt.tight_layout(); plt.savefig('/mnt/f/a1975/sim_rebuild/trig/fiducial_result.png', dpi=90)
