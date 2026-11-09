# 16C(d,p)17C bound states: writes the fiducial result adopted 2026-11-09.
# r_res >= 28 mm, vz 200-900, theta_cm 25-45, acceptance x Spyral validity, environment loss 0.84.
# Usage: python -I adopt_fid_dp.py OUT_PREFIX    (writes OUT_PREFIX.csv and OUT_PREFIX.png)
import sys, io, contextlib, runpy
import numpy as np, uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = sys.argv[1]
R, ENV = 28, 0.84
TRIG = '/mnt/f/a1975/sim_rebuild/trig/'
with contextlib.redirect_stdout(io.StringIO()):
    G = runpy.run_path(TRIG + 'dsigma_fid.py')
run, L2, L0, TH, TB, rres = G['run'], G['L2'], G['L0'], G['TH'], G['TB'], G['rres']
mid = .5 * (TB[:-1] + TB[1:])
X = uproot.open('/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root')['c'].arrays(['thl', 'vz', 'ke', 'ex', 'tcm', 'sv'], library='np')


def fit1(x, er, m, u):
    S = np.sum((x * m / er**2)[u]) / np.sum((m * m / er**2)[u])
    return S, np.sqrt(1 / np.sum((m * m / er**2)[u])), np.sum((((x - S * m) / er)**2)[u])


x, er, Y, B, e, Le = run(R)
b = (X['ex'] > -1.5) & (X['ex'] < 1.0) & (X['vz'] > 200) & (X['vz'] < 900) & (X['thl'] >= 89) & (X['thl'] <= 161) & (X['ke'] >= 1.2) & (rres(X['ke'], X['thl'], X['vz']) >= R)
sv, nsv = [], []
for lo, hi in zip(TB[:-1], TB[1:]):
    c = b & (X['tcm'] >= lo) & (X['tcm'] < hi)
    nsv.append(int(c.sum()))
    sv.append(np.mean(X['sv'][c] > 0) if c.sum() > 20 else 1.0)
sv = np.array(sv)
use = (Y >= 3) & (TB[:-1] >= 25) & (TB[1:] <= 45)

xs, es = x / sv, er / sv            # Spyral term only
xa, ea = xs / ENV, es / ENV         # adopted
S0 = fit1(xs, es, L2, use)
Sa = fit1(xa, ea, L2, use)
Sl0 = fit1(xa, ea, L0, use)
ndf = int(use.sum()) - 1

with open(OUT + '.csv', 'w') as f:
    f.write('# ADOPTED 2026-11-09 (Yassid): 16C(d,p)17C bound states (g.s.+0.217+0.331 unresolved), trigger-safe fiducial r_res>=28 mm, vz 200-900, theta_cm 25-45 ONLY.\n')
    f.write('# Producer D2_UKF/trigger/adopt_fid_dp.py. eps includes the Spyral validity spy_valid; dsdo includes environment loss %.2f; dsdo_noenv does not. Stat errors only; absolute scale open.\n' % ENV)
    f.write(f'# l=2 sum fit on use=1: S(no env) {S0[0]:.3f}({S0[1]:.3f}), S(adopted) {Sa[0]:.3f}({Sa[1]:.3f}), chi2 {Sa[2]:.2f}/{ndf}\n')
    f.write('thcm_lo,thcm_hi,thcm_mid,Y,bg,eps_genfit,spy_valid,n_sim_spy,Leff_mb-1,dsdo_noenv_mb_sr,err_noenv_mb_sr,dsdo_mb_sr,err_mb_sr,use\n')
    for i in range(len(mid)):
        if Y[i] == 0:
            f.write(f'{TB[i]:.0f},{TB[i+1]:.0f},{mid[i]:.1f},0,0.00,nan,nan,{nsv[i]},nan,nan,nan,nan,nan,0\n')
            continue
        f.write(f'{TB[i]:.0f},{TB[i+1]:.0f},{mid[i]:.1f},{Y[i]:.0f},{B[i]:.2f},{e[i]:.4f},{sv[i]:.4f},{nsv[i]},{Le[i]:.3f},{xs[i]:.4f},{es[i]:.4f},{xa[i]:.4f},{ea[i]:.4f},{int(use[i])}\n')

print('theta_cm   Y     B    eps   spy   dsdo_noenv   dsdo_adopted')
for i in np.where(use)[0]:
    print(f'{TB[i]:.0f}-{TB[i+1]:.0f}  {Y[i]:4.0f} {B[i]:5.2f} {e[i]:.3f} {sv[i]:.3f}  {xs[i]:.2f}({es[i]:.2f})  {xa[i]:.2f}({ea[i]:.2f})')
print(f'l=2 S no env {S0[0]:.3f}({S0[1]:.3f}) chi2 {S0[2]:.2f}/{ndf}')
print(f'l=2 S adopted {Sa[0]:.3f}({Sa[1]:.3f}) chi2 {Sa[2]:.2f}/{ndf}')
print(f'l=0 S adopted {Sl0[0]:.3f}({Sl0[1]:.3f}) chi2 {Sl0[2]:.2f}/{ndf}')

plb = TH @ np.array([0.03, 0.64, 0.62])
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(mid[use], xa[use], ea[use], fmt='D', color='k', capsize=2, label=r'r$_{res}\geq$28 mm, Spyral term, environment 0.84')
tt = np.r_[TB[:-1], TB[-1]]
ax.step(tt, np.r_[Sa[0] * L2, Sa[0] * L2[-1]], where='post', color='tab:red', label=f'l=2 sum x {Sa[0]:.3f}({Sa[1]:.3f}), chi2 {Sa[2]:.1f}/{ndf}')
ax.step(tt, np.r_[plb, plb[-1]], where='post', color='tab:blue', ls='--', label='PLB 811 C2S (.03,.64,.62) on our FRESCO curves')
ax.set_yscale('log')
ax.set_xlim(22, 48)
lo = np.min((xa - ea)[use])
ax.set_ylim(max(lo, .05 * np.max(xa[use])) * .6, np.max((xa + ea)[use]) * 2.5)
ax.grid(alpha=.3, which='both')
ax.legend(fontsize=8)
ax.set_xlabel(r'$\theta_{cm}$ [deg]')
ax.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]')
ax.set_title(r'$^{16}$C(d,p)$^{17}$C bound states, adopted 2026-11-09 (stat. errors; scale open)', fontsize=10)
plt.tight_layout()
plt.savefig(OUT + '.png', dpi=100)
