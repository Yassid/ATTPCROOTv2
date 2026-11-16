# 2026-11-16: what the (d,p) calibration does (Ex vs theta_cm, bound ridge) + bound AD uncalibrated vs calibrated (+adopted).
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("K = np.loadtxt")[0])
T = '/mnt/f/a1975/sim_rebuild/trig_hdbtau/'
d = np.loadtxt(T + 'gf_data.txt'); CAL = (184.25, 11.0, 0.18, -3, -0.05, -0.6947)
def ex_tc(cal):
    ke, th, vz, c2 = d.T.copy(); eb = Eb(vz)
    if cal: ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); eb = CAL[0] - CAL[1] / 1000 * (vz - 500)
    ex, tc = kine(eb, np.radians(th), ke); ex = ex + (CAL[5] if cal else 0)
    ok = (c2 <= 13.5) & (th >= 89) & (th <= 161) & (ke >= 1.2) & np.isfinite(ex) & (vz >= -100) & (vz <= 940); return ex[ok], tc[ok]
def ridge(ex, tc):
    out = []
    for a in range(10, 45, 5):
        x = ex[(tc >= a) & (tc < a + 5) & (ex > -1.2) & (ex < 1.2)]; m = np.median(x)
        for _ in range(20): s = x[np.abs(x - m) < 0.45]; m = s.mean()
        out.append((a + 2.5, m, s.std() / np.sqrt(len(s))))
    return np.array(out)
def csv(p):
    L = [l for l in open(p) if not l.startswith('#')]; return np.genfromtxt(L, delimiter=',', names=True)
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
for j, (cal, lab) in enumerate(((False, 'NO calibration'), (True, 'CALIBRATED'))):
    ex, tc = ex_tc(cal); a = ax[0, j]
    a.hist2d(tc, ex, bins=[np.arange(5, 61, 1), np.arange(-2, 6, 0.05)], cmap='viridis', cmin=1)
    R = ridge(ex, tc); a.errorbar(R[:, 0], R[:, 1], R[:, 2], fmt='o', color='r', ms=5, capsize=2, label='bound peak centroid')
    a.axhline(0, color='w', ls='--', lw=1)
    for e in (2.763, 3.661, 4.231): a.axhline(e, color='w', ls=':', lw=.8)
    a.set_xlim(5, 60); a.set_ylim(-2, 6); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'$E_x$ [MeV]')
    a.set_title(f'{lab}: bound centroid {R[:,1].min():+.2f} to {R[:,1].max():+.2f} MeV'); a.legend(loc="upper right", fontsize=8)
    ins = a.inset_axes([0.62, 0.08, 0.35, 0.27]); ins.errorbar(R[:, 0], R[:, 1], R[:, 2], fmt='o-', color='r', ms=3); ins.axhline(0, color='k', lw=.6)
    ins.set_ylim(-0.35, 0.3); ins.set_title('centroid vs angle', fontsize=7); ins.tick_params(labelsize=7); ins.patch.set_alpha(.85)
A = csv('/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/dp_bound_angdist_ADOPTED.csv')
N = csv(T + 'bound_hdbtau_nocal.csv'); C = csv(T + 'bound_hdbtau_cal.csv')
# DWBA l=2 sum shape scaled by each S, bin-averaged
G = [np.loadtxt('/home/yassid/fair_install/16C_dp/fresco/outputs/dp17C_bound_adwa_state%d.dat' % i, comments=['#', '&', '@']) for i in (1, 3)]
tt = np.linspace(20, 50, 300); L2 = sum(np.interp(tt, g[:, 0], g[:, 1]) for g in G) * 0.84
for k, (D, lab, c, S, dx) in enumerate(((N, 'NO calibration  S = 0.295(23), $\\chi^2$ 2.29/3', 'C0', .295, -0.3), (C, 'CALIBRATED  S = 0.288(23), $\\chi^2$ 4.49/3', 'C3', .288, 0.3))):
    for a in (ax[1, 0], ax[1, 1]):
        u = D['use'] == 1
        if a is ax[1, 0] or True: a.errorbar(D['thcm_mid'][u] + dx, D['dsdo_mb_sr'][u], D['err_mb_sr'][u], fmt='o', color=c, capsize=3, label=lab)
        a.plot(tt, S * L2, color=c, lw=1, alpha=.7)
u = A['use'] == 1; ax[1, 1].errorbar(A['thcm_mid'][u], A['dsdo_mb_sr'][u], A['err_mb_sr'][u], fmt='s', mfc='none', color='k', capsize=3, label='adopted 2026-11-09  S = 0.294(24), $\\chi^2$ 0.90/3')
for a, t, lg in ((ax[1, 0], 'bound AD: no calibration vs calibrated (HDBSCAN + float tau)', False), (ax[1, 1], 'same, log scale, with the adopted points', True)):
    a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); a.set_title(t, fontsize=10); a.set_xlim(22, 48); a.grid(alpha=.3); a.legend(fontsize=8)
    if lg: a.set_yscale('log')
ax[1, 0].set_ylim(0, 12)
fig.text(0.5, 0.005, r'lines: summed $\ell$=2 ADWA (g.s. + 0.331), $\times$ each S; $r_{res}\geq$28 mm, vz 200-900, Spyral validity, env. 0.84. Dotted white: literature 2.763 / 3.661 / 4.231.', ha='center', fontsize=8)
fig.tight_layout(rect=(0, 0.02, 1, 1)); fig.savefig(T + 'fig/dp_calib_and_AD_2026-11-16.png', dpi=120); print('saved')
