# 2026-11-16 summary figure: (d,p) redo on HDBSCAN + float tau vs adopted. Bound AD (adopt_fid_dp.py CSVs) + unbound spectra.
import numpy as np, uproot, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py').read().split("K = np.loadtxt")[0])
T = '/mnt/f/a1975/sim_rebuild/trig_tctau/'
def csv(p):
    L = [l for l in open(p) if not l.startswith('#')]; return np.genfromtxt(L, delimiter=',', names=True)
A = csv('/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/dp_bound_angdist_ADOPTED.csv')
N = csv(T + 'bound_hdbtau_nocal.csv'); C = csv(T + 'bound_hdbtau_cal.csv')
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for d, lab, c, dx in ((A, 'adopted 2026-11-09 (HDBSCAN fixed tau)  S .294(24)', 'k', -0.4), (N, 'HDBSCAN + float tau, no calib  S .295(23)', 'C0', 0), (C, 'HDBSCAN + float tau, calibrated  S .288(23)', 'C3', 0.4)):
    u = d['use'] == 1; ax[0].errorbar(d['thcm_mid'][u] + dx, d['dsdo_mb_sr'][u], d['err_mb_sr'][u], fmt='o', color=c, label=lab, capsize=3)
ax[0].set_xlabel(r'$\theta_{cm}$ [deg]'); ax[0].set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); ax[0].set_title(r'bound (g.s.+0.217+0.331), $r_{res}\geq$28, vz 200-900'); ax[0].legend(fontsize=8); ax[0].set_yscale('log'); ax[0].grid(alpha=.3)
CAL = (184.25, 11.0, 0.18, -3, -0.05, -0.6947)
def spec(p, cal=False):
    t = uproot.open(p)['pk'].arrays(['ke', 'theta', 'vertexz', 'chi2ndf'], library='np'); ke, th, vz, c2 = t['ke'].astype(float), t['theta'].astype(float), t['vertexz'].astype(float), t['chi2ndf']
    eb = 183.873 - 0.011099 * vz
    if cal: ke = ke + CAL[4]; th = th - CAL[2] * (ke - CAL[3]); eb = CAL[0] - CAL[1] / 1000 * (vz - 500)
    ex, tc = kine(eb, np.radians(th), ke); ex = ex + (CAL[5] if cal else 0)
    w = (vz >= 100) & (vz < 900) & (c2 <= 13.5) & (th >= 89) & (th < 113) & (ke >= 1.2) & np.isfinite(ex); return ex[w]
bins = np.arange(-1.5, 5.0001, 0.1)
for p, cal, lab, c in (('/mnt/f/a1975/caches/.explorer/exp_dp_find805_r280.root', False, 'adopted products (no calib)', 'k'), (T + 'exp_dp_tctau.root', False, 'HDBSCAN + float tau (no calib)', 'C0'), (T + 'exp_dp_tctau.root', True, 'HDBSCAN + float tau, calibrated', 'C3')):
    ax[1].hist(spec(p, cal), bins, histtype='step', color=c, lw=1.5, label=lab)
for e in (2.763, 3.661, 4.231, 4.841): ax[1].axvline(e, color='gray', ls=':', lw=1)
ax[1].set_xlabel(r'$E_x$ [MeV]'); ax[1].set_ylabel('counts / 100 keV'); ax[1].set_title(r'$\theta_{lab}$ 89-113, vz 100-900 (literature lines dotted)'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(T + 'fig/dp_redo_hdbtau_2026-11-16.png', dpi=130); print('saved')
