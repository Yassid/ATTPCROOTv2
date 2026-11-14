# (p,d) fiducial reweighting: keep tracks whose residual clears r_res >= R at the pad plane, and weight each by 500 / (length of 0-500 over which its residual would clear R)
# (r_res is linear in 1000-vz at fixed kinematics). Assumes trigger efficiency 1 above R, uniform luminosity per mm and flat acceptance in vz (sim: 0.84-0.89 over 0-800 at tcm 25-50).
# Output: correction factor to the ADOPTED vz 0-500 yields, per level and theta_cm bin. Stability vs R is the test.
import numpy as np
D = '/mnt/f/a1975/sim_rebuild/trig/'; z = np.load(D + 'pd/joined_pd.npz'); T = z['T']; ex, tcm, rres, base = z['ex'], z['tcm'], z['rres'], z['base']; vz = T[:, 4]
tanr = rres/np.maximum(1000 - vz, 1)  # tan(theta_res) per track
sel = base & (vz >= 0) & (vz <= 500)
TB = [15, 20, 25, 30, 35, 40, 45, 50, 55]
for nm, (a, b) in {'g.s. (-0.6..0.35)': (-0.6, 0.35), '0.740 (0.45..1.2)': (0.45, 1.2), '3.103 (2.7..3.6)': (2.7, 3.6)}.items():
    print(f'\n== {nm}: adopted counts | correction factor (fiducial reweighted / adopted) for R = 30 35 40 45 [kept tracks]')
    wl = sel & (ex > a) & (ex < b)
    for p, q in zip(TB[:-1], TB[1:]):
        w = wl & (tcm >= p) & (tcm < q); N = w.sum(); s = f'  tcm {p}-{q}: N {N:4d} |'
        for R in (30, 35, 40, 45):
            vmax = 1000 - R/np.maximum(tanr, 1e-9); L = np.clip(vmax, 0, 500); k = w & (rres >= R) & (L > 50)
            c = np.sum(500/L[k])/max(N, 1); err = np.sqrt(np.sum((500/L[k])**2))/max(N, 1); s += f' {c:5.2f}({err:4.2f})[{k.sum():3d}]'
        print(s)
