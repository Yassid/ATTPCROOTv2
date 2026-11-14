# (p,d): is the vz falloff organised by the residual radius at the pad plane (trigger) or by vz itself? Sim acceptance in theta_cm 25-50 is flat 0.84-0.89 over vz 0-800, so data density / (density at vz 0-100, same theta_cm bin) = relative efficiency.
import numpy as np
D = '/mnt/f/a1975/sim_rebuild/trig/'; z = np.load(D + 'pd/joined_pd.npz'); T = z['T']; ex, tcm, rres, base, early = z['ex'], z['tcm'], z['rres'], z['base'], z['early']; vz = T[:, 4]
w0 = base & (ex > -0.6) & (ex < 1.2)
ZB = np.arange(0, 800, 100); TB = [20, 25, 30, 35, 40, 45, 50]
rows = []
print('g.s.+0.740: relative efficiency (density / vz 0-100 density, same theta_cm) [median r_res]')
for a, b in zip(TB[:-1], TB[1:]):
    ref = (w0 & (tcm >= a) & (tcm < b) & (vz >= 0) & (vz < 100)).sum() / 100.; s = f'  tcm {a}-{b} ref {ref*100:4.0f}/100mm:'
    for z0 in ZB:
        w = w0 & (tcm >= a) & (tcm < b) & (vz >= z0) & (vz < z0 + 100); e = w.sum()/100./ref; r = np.median(rres[w]) if w.sum() else np.nan
        # expected r_res at this vz for this theta_cm bin even if no events survive: from vz 0-100 tracks scaled by lever arm
        r_exp = np.median(rres[w0 & (tcm >= a) & (tcm < b) & (vz < 100)]) * (1000 - z0 - 50) / 950.
        s += f' {e:4.2f}[{r_exp:3.0f}]'; rows.append((r_exp, e, np.sqrt(max(w.sum(), 1))/100./ref, a, z0))
    print(s)
R = np.array(rows); print('\nrelative efficiency vs EXPECTED r_res, pooled over theta_cm and vz cells:')
for lo, hi in [(0, 10), (10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 50), (50, 70)]:
    k = (R[:, 0] >= lo) & (R[:, 0] < hi)
    if k.any(): print(f'  r_res {lo:2d}-{hi:2d}: cells {k.sum():2d} eff mean {R[k,1].mean():.2f} min {R[k,1].min():.2f} max {R[k,1].max():.2f}  (theta,vz):', ' '.join(f'({int(a)},{int(z0)})={e:.2f}' for r, e, _, a, z0 in R[k]))
