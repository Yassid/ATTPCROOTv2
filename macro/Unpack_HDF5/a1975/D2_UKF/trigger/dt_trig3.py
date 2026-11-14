# (d,t): (a) does the triton-fired (late) class sit on the mesh hardware edge? (b) is the vz dependence of some windows an Ex drift (constant Ebeam 184.0) rather than a loss?
import numpy as np
D = '/mnt/f/a1975/sim_rebuild/trig/dt/'; z = np.load(D + 'joined_dt.npz'); T = z['T']; ex = z['ex']; tcm = z['tcm']; early = z['early']; base = z['base']; pk = z['pk']; rel = z['rel']; vz = T[:, 5]
sel = base & (vz >= 25) & (vz <= 600) & (tcm < 45)
E = [0, 200, 260, 280, 300, 320, 340, 360, 400, 500, 700, 1000, 1500, 1700, 5000]
for lab, w in [('early (residual path)', sel & early), ('late (triton path)', sel & ~early & (rel < 1e5))]:
    h = np.histogram(pk[w], E)[0]; print(f'{lab:24s} n {w.sum():5d} mesh peak p1 p5 p16 p50', np.round(np.percentile(pk[w], [1, 5, 16, 50])), '\n   hist', ' '.join(f'{a}-{b}:{c}' for a, b, c in zip(E[:-1], E[1:], h)))
# late class, by level window, how close to the edge?
for nm, wl in {'g.s.+0.740': (ex > -1) & (ex < 1.3), '3.103': (ex > 2.7) & (ex < 3.6), 'cont 8-14': (ex > 8) & (ex < 14)}.items():
    w = sel & ~early & (rel < 1e5) & wl; print(f'late {nm:12s} n {w.sum():4d} frac with peak < 400: {np.mean(pk[w] < 400):.3f}  < 500: {np.mean(pk[w] < 500):.3f}  p5 p50 {np.round(np.percentile(pk[w], [5, 50]))}')
# (b) Ex centroid per vz slice, bound states, adopted constant Ebeam
ZB = [25, 125, 225, 325, 425, 525, 600]
print('\nmedian Ex in window by vz slice (constant Ebeam 184.0, shift +0.2828)')
for nm, (a, b) in {'0.740 (0.3-1.3)': (0.3, 1.3), '3.103 (2.6-3.7)': (2.6, 3.7), '4.657 (4.1-5.3)': (4.1, 5.3), '6.700 (6.0-7.5)': (6.0, 7.5)}.items():
    print(f'  {nm:16s}', ' '.join(f'{np.median(ex[sel & (ex > a) & (ex < b) & (vz >= z0) & (vz < z1)]):6.3f}' for z0, z1 in zip(ZB[:-1], ZB[1:])))
