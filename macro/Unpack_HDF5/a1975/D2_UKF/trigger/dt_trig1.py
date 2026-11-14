# (d,t) under the mesh trigger: join adopted-selection tritons to raw FRIB mesh features; early-fire fraction vs predicted residual radius at the pad plane, vz, theta_cm; and yield per vz slice per theta_cm bin.
import numpy as np, glob, os, sys
D = '/mnt/f/a1975/sim_rebuild/trig/'; U = 931.49401; M1, M2, M3, M4 = 16.0147013*U, 2.0135532*U, 3.01550072*U, 15.0105993*U
def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))
def kine(Kp, th, Ke):  # same two-body inversion as ex_dt_a1975.C kine_2b
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
# mesh features
rows = []; names = None
for f in sorted(glob.glob(D + 'mesh/run_*.npz')):
    z = np.load(f); names = list(z['names']); rows.append(np.c_[np.full(len(z['d']), int(os.path.basename(f)[4:8])), z['d']])
A = np.vstack(rows); col = {n: i + 1 for i, n in enumerate(names)}; g = lambda n: A[:, col[n]]
key = {(int(r), int(e)): i for i, (r, e) in enumerate(zip(A[:, 0], g('ev')))}
T = np.vstack([np.loadtxt(f, ndmin=2) for f in sorted(glob.glob(D + 'dt/tracks_*.txt')) if os.path.getsize(f)])
print('dumped triton rows', len(T))
run, ent, ke, th, vz, vr, c2, ic = T[:, 0], T[:, 1], T[:, 3], T[:, 4], T[:, 5], T[:, 6], T[:, 7], T[:, 8]
# adopted selection (explorer_params_dt_2026-08-22): Ebeam 184.0, chi2 < 4.5, vz 25-600, IC 900-1400, theta - 0.046*ke, shift +0.2828 (g.s. anchored)
thc = th - 0.046*ke; ex, tcm = kine(184.0, np.radians(thc), ke); ex = ex + 0.2828
base = (c2 < 4.5) & (ic >= 900) & (ic <= 1400) & np.isfinite(ex)
sel = base & (vz >= 25) & (vz <= 600)
print('adopted selection', sel.sum(), ' theta_cm<40:', (sel & (tcm < 40)).sum(), '(memory: 6157)')
idx = np.array([key.get((int(r), int(e) + (1 if r >= 31 else 0)), -1) for r, e in zip(run, ent)])
print('joined to FRIB', (idx >= 0).sum(), 'of', len(T), '; adopted joined', (sel & (idx >= 0)).sum())
# residual radius at the pad plane (vz frame of the (d,p) study: beam enters at 0, pad plane at 1000)
p = np.sqrt(ke**2 + 2*ke*M3); t = np.radians(thc); Eb = 183.871 - 0.011071*vz; pb = np.sqrt(Eb**2 + 2*Eb*M1)
rres = (1000 - vz)*np.tan(np.arctan2(p*np.sin(t), pb - p*np.cos(t)))
j = np.where(idx >= 0, idx, 0); m = A[j, col['m60']]; icb = A[j, col['icbin']]; rel = np.where(m >= 0, m - icb, 1e6); early = (rel < 100) & (idx >= 0); nocr = (m < 0) & (idx >= 0); pk = A[j, col['mpk']]
np.savez(D + 'dt/joined_dt.npz', T=T, idx=idx, ex=ex, tcm=tcm, rres=rres, early=early, rel=rel, pk=pk, sel=sel, base=base)
J = idx >= 0
def line(w, lab):
    w = w & J; return f'{lab:22s} n {w.sum():6d} early {early[w].mean() if w.sum() else -1:.3f} nocross {nocr[w].mean() if w.sum() else -1:.3f} med r_res {np.median(rres[w]) if w.sum() else -1:6.1f} med mesh pk {np.median(pk[w]) if w.sum() else -1:6.0f}'
print('\nearly fire = mesh above 60 within 100 bins of the IC (residual path). Adopted selection, ALL Ex:')
for a, b in [(0,5),(5,10),(10,14),(14,17),(17,20),(20,23),(23,26),(26,30),(30,40),(40,60),(60,100),(100,400)]:
    print(line(sel & (rres >= a) & (rres < b), f'r_res {a}-{b}'))
L = {'g.s.+0.740 (-1..1.3)': (ex > -1) & (ex < 1.3), '3.103 (2.7-3.6)': (ex > 2.7) & (ex < 3.6), '4.657 (4.2-5.2)': (ex > 4.2) & (ex < 5.2), '6.700 (6.2-7.3)': (ex > 6.2) & (ex < 7.3)}
TB = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]; ZB = [25, 100, 200, 300, 400, 500, 600]
for nm, wl in L.items():
    print(f'\n{nm}: early fraction [n] by vz (rows) x theta_cm (cols {TB})')
    for z0, z1 in zip(ZB[:-1], ZB[1:]):
        s = f'  vz {z0:3d}-{z1:3d}'
        for a, b in zip(TB[:-1], TB[1:]):
            w = sel & J & wl & (vz >= z0) & (vz < z1) & (tcm >= a) & (tcm < b); s += f' {early[w].mean() if w.sum() >= 5 else -1:5.2f}[{w.sum():3d}]'
        print(s)
    print('  median r_res by theta_cm (vz 25-600):', ' '.join(f'{np.median(rres[sel & wl & (tcm >= a) & (tcm < b)]) if (sel & wl & (tcm >= a) & (tcm < b)).sum() else -1:5.1f}' for a, b in zip(TB[:-1], TB[1:])))
# outside the adopted vz window too: does the early fraction collapse further downstream like (d,p)?
print('\nall vz (base cuts), early fraction by vz, g.s.+0.740 / all Ex:')
for z0 in range(0, 1000, 100):
    w = base & J & (vz >= z0) & (vz < z0 + 100); w2 = w & L['g.s.+0.740 (-1..1.3)']
    print(f'  vz {z0:3d}-{z0+100:3d}: all n {w.sum():5d} early {early[w].mean() if w.sum() else -1:.3f} | bound n {w2.sum():4d} early {early[w2].mean() if w2.sum() else -1:.3f} med r_res {np.median(rres[w2]) if w2.sum() else -1:5.1f}')
