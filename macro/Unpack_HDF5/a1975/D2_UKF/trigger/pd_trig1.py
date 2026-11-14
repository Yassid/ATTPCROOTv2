# (p,d) under the mesh trigger. H2 runs 106-189, adopted selection (run_pd_chain.sh): Ebeam 183.5, theta-(360/4864.9)(KE-0.1-15.5), KE-0.1, chi2<5, vz 0-500, IC 950-1350, EXSHIFT -0.6309.
# 1) alignment reco entry <-> FRIB evt, per run, from the IC amplitude correlation; 2) early-fire vs r_res; 3) yield per mm vs vz per level x theta_cm, out to vz 1000.
import numpy as np, glob, os
from scipy.stats import chi2 as C2
D = '/mnt/f/a1975/sim_rebuild/trig/'; U = 931.49401; M1, M2, M3, M4 = 16.0147013*U, 1.00727646688*U, 2.0135532*U, 15.0105993*U
def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))
def kine(Kp, th, Ke):
    Et1 = Kp + M1; Et3 = Ke + M3; Et4 = Et1 + M2 - Et3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    a = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    a = np.where(a > 0, a, np.nan); m4 = np.sqrt(a); ex = m4 - M4; t = M2*M2 + m4*m4 - 2*M2*Et4
    arg = (s*s + s*(2*t - M1*M1 - M2*M2 - M3*M3 - m4*m4) + (M1*M1 - M2*M2)*(M3*M3 - m4*m4))/(om(s, M1*M1, M2*M2)*om(s, M3*M3, m4*m4))
    return ex, np.degrees(np.pi - np.arccos(np.clip(arg, -1, 1)))
F = {}; names = None
for f in sorted(glob.glob(D + 'mesh_h2/run_*.npz')):
    z = np.load(f); names = list(z['names']); F[int(os.path.basename(f)[4:8])] = z['d']
col = {n: i for i, n in enumerate(names)}
T = np.vstack([np.loadtxt(f, ndmin=2) for f in sorted(glob.glob(D + 'pd/tracks_*.txt')) if os.path.getsize(f)]); print('dumped deuteron rows', len(T), 'runs with mesh', len(F))
run, ent, ke, th, vz, c2, ic = T.T
# alignment per run: corr(ic_cache, icamp_frib) at offsets
off = {}; print('run  n  corr@-1  0  +1  +2')
for r in np.unique(run).astype(int):
    if r not in F: continue
    d = F[r]; idx = {int(e): i for i, e in enumerate(d[:, col['ev']])}; w = run == r; cs = []
    for o in (-1, 0, 1, 2):
        j = np.array([idx.get(int(e) + o, -1) for e in ent[w]]); k = (j >= 0) & (ic[w] > 0); cs.append(np.corrcoef(ic[w][k], d[j[k], col['icamp']])[0, 1] if k.sum() > 20 else np.nan)
    off[r] = (-1, 0, 1, 2)[int(np.nanargmax(cs))] if np.isfinite(cs).any() else None
    if r % 10 == 0 or max(cs) < .5: print(r, w.sum(), ' '.join(f'{c:5.2f}' for c in cs), '<- LOW' if max(cs) < .5 else '')
print('offset histogram', {o: sum(1 for v in off.values() if v == o) for o in (-1, 0, 1, 2, None)})
idx = np.full(len(T), -1); bad = set(r for r, o in off.items() if o is None)
for r, o in off.items():
    if o is None: continue
    d = F[r]; m = {int(e): i for i, e in enumerate(d[:, col['ev']])}; w = np.nonzero(run == r)[0]; idx[w] = [m.get(int(e) + o, -1) for e in ent[w]]
rr = run.astype(int); J = idx >= 0; RI = np.zeros(len(T), int); A = np.zeros((len(T), len(names)))
for r in F:
    w = (rr == r) & J; A[w] = F[r][idx[w]]
KE = ke - 0.1; thc = th - (360/4864.9)*(KE - 15.5); ex, tcm = kine(183.5, np.radians(thc), KE); ex = ex - 0.6309
base = (c2 < 5) & (ic >= 950) & (ic <= 1350) & (KE > 0) & (KE < 30) & np.isfinite(ex); sel = base & (vz >= 0) & (vz <= 500)
print('adopted selection', sel.sum(), ' joined', (sel & J).sum())
p = np.sqrt(KE**2 + 2*KE*M3); t = np.radians(thc); pb = np.sqrt(183.5**2 + 2*183.5*M1); rres = (1000 - vz)*np.tan(np.arctan2(p*np.sin(t), pb - p*np.cos(t)))
m60 = A[:, col['m60']]; rel = np.where(m60 >= 0, m60 - A[:, col['icbin']], 1e6); early = (rel < 100) & J; nocr = (m60 < 0) & J; pk = A[:, col['mpk']]
np.savez(D + 'pd/joined_pd.npz', T=T, idx=idx, ex=ex, tcm=tcm, rres=rres, early=early, rel=rel, pk=pk, base=base)
print('\nadopted, early fraction by r_res')
for a, b in [(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,40),(40,60),(60,100),(100,400)]:
    w = sel & J & (rres >= a) & (rres < b)
    if w.sum() > 5: print(f'  r_res {a:3d}-{b:3d}: n {w.sum():5d} early {early[w].mean():.3f} nocross {nocr[w].mean():.3f} mesh pk med {np.median(pk[w]):5.0f}')
ZB = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000.]); dz = np.diff(ZB); TB = [10, 15, 20, 25, 30, 35, 40, 45, 50]
L = {'g.s. (-0.6..0.35)': (-0.6, 0.35), '0.740 (0.45..1.2)': (0.45, 1.2), 'g.s.+0.740': (-0.6, 1.2), '3.103 (2.7..3.6)': (2.7, 3.6), '~4.7 (4.2..5.2)': (4.2, 5.2), 'cont 6-10': (6, 10)}
def test(w, lab, Z):
    n = np.histogram(vz[w], Z)[0].astype(float); d_ = np.diff(Z); lam = n.sum()*d_/d_.sum(); c = np.sum((n - lam)**2/np.maximum(lam, 1))
    x = .5*(Z[1:] + Z[:-1]); y = n/d_; s = np.sqrt(np.maximum(n, 1))/d_; W = 1/s**2; xm = np.sum(W*x)/W.sum(); ym = np.sum(W*y)/W.sum(); b = np.sum(W*(x - xm)*(y - ym))/np.sum(W*(x - xm)**2); be = 1/np.sqrt(np.sum(W*(x - xm)**2))
    ef = [early[w & (vz >= a) & (vz < b_)].mean() if (w & (vz >= a) & (vz < b_)).sum() >= 5 else np.nan for a, b_ in zip(Z[:-1], Z[1:])]
    print(f'  {lab:12s} N {int(n.sum()):5d} per100', ' '.join(f'{v*100:5.0f}' for v in n/d_), f'| chi2 {c:5.1f}/{len(n)-1} | slope/{Z[-1]-Z[0]:.0f}mm {b*(Z[-1]-Z[0])/ym:+.2f}({be*(Z[-1]-Z[0])/ym:.2f}) | early', ' '.join(f'{v:.2f}' for v in ef))
for Z, tag in [(ZB[:6], 'ADOPTED vz 0-500'), (ZB, 'FULL vz 0-1000')]:
    print(f'\n######## {tag}')
    for nm, (a, b) in L.items():
        wl = base & J & (ex > a) & (ex < b); print(f'== {nm}'); test(wl & (tcm < 50), 'tcm<50', Z)
        for p_, q in zip(TB[:-1], TB[1:]): test(wl & (tcm >= p_) & (tcm < q), f'{p_}-{q}', Z)
