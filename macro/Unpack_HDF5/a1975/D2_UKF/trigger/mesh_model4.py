# Two-path trigger model: the event triggers if the residual is outside the hole (Phi((r_res-R0)/sr)) or the proton's own pulse is high enough (Phi(ln(Pp/T)/w)), plus a small floor. Cells as mesh_model2. Then the corrected dsigma per vz range.
import numpy as np, runpy, io, contextlib
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'
with contextlib.redirect_stdout(io.StringIO()): G_ = runpy.run_path(D + 'dsigma_fid.py')
rres = G_['rres']; SG = [10, 20, 30, 40, 60, 80]
T = np.loadtxt(D + 'tracks.txt'); T = T[T[:, 2] == 1]; dth, dvz, dke, dex, dtc = T[:, 3], T[:, 4], T[:, 5], T[:, 6], T[:, 7]
S = np.loadtxt(D + 'sim_tracks.txt'); qd = {}
for s, f in ((0, 'qp_all'), (1, 'qp_s1')):
    for r in np.loadtxt(D + f'qhits/{f}.txt'): qd[(s, int(r[1]))] = r[2:]
keep = np.array([(int(a), int(b)) in qd for a, b in S[:, :2]]); print('sim rows with proton pulse', keep.sum(), 'of', len(S)); S = S[keep]; P = np.array([qd[(int(a), int(b))] for a, b in S[:, :2]]); sam, sth, svz, ske, stc = S[:, 0], S[:, 2], S[:, 3], S[:, 4], S[:, 6]; srr = rres(ske, sth, svz)
ZB = np.array([200, 300, 400, 500, 600, 700, 800, 900.]); TC = np.array([10, 20, 25, 30, 35, 40.]); TL = np.array([95, 101, 110, 120, 135.])
def counts(z, a, ab, w=None):
    iz = np.digitize(z, ZB) - 1; ia = np.digitize(a, ab) - 1; ok = (iz >= 0) & (iz < len(ZB) - 1) & (ia >= 0) & (ia < len(ab) - 1); h = np.zeros((len(ZB) - 1, len(ab) - 1)); np.add.at(h, (iz[ok], ia[ok]), 1. if w is None else w[ok]); return h
db = (dex > -1.5) & (dex < 0.729); du = (dex > 3) & (dex < 5); DB = counts(dvz[db], dtc[db], TC); DU = counts(dvz[du], dth[du], TL); b0 = sam == 0; u0 = sam == 1
def dev(d, m): d = d.ravel(); m = np.maximum(m.ravel(), 1e-6); return 2 * np.sum(m - d + np.where(d > 0, d * np.log(np.where(d > 0, d, 1) / m), 0))
def fit(d, s): F = d.sum(axis=0) / np.maximum(s.sum(axis=0), 1e-9); m = s * F; return dev(d, m), F, m
res = []
for k in (0, 1, 2, 3):
    x = np.log(P[:, k])
    for lt in np.linspace(*np.percentile(x, [5, 92]), 28):
        for w_ in (0.1, 0.15, 0.2, 0.3, 0.4, 0.6):
            pp = ndtr((x - lt) / w_)
            for R0 in (14, 18, 20, 22, 24, 26, 28, 30, 34):
                for sr in (1, 2, 3, 4, 6, 9, 13):
                    pr = ndtr((srr - R0) / sr)
                    for e0 in (0, 0.03, 0.06, 0.1, 0.13, 0.16, 0.2, 0.25):
                        p = 1 - (1 - e0) * (1 - pp) * (1 - pr); a = fit(DB, counts(svz[b0], stc[b0], TC, p[b0]))[0]; b = fit(DU, counts(svz[u0], sth[u0], TL, p[u0]))[0]; res.append((a + b, a, b, k, np.exp(lt), w_, R0, sr, e0))
res.sort(key=lambda r: r[0]); print('no trigger bound', round(fit(DB, counts(svz[b0], stc[b0], TC))[0], 1), 'Ex3-5', round(fit(DU, counts(svz[u0], sth[u0], TL))[0], 1), ' ndf', DB.size - 5, DU.size - 4)
for k in (0, 1, 2, 3):
    r = next(q for q in res if q[3] == k); print(f'sigma {SG[k]:3d}: total {r[0]:6.1f} bound {r[1]:6.1f} Ex3-5 {r[2]:6.1f}  T {r[4]:8.0f} w {r[5]} R0 {r[6]} sr {r[7]} floor {r[8]}')
r = res[0]; k = r[3]; p = 1 - (1 - r[8]) * (1 - ndtr((np.log(P[:, k]) - np.log(r[4])) / r[5])) * (1 - ndtr((srr - r[6]) / r[7])); _, Fb, mb = fit(DB, counts(svz[b0], stc[b0], TC, p[b0])); _, Fu, mu = fit(DU, counts(svz[u0], sth[u0], TL, p[u0]))
print(f'BEST sigma {SG[k]} T {r[4]:.0f} w {r[5]} R0 {r[6]} sr {r[7]} floor {r[8]}\nbound data(model), columns theta_cm 10-20 20-25 25-30 30-35 35-40')
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{DB[i,j]:4.0f}({mb[i,j]:5.1f})' for j in range(DB.shape[1])))
print('Ex3-5 data(model), columns theta_lab 95-101 101-110 110-120 120-135')
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{DU[i,j]:4.0f}({mu[i,j]:5.1f})' for j in range(DU.shape[1])))
pc = counts(svz[b0], stc[b0], TC, p[b0]) / np.maximum(counts(svz[b0], stc[b0], TC), 1); print('modelled trigger efficiency, bound')
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{pc[i,j]:6.3f}' for j in range(DB.shape[1])))
print('F rel 35-40', np.round(Fb / Fb[-1], 3)); np.savez(D + 'best_model4.npz', r=np.array(r)); print('top 5:'); [print(' ', [round(float(v), 2) for v in q]) for q in res[:5]]
