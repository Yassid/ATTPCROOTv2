# Soft trigger model: p(pass) = e0 + (1-e0) * Phi((ln P - ln T)/w). e0 floor: triggers not caused by this track (data has triggers with no mesh pulse). Weighted sim counts per cell, F free per angle. Bound map on theta_cm, Ex 3-5 map on theta_lab.
import numpy as np
from scipy.special import ndtr
D = '/mnt/f/a1975/sim_rebuild/trig/'; SG = [10, 20, 30, 40, 60, 80]
T = np.loadtxt(D + 'tracks.txt'); T = T[T[:, 2] == 1]; dth, dvz, dex, dtc = T[:, 3], T[:, 4], T[:, 6], T[:, 7]
S = np.loadtxt(D + 'sim_tracks.txt'); qd = {}
for s in (0, 1):
    for r in np.loadtxt(D + f'qhits/q_s{s}.txt'): qd[(s, int(r[1]))] = r[2:]
P = np.array([qd[(int(a), int(b))] for a, b in S[:, :2]]); sam, sth, svz, sex, stc = S[:, 0], S[:, 2], S[:, 3], S[:, 5], S[:, 6]
ZB = np.array([200, 300, 400, 500, 600, 700, 800, 900.]); TC = np.array([10, 20, 25, 30, 35, 40.]); TL = np.array([95, 101, 110, 120, 135.])
def counts(z, a, ab, w=None):
    iz = np.digitize(z, ZB) - 1; ia = np.digitize(a, ab) - 1; ok = (iz >= 0) & (iz < len(ZB) - 1) & (ia >= 0) & (ia < len(ab) - 1); h = np.zeros((len(ZB) - 1, len(ab) - 1)); np.add.at(h, (iz[ok], ia[ok]), 1. if w is None else w[ok]); return h
db = (dex > -1.5) & (dex < 0.729); du = (dex > 3) & (dex < 5); DB = counts(dvz[db], dtc[db], TC); DU = counts(dvz[du], dth[du], TL); b0 = sam == 0; u0 = sam == 1
def dev(d, m): d = d.ravel(); m = np.maximum(m.ravel(), 1e-6); return 2 * np.sum(m - d + np.where(d > 0, d * np.log(np.where(d > 0, d, 1) / m), 0))
def fit(d, s, mask=None):
    if mask is None: mask = np.ones_like(d, bool)
    F = np.array([d[:, j][mask[:, j]].sum() / max(s[:, j][mask[:, j]].sum(), 1e-9) for j in range(d.shape[1])]); m = s * F; return dev(d[mask], m[mask]), F, m
SB = counts(svz[b0], stc[b0], TC); SU = counts(svz[u0], sth[u0], TL)
print('no trigger: bound', round(fit(DB, SB)[0], 1), 'ndf', DB.size - 5, ' Ex3-5', round(fit(DU, SU)[0], 1), 'ndf', DU.size - 4)
mb2 = np.ones_like(DB, bool); mb2[:, 0] = False
res = []
for k, sg in enumerate(SG):
    x = np.log(P[:, k]); lo, hi = np.percentile(x, [2, 95])
    for lt in np.linspace(lo, hi, 28):
        for w_ in (0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8):
            ph = ndtr((x - lt) / w_)
            for e0 in (0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3):
                p = e0 + (1 - e0) * ph; a = fit(DB, counts(svz[b0], stc[b0], TC, p[b0]))[0]; b = fit(DU, counts(svz[u0], sth[u0], TL, p[u0]))[0]; res.append((a + b, a, b, k, np.exp(lt), w_, e0))
res.sort(key=lambda r: r[0])
for k, sg in enumerate(SG):
    r = next(r for r in res if r[3] == k); print(f'sigma {sg:3d} TB: dev total {r[0]:6.1f} (ndf {DB.size + DU.size - 9 - 3}) bound {r[1]:6.1f} Ex3-5 {r[2]:6.1f}  T {r[4]:8.0f} w {r[5]:.2f} floor {r[6]:.2f}')
r = res[0]; k = r[3]; x = np.log(P[:, k]); ph = ndtr((x - np.log(r[4])) / r[5]); p = r[6] + (1 - r[6]) * ph
_, Fb, mb = fit(DB, counts(svz[b0], stc[b0], TC, p[b0])); _, Fu, mu = fit(DU, counts(svz[u0], sth[u0], TL, p[u0])); print(f'\nBEST sigma {SG[k]} T {r[4]:.0f} w {r[5]} floor {r[6]}')
print('bound data(model), rows vz, columns theta_cm 10-20 20-25 25-30 30-35 35-40')
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{DB[i,j]:4.0f}({mb[i,j]:5.1f})' for j in range(DB.shape[1])))
print('Ex3-5 data(model), columns theta_lab 95-101 101-110 110-120 120-135')
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{DU[i,j]:4.0f}({mu[i,j]:5.1f})' for j in range(DU.shape[1])))
pc = counts(svz[b0], stc[b0], TC, p[b0]) / np.maximum(SB, 1); print('modelled trigger efficiency, bound'); 
for i in range(len(ZB) - 1): print(f' {ZB[i]:.0f}-{ZB[i+1]:.0f}', ' '.join(f'{pc[i,j]:6.3f}' for j in range(DB.shape[1])))
print('F(theta_cm) relative to 35-40, i.e. trigger-corrected shape of counts per sim event:', np.round(Fb / Fb[-1], 3), ' l=2 for 20-25,25-30,30-35: 4.9 3.0 1.7')
# stability: each map alone, and the data/sim-consistency control of a random variable in place of P
for lab, mk in [('bound alone', 0), ('Ex3-5 alone', 1)]:
    rr = min(res, key=lambda q: q[1 + mk]); print(lab, ': best dev', round(rr[1 + mk], 1), 'sigma', SG[rr[3]], 'T', round(rr[4]), 'w', rr[5], 'floor', rr[6])
rng = np.random.default_rng(1); x = rng.permutation(np.log(P[:, 3])); best = 1e9
for lt in np.linspace(*np.percentile(x, [2, 95]), 20):
    for w_ in (0.1, 0.3, 0.6):
        for e0 in (0, 0.1, 0.3):
            pz = e0 + (1 - e0) * ndtr((x - lt) / w_); d = fit(DB, counts(svz[b0], stc[b0], TC, pz[b0]))[0] + fit(DU, counts(svz[u0], sth[u0], TL, pz[u0]))[0]; best = min(best, d)
print('control, P shuffled among events: best total dev', round(best, 1))
np.savez(D + 'best_model2.npz', sigma=SG[k], T=r[4], w=r[5], floor=r[6])
