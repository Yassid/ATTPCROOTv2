# d/t separation, peak amplitude (a1975 default) vs trace integral (C15d option), same tracks.
# Input: pid_width_NNNN.root from reco_tau_d2.sh (float-tau reco + pid_width_dp.C).
# Metric y = ln(sqrt(dE/dx)*Brho) ~ ln m: ladder p 0, d +.693, t +1.099, 4He +1.386 up to one constant per plane.
# Species ANCHORED, not "tallest peak = d" (that mislabelled every window, 2026-11-07): the d band position is
# taken per plane from Brho 0.30-0.50 (pure p/d region, verified on the y-vs-Brho plot), then searched within
# +-0.15 in each window; t and 4He are fitted with means seeded at d+0.405 / d+0.693 (free within +-0.15),
# p seeded at d-0.693, plus a linear background, on y in [d-0.95, d+1.05].
# R = |mu2-mu1|/(sigma1+sigma2); valley = smoothed minimum between the two means / smaller of the two peak heights.
#   python -I pid_integral_dt_d2.py            (all finished runs except the 0016 8000-event test)
import os, glob, uproot, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
fs = sorted(f for f in glob.glob('/mnt/f/a1975/reco_d2_tau/pid_width_00[0-9][0-9].root')
            if not f.endswith('0016.root') and os.path.exists(f + '.done'))
cols = ['va', 'br', 'sqa', 'sqw']
A = {c: np.concatenate([uproot.open(f)['w'].arrays(cols, library='np')[c] for f in fs]) for c in cols}
v = (A['va'] > 0) & (A['br'] > 0.05) & (A['sqa'] > 0) & (A['sqw'] > 0)
br = A['br'][v]; Y = {q: np.log(A[q][v] * br) for q in ['sqa', 'sqw']}
print(len(fs), 'runs:', ' '.join(os.path.basename(f)[10:14] for f in fs), '|', v.sum(), 'tracks')
OFF = {'p': -0.693, 'd': 0.0, 't': 0.405, '4He': 0.693}
def g(x, a, m, s): return a * np.exp(-.5 * ((x - m) / s) ** 2)
def model(x, *p): return sum(g(x, *p[3*i:3*i+3]) for i in range(4)) + p[12] + p[13] * x
def dpos(y, sel, c, w):
    h, e = np.histogram(y[sel], 200, range=(c - w, c + w)); xc = .5*(e[1:] + e[:-1]); hs = uniform_filter1d(h.astype(float), 7)
    return xc[np.argmax(hs)]
# global anchor: d band in Brho 0.30-0.50 (search around the mode of y there, which is the d band: p is ~10x weaker)
anchor = {q: dpos(Y[q], (br > 0.30) & (br < 0.50), np.median(Y[q][(br > 0.30) & (br < 0.50)]), 0.5) for q in Y}
print('d anchor (Brho 0.30-0.50):', {q: round(a, 3) for q, a in anchor.items()})
WIN = [(0.30, 0.60), (0.60, 0.90), (0.90, 1.20), (1.20, 1.50), (0.90, 1.50)]
LAB = {'sqa': 'peak amplitude', 'sqw': 'trace integral'}
fig, ax = plt.subplots(len(WIN), 2, figsize=(14, 4.2 * len(WIN))); rows = []
for i, (b0, b1) in enumerate(WIN):
    sel = (br > b0) & (br < b1)
    for j, q in enumerate(['sqa', 'sqw']):
        y = Y[q][sel]; md = dpos(y, np.ones(len(y), bool), anchor[q], 0.15); x = y - md
        h, e = np.histogram(x, 200, range=(-0.95, 1.05)); xc = .5*(e[1:] + e[:-1]); hs = uniform_filter1d(h.astype(float), 5)
        p0, lo, hi = [], [], []
        for nm, o in OFF.items():
            a0 = max(hs[np.abs(xc - o) < 0.03].max(), 1.0)
            p0 += [a0, o, 0.07]; lo += [0, o - (0.03 if nm == 'd' else 0.15), 0.02]; hi += [np.inf, o + (0.03 if nm == 'd' else 0.15), 0.25]
        p0 += [1.0, 0.0]; lo += [0, -np.inf]; hi += [np.inf, np.inf]
        try:
            p, cv = curve_fit(model, xc, h, p0=p0, bounds=(lo, hi), sigma=np.sqrt(np.maximum(h, 1)), maxfev=40000); ok = True
        except Exception as ex:
            p, ok = np.array(p0), False
        res = {nm: (p[3*k+1], p[3*k+2], p[3*k] * p[3*k+2] * np.sqrt(2*np.pi) / (e[1] - e[0])) for k, nm in enumerate(OFF)}
        R = lambda a, b: abs(res[b][0] - res[a][0]) / (res[a][1] + res[b][1])
        def valley(a, b):
            l, u = sorted([res[a][0], res[b][0]]); w = (xc > l) & (xc < u); top = lambda c: max(hs[np.abs(xc - c) < 0.03].max(), 1e-9)
            return hs[w].min() / min(top(res[a][0]), top(res[b][0])) if w.any() else np.nan
        r = dict(win=f'{b0:.2f}-{b1:.2f}', plane=LAB[q], ok=ok, n=len(y), Rpd=R('p', 'd'), Rdt=R('d', 't'), Rta=R('t', '4He'),
                 vdt=valley('d', 't'), vta=valley('t', '4He'), Nt=res['t'][2], tpos=res['t'][0], apos=res['4He'][0],
                 sd=res['d'][1], st=res['t'][1]); rows.append(r)
        b = ax[i, j]; b.stairs(h, e, color='k'); xx = np.linspace(-0.95, 1.05, 800); b.plot(xx, model(xx, *p), 'C1', lw=1)
        for k, nm in enumerate(OFF):
            b.plot(xx, g(xx, *p[3*k:3*k+3]), 'C3', lw=.8, ls='--'); b.axvline(OFF[nm], color='C0', ls=':', alpha=.6)
            b.text(p[3*k+1], max(p[3*k], 1) * 1.15, nm, ha='center', color='C3')
        b.set_title(f'Brho {b0:.2f}-{b1:.2f}  {LAB[q]}:  d/t R {r["Rdt"]:.2f} valley {r["vdt"]:.2f} | t/4He R {r["Rta"]:.2f} | p/d R {r["Rpd"]:.2f}'
                    + ('' if ok else '  FIT FAILED'), fontsize=9)
        b.set_xlabel('ln(sqrt(dE/dx)*Brho) - d   (blue dotted: ideal ln(m/m_d))'); b.set_xlim(-0.95, 1.05)
print(f"\n{'Brho':10s} {'plane':15s} {'tracks':>7s} {'p/d R':>6s} {'d/t R':>6s} {'d/t val':>7s} {'t/a R':>6s} {'t/a val':>7s} {'N_t':>6s} {'t-d':>6s} {'a-d':>6s} {'sig_d':>6s} {'sig_t':>6s}")
for r in rows:
    print(f"{r['win']:10s} {r['plane']:15s} {r['n']:7d} {r['Rpd']:6.2f} {r['Rdt']:6.2f} {r['vdt']:7.2f} {r['Rta']:6.2f} {r['vta']:7.2f} {r['Nt']:6.0f} {r['tpos']:+6.3f} {r['apos']:+6.3f} {r['sd']:6.3f} {r['st']:6.3f}" + ('' if r['ok'] else '  FIT FAILED'))
plt.tight_layout(); fo = '/mnt/f/a1975/reco_d2_tau/pid_dt_peak_vs_integral_d2_windows.png'; plt.savefig(fo, dpi=95); print('\n->', fo)
