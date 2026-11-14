# C15d claims dE/dx from the INTEGRATED pulse (amp * fitted tau) separates d-t 31 % better. Same metric on a1975 D2 run_0016 (reco_d2_tau, float tau):
# ln(sqrt(dEdx)*Brho) ~ ln m, peaks p/d/t/alpha, spacings vs d and the parameter-free alpha check ln 2 = 0.693 (same beta at same Brho).
import uproot, numpy as np
from scipy.optimize import curve_fit
t = uproot.open('/mnt/f/a1975/reco_d2_tau/pid_width_0016_axis.root')['w'].arrays(library='np')
v = (t['va'] > 0) & (t['br'] > 0.05) & (t['sqa'] > 0) & (t['sqw'] > 0); print('valid Spyral tracks', v.sum(), 'of', len(v))
def g(x, a, m, s): return a*np.exp(-.5*((x - m)/s)**2)
def fitpk(x, c, w=0.12):
    s = (x > c - w) & (x < c + w); h, e = np.histogram(x[s], 40); xc = .5*(e[1:] + e[:-1])
    try: p, cv = curve_fit(g, xc, h, p0=[h.max(), c, w/2], sigma=np.sqrt(np.maximum(h, 1))); return p[1], abs(p[2]), np.sqrt(cv[1, 1]), s.sum()
    except Exception: return np.nan, np.nan, np.nan, s.sum()
for lab, q in [('PEAK amplitude (production)', 'sqa'), ('WIDTH-weighted = integral', 'sqw')]:
    x = np.log(t[q][v]*t['br'][v]); h, e = np.histogram(x, 300, range=(np.percentile(x, .5), np.percentile(x, 99.5))); xc = .5*(e[1:] + e[:-1])
    d0 = xc[np.argmax(h)]  # dominant = deuteron in D2? print to check
    md, sd, ed, nd = fitpk(x, d0)
    print(f'\n{lab}: dominant peak {d0:.3f}')
    res = {}
    for nm, off in [('p', -0.578), ('d', 0.0), ('t', 0.37), ('4He', 0.693)]:
        # search the local max near the expected spacing, then fit
        w = (xc > md + off - .15) & (xc < md + off + .15); c = xc[w][np.argmax(h[w])] if w.any() else md + off
        m, s, em, n = fitpk(x, c, 0.10 if nm != 'd' else 0.12); res[nm] = (m, s); print(f'  {nm:4s} mu-d {m - md:+.3f}({em:.3f}) sigma {s:.4f} n_in_window {n}')
    (mdd, sdd), (mt, st) = res['d'], res['t']; print(f'  d-t: separation {mt - mdd:.3f}, R = sep/(s_d+s_t) = {(mt - mdd)/(sdd + st):.2f};  p-d R = {(mdd - res["p"][0])/(sdd + res["p"][1]):.2f};  4He-d {res["4He"][0] - mdd:+.3f} (parameter-free 0.693)')
