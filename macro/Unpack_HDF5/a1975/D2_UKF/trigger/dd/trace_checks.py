# The three trace checks for the trigger questions (2026-11-10), on the runs with pulses_NNNN.npz and dd caches:
#  1. mesh-minus-IC delay vs reconstructed vz for 16C(d,d) elastic events: hard edge = coincidence window,
#     events with a later mesh pulse but no track = readout window;
#  2. elastic events with a low IC amplitude (pile-up candidates): IC pulses elsewhere in the trace? their vz?
#  3. downscaled-IC triggers (word 2): fraction with no IC pulse at the trigger position / anywhere.
# Trigger word: runs 16-27 pair word(i+1) with waveform(i) (checked: then 97-98 % of word-2 events have no mesh pulse,
# as in run 31); run 22 fits neither pairing (58839 FRIB events for 39978 entries) and is left out.
# Join: dd cache row ev = reco tree entry; FRIB id = entry (runs <= 27) or entry + 1 (runs >= 31).
# Usage: python -I trace_checks.py CACHE_DIR PULSE_DIR OUT_PREFIX RUNS(comma)
import sys
import numpy as np, uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CDIR, PDIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
RUNS = [int(x) for x in sys.argv[4].split(',')]
U = 931.49401; M1 = 16.0147013 * U; Md = 2.0135532 * U
Eb = lambda z: 183.871 - 0.011071 * z
MMPB = 1.10424 * 10 * 0.08          # mm of drift per 80 ns bin (dv 1.10424 cm/us)
ZPP = 971.73                        # pad plane in z_true (par ZPadPlane); vz_fit ~ z_true + 19.1 (sim)
out = []
P = lambda s='': (print(s), out.append(s))


def om(x, y, z): return np.sqrt(np.maximum(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z, 0))


def ex_el(Kp, th, Ke):
    M2 = M3 = Md; M4 = M1
    Et1 = Kp + M1; Et3 = Ke + M3; s = M1*M1 + M2*M2 + 2*M2*Et1; uu = M2*M2 + M3*M3 - 2*M2*Et3
    A = (np.cos(th)*om(s, M1*M1, M2*M2)*om(uu, M2*M2, M3*M3) - (s - M1*M1 - M2*M2)*(M2*M2 + M3*M3 - uu))/(2*M2*M2) + s + uu - M2*M2
    return np.sqrt(np.where(A > 0, A, np.nan)) - M4


rows = []
PUL = {}
for r in RUNS:
    z = np.load(f'{PDIR}/pulses_{r:04d}.npz'); c = list(z['cols']); d = z['d']
    if r <= 27:  # runs 16-27: the trigger word stored with FRIB event i belongs to the waveform of event i+1
        d = d.copy(); d[:-1, c.index('word')] = d[1:, c.index('word')]; d[-1, c.index('word')] = np.nan
    PUL[r] = (c, d)
    t = uproot.open(f'{CDIR}/run_{r:04d}.root')['d'].arrays(library='np')
    fid = t['ev'].astype(int) + (1 if r >= 31 else 0)
    idx = {int(e): k for k, e in enumerate(d[:, 0])}
    k = np.array([idx.get(int(e), -1) for e in fid])
    ok = k >= 0
    rec = {kk: t[kk][ok] for kk in t}
    for kk in c:
        rec['p_' + kk] = d[k[ok], c.index(kk)]
    rows.append(rec)
D = {kk: np.concatenate([x[kk] for x in rows]) for kk in rows[0]}
ICP = np.stack([D[f'p_icpos{j}'] for j in range(8)], 1); ICA = np.stack([D[f'p_icamp{j}'] for j in range(8)], 1)
MP = np.stack([D[f'p_mpos{j}'] for j in range(8)], 1); MA = np.stack([D[f'p_mamp{j}'] for j in range(8)], 1)

# trigger position of the IC pulse: the mode of the first pulse in coincidence events
w = D['p_word']
first = ICP[:, 0]
h, b = np.histogram(first[(w == 1) & (first >= 0)], np.arange(0, 2049, 4))
T0 = b[h.argmax()] + 2
P(f'Runs {RUNS}; joined {len(w)} tracks. IC pulse at the trigger: mode of first IC pulse = bin {T0:.0f}')
intrig = (ICP >= T0 - 30) & (ICP <= T0 + 30)
amp_trig = np.where(intrig, ICA, 0).max(1)
# join check: cache ic (unpackIC, entry-aligned) against the amplitude at the trigger position
sel = (D['ic'] > 0) & (amp_trig > 0)
for r in RUNS:
    s = sel & (D['run'] == r)
    if s.sum() > 50:
        P(f'  join check run {r}: corr(cache ic, trace IC at trigger) = {np.corrcoef(D["ic"][s], amp_trig[s])[0,1]:.3f}  (n {s.sum()})')

ex = ex_el(Eb(D['vz']), np.radians(D['th']), D['ke'])
el = (D['th'] > 60) & (D['th'] < 89.5) & (D['ke'] > 0.3) & (np.abs(ex + 0.1) < 0.6)
carbon = (D['ic'] >= 1000) & (D['ic'] <= 1350)

# ---------------- check 1
P('\n=== Check 1: mesh pulse time vs vertex, elastic deuterons (theta_lab 60-89.5, |Ex+0.1|<0.6) ===')
P(f'Expected drift for a vertex at vz: (ZPP - (vz - 19.1)) / {MMPB:.4f} mm per bin, ZPP {ZPP}')
# the mesh pulse belonging to the deuteron: the mesh pulse closest to the expected drift time
exp_bins = (ZPP - (D['vz'] - 19.1)) / MMPB
dm = np.where(MP >= 0, MP - T0, np.nan)
diff = np.abs(dm - exp_bins[:, None])
jbest = np.nanargmin(np.where(np.isnan(diff), 1e9, diff), 1)
mbest = dm[np.arange(len(dm)), jbest]
mfirst = np.where(MP[:, 0] >= 0, MP[:, 0] - T0, np.nan)
ZB = np.arange(50, 951, 50.)
P('vz bin    n_el  frac word1  frac no mesh pulse | median first-mesh - IC (bins)  expected  | median matched  p90 matched | last mesh bin in trace - IC')
for z0, z1 in zip(ZB[:-1], ZB[1:]):
    s = el & carbon & (D['vz'] >= z0) & (D['vz'] < z1)
    if s.sum() < 5:
        P(f'{z0:4.0f}-{z1:<4.0f} {s.sum():6d}'); continue
    e = np.median(exp_bins[s])
    P(f'{z0:4.0f}-{z1:<4.0f} {s.sum():6d}   {np.mean(w[s]==1):.3f}      {np.mean(D["p_mn"][s]==0):.3f}        | {np.nanmedian(mfirst[s]):7.0f}   {e:7.0f}  | {np.nanmedian(mbest[s]):7.0f}  {np.nanpercentile(mbest[s],90):7.0f} | {2047-T0:.0f}')

# all coincidence events (not only reconstructed): latest mesh pulse onset relative to IC
lastm = np.where(MP >= 0, MP, -1).max(1)
P('\nAll joined tracks, word 1: distribution of (first mesh pulse - IC trigger) in bins of 50, and of the latest mesh onset:')
hb = np.arange(-1200, 951, 50)
hf = np.histogram(mfirst[(w == 1)], hb)[0]
P('  ' + ' '.join(f'{a:.0f}:{n}' for a, n in zip(hb[:-1], hf)))

# ---------------- check 2
P('\n=== Check 2: elastic events with a low IC amplitude (pile-up candidates) ===')
low = el & (D['ic'] >= 0) & (D['ic'] < 1000)
norm = el & carbon
P(f'elastic peak: carbon {norm.sum()}, low-IC (<1000) {low.sum()} ({low.sum()/(low.sum()+norm.sum()+np.sum(el&(D["ic"]>1350))):.3f} of the peak)')
for name, s in [('carbon', norm), ('low-IC', low)]:
    n_any = np.mean(D['p_icn'][s] > 0); n_trig = np.mean(amp_trig[s] > 300)
    other = np.where(~intrig & (ICP >= 0), ICP - T0, np.nan)
    P(f'  {name:7s}: n {s.sum():6d}  IC pulse at trigger position {n_trig:.3f}  any IC pulse {n_any:.3f}  median IC pulses {np.median(D["p_icn"][s]):.0f}  word1 {np.mean(w[s]==1):.3f}')
lowp = low & (amp_trig <= 300)
oth = np.where(~intrig & (ICP >= 0), ICP - T0, np.nan)[lowp]
oth = oth[np.isfinite(oth)]
if oth.size:
    P(f'  low-IC with no pulse at the trigger: other IC pulses at (bin - trigger) p10/50/90 {np.percentile(oth,[10,50,90]).round(0)}; '
      f'frac before trigger {np.mean(oth<0):.2f}; amplitude median {np.median(ICA[lowp][(~intrig[lowp])&(ICP[lowp]>=0)]):.0f}')
vb = np.arange(0, 1001, 100.)
P('  vz distribution (fraction per 100 mm), carbon vs low-IC-no-trigger-pulse:')
P('   carbon ' + ' '.join(f'{x:.3f}' for x in np.histogram(D['vz'][norm], vb)[0] / norm.sum()))
if lowp.sum():
    P('   low-IC ' + ' '.join(f'{x:.3f}' for x in np.histogram(D['vz'][lowp], vb)[0] / lowp.sum()))
    P(f'   vz outside 0-1000: carbon {np.mean((D["vz"][norm]<0)|(D["vz"][norm]>1000)):.3f}, low-IC {np.mean((D["vz"][lowp]<0)|(D["vz"][lowp]>1000)):.3f}')
    P(f'   Ex median/sigma: carbon {np.median(ex[norm]):+.3f}/{np.std(ex[norm]):.3f}, low-IC {np.median(ex[lowp]):+.3f}/{np.std(ex[lowp]):.3f}')

# ---------------- check 3: all events of the runs (not only reconstructed)
P('\n=== Check 3: downscaled-IC triggers (word 2) vs coincidence triggers (word 1), every recorded event ===')
for r in RUNS:
    c, d = PUL[r]
    ww = d[:, c.index('word')]
    icp = d[:, [c.index(f'icpos{j}') for j in range(8)]]; ica = d[:, [c.index(f'icamp{j}') for j in range(8)]]
    it = (icp >= T0 - 30) & (icp <= T0 + 30)
    at = np.where(it, ica, 0).max(1); n = d[:, c.index('icn')]; mn = d[:, c.index('mn')]
    s = []
    for wv in (1, 2, 3):
        m = ww == wv
        if m.sum():
            s.append(f'w{wv} n {m.sum():5d} noICatTrig {np.mean(at[m]<=300):.4f} noICanywhere {np.mean(n[m]==0):.4f} noMesh {np.mean(mn[m]==0):.3f}')
    P(f'run {r}: ' + ' | '.join(s))
allc = np.concatenate([PUL[r][1] for r in RUNS]); c = PUL[RUNS[0]][0]
ww = allc[:, c.index('word')]; icp = allc[:, [c.index(f'icpos{j}') for j in range(8)]]; ica = allc[:, [c.index(f'icamp{j}') for j in range(8)]]
at = np.where((icp >= T0 - 30) & (icp <= T0 + 30), ica, 0).max(1); n = allc[:, c.index('icn')]
for wv in (1, 2):
    m = ww == wv
    P(f'ALL word {wv}: n {m.sum()}  no IC pulse at trigger {np.sum(at[m]<=300)} ({np.mean(at[m]<=300):.4f})  no IC pulse anywhere {np.sum(n[m]==0)} ({np.mean(n[m]==0):.4f})')
m2 = (ww == 2) & (at <= 300) & (n > 0)
if m2.sum():
    o = np.where((icp >= 0) & ~((icp >= T0 - 30) & (icp <= T0 + 30)), icp - T0, np.nan)[m2]
    o = o[np.isfinite(o)]
    P(f'  word 2 with no pulse at trigger but a pulse elsewhere: {m2.sum()}; pulse position - trigger p10/50/90 {np.percentile(o,[10,50,90]).round(0)}')

open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
s = el & carbon
ax[0].plot(D['vz'][s], mfirst[s], ',', alpha=.4, label='first mesh pulse')
zz = np.linspace(0, 1000, 50); ax[0].plot(zz, (ZPP - (zz - 19.1)) / MMPB, 'r-', lw=1, label='drift from vertex (dv 1.104)')
ax[0].axhline(2047 - T0, color='k', ls=':', lw=1, label='end of trace')
ax[0].set_xlabel('vz [mm]'); ax[0].set_ylabel('mesh onset - IC trigger [80 ns bins]'); ax[0].set_title('elastic d (carbon): mesh onset vs vertex', fontsize=10)
ax[0].legend(fontsize=7, markerscale=10); ax[0].set_ylim(-1200, 950)
ax[1].hist(D['vz'][s], np.arange(0, 1001, 20), histtype='step', label='elastic d, carbon')
ax[1].set_xlabel('vz [mm]'); ax[1].set_title('elastic vertex distribution', fontsize=10); ax[1].legend(fontsize=7)
ax[2].hist(mfirst[w == 1], np.arange(-1200, 951, 10), histtype='step', label='word 1 (coincidence)', density=True)
ax[2].hist(mfirst[w == 2], np.arange(-1200, 951, 10), histtype='step', label='word 2 (downscaled IC)', density=True)
ax[2].set_xlabel('first mesh onset - IC trigger [bins]'); ax[2].set_title('reconstructed tracks: mesh onset', fontsize=10); ax[2].legend(fontsize=7); ax[2].set_yscale('log')
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=100)
