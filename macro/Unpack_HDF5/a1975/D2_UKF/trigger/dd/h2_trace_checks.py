# The trace checks of trace_checks.py repeated on the H2 (proton-target) runs (2026-11-10).
#  a. trigger-word pairing per run: word(i) or word(i+1) with waveform(i), judged by the downscaled-IC events having
#     no mesh pulse;
#  b. coincidence ceiling: latest mesh onset after the IC pulse, over every coincidence trigger;
#  c. IC pulse missing at the trigger position, per trigger type;
#  d. (p,d) deuterons (adopted selection, joined_pd.npz of 2026-11-08): mesh onset vs vz and vs r_res.
# Join for (d): reco entry i = FRIB event i+1 in all H2 runs (pd_trig note); checked here with the IC amplitude.
# Usage: python -I h2_trace_checks.py PULSE_DIR OUT_PREFIX
import sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PDIR, OUT = sys.argv[1], sys.argv[2]
out = []
P = lambda s='': (print(s), out.append(s))
files = sorted(glob.glob(PDIR + '/pulses_*.npz'))
c = list(np.load(files[0])['cols'])
ci = lambda n: c.index(n)
ICP = [ci(f'icpos{j}') for j in range(8)]; ICA = [ci(f'icamp{j}') for j in range(8)]; MPOS = [ci(f'mpos{j}') for j in range(8)]

# trigger position of the IC pulse from all events
allfirst = np.concatenate([np.load(f)['d'][:, ICP[0]] for f in files[:10]])
h, b = np.histogram(allfirst[allfirst >= 0], np.arange(0, 2049, 4)); T0 = b[h.argmax()] + 2
P(f'H2 runs: {len(files)}. IC pulse at the trigger: bin {T0:.0f}')

P('\n=== a. trigger-word pairing per run (frac of word-2 events with no mesh pulse, pairing word(i+s) with waveform(i)) ===')
pair = {}
D = {}
for f in files:
    r = int(f[-8:-4]); d = np.load(f)['d']; w = d[:, ci('word')]; mn = d[:, ci('mn')]
    fr = {}
    for s in (0, 1):
        ws = np.r_[w[s:], [np.nan] * s]
        fr[s] = (np.mean(mn[ws == 2] == 0), np.mean(mn[ws == 1] == 0))
    best = max(fr, key=lambda s: fr[s][0] - fr[s][1])
    good = fr[best][0] > 0.9 and fr[best][1] < 0.05
    pair[r] = best if good else None
    if good:
        dd = d.copy(); dd[:, ci('word')] = np.r_[w[best:], [np.nan] * best]; D[r] = dd
    P(f'  run {r}: s=0 w2 noMesh {fr[0][0]:.3f} w1 noMesh {fr[0][1]:.3f} | s=1 w2 {fr[1][0]:.3f} w1 {fr[1][1]:.3f} -> ' + (f'pair +{best}' if good else 'NEITHER, left out'))
P(f'pairing: +0 in {sum(1 for v in pair.values() if v == 0)} runs, +1 in {sum(1 for v in pair.values() if v == 1)}, neither in {sum(1 for v in pair.values() if v is None)} ({[r for r, v in pair.items() if v is None]})')

A = np.concatenate(list(D.values()))
w = A[:, ci('word')]; mp = A[:, MPOS]; icp = A[:, ICP]; ica = A[:, ICA]
first_m = np.where(mp[:, 0] >= 0, mp[:, 0] - T0, np.nan)
P('\n=== b. coincidence ceiling: first mesh onset - IC trigger, word-1 events, bins of 50 ===')
hb = np.arange(-1200, 951, 50)
P('  ' + ' '.join(f'{a:.0f}:{n}' for a, n in zip(hb[:-1], np.histogram(first_m[w == 1], hb)[0])))
lastm = np.where(mp >= 0, mp - T0, -1e9).max(1)
pos = first_m[(w == 1) & (first_m > 0)]
P(f'  word 1: p99 / p99.9 / max of the first mesh onset after the IC: {np.percentile(pos, 99):.0f} / {np.percentile(pos, 99.9):.0f} / {pos.max():.0f} bins; trace end {2047 - T0:.0f}')

P('\n=== c. IC pulse at the trigger position (+-30 bins, > 300 ADC) ===')
at = np.where((icp >= T0 - 30) & (icp <= T0 + 30), ica, 0).max(1); n = A[:, ci('icn')]
for wv in (1, 2):
    m = w == wv
    P(f'  word {wv}: n {m.sum()}  no IC pulse at trigger {np.mean(at[m] <= 300):.4f}  none anywhere {np.mean(n[m] == 0):.4f}  no mesh {np.mean(A[m, ci("mn")] == 0):.4f}')

P('\n=== d. (p,d) deuterons: mesh onset vs vz and r_res ===')
z = np.load('/mnt/f/a1975/sim_rebuild/trig/pd/joined_pd.npz'); T, ex, tcm, rres, base = z['T'], z['ex'], z['tcm'], z['rres'], z['base']
rows = []
for r in np.unique(T[:, 0]).astype(int):
    if r not in D:
        continue
    dd = D[r]; idx = {int(e): k for k, e in enumerate(dd[:, 0])}
    m = np.nonzero(T[:, 0] == r)[0]
    for k in m:
        j = idx.get(int(T[k, 1]) + 1, -1)
        if j >= 0:
            rows.append((k, j, r))
K = np.array([x[0] for x in rows]); J = np.array([x[1] for x in rows]); RR = np.array([x[2] for x in rows])
Dcat = {r: D[r] for r in D}
icamp_trig = np.array([np.where((Dcat[r][j, ICP] >= T0 - 30) & (Dcat[r][j, ICP] <= T0 + 30), Dcat[r][j, ICA], 0).max() for j, r in zip(J, RR)])
P(f'  joined {len(K)} of {len(T)} tracks; corr(cache IC, trace IC at trigger) = {np.corrcoef(T[K, 6], icamp_trig)[0,1]:.3f}')
fm = np.array([(Dcat[r][j, MPOS[0]] - T0) if Dcat[r][j, MPOS[0]] >= 0 else np.nan for j, r in zip(J, RR)])
wd = np.array([Dcat[r][j, ci('word')] for j, r in zip(J, RR)])
sel = base[K] & (T[K, 4] >= 0) & (T[K, 4] <= 900)
P(f'  word-1 fraction among (p,d) tracks: {np.mean(wd[sel] == 1):.4f}')
early = fm < 150
P('  vz slice: n  early mesh (<150 bins) frac | median first onset (bins)')
for z0 in range(0, 900, 100):
    s = sel & (T[K, 4] >= z0) & (T[K, 4] < z0 + 100)
    if s.sum() > 10:
        P(f'   {z0:3d}-{z0+100:<4d} {s.sum():5d}  {np.mean(early[s]):.3f} | {np.nanmedian(fm[s]):6.0f}')
P('  r_res slice: n early frac')
for r0, r1 in [(0, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 45), (45, 200)]:
    s = sel & (rres[K] >= r0) & (rres[K] < r1)
    if s.sum() > 10:
        P(f'   {r0:3d}-{r1:<4d} {s.sum():5d}  {np.mean(early[s]):.3f}')
open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].hist(first_m[w == 1], np.arange(-1200, 951, 10), histtype='step', label='word 1 (coincidence)')
ax[0].hist(first_m[w == 2], np.arange(-1200, 951, 10), histtype='step', label='word 2 (downscaled IC)')
ax[0].axvline(2047 - T0, color='k', ls=':', lw=1, label='end of trace'); ax[0].set_yscale('log'); ax[0].legend(fontsize=7)
ax[0].set_xlabel('first mesh onset - IC trigger [80 ns bins]'); ax[0].set_title('H2 runs, every recorded event', fontsize=10)
s = sel
ax[1].plot(T[K, 4][s], fm[s], ',', alpha=.5); ax[1].set_xlabel('vz [mm]'); ax[1].set_ylabel('first mesh onset - IC [bins]')
ax[1].set_title('(p,d) deuterons (adopted selection, vz 0-900)', fontsize=10); ax[1].set_ylim(-1200, 950)
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=100)
