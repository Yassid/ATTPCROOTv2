# a1975 raw FRIB traces -> per event: trigger word (evt<N>_977: 1 = IC x mesh coincidence, 2 = downscaled IC,
# 3 = both; checked against the ic_ds scaler) and the list of IC and mesh pulses (position, amplitude).
# 12.5 MHz (80 ns/bin), 2048 bins; col 0 IC, col 1 mesh. Baseline = median of the whole trace (robust when a
# pulse sits in the first bins), 5-bin running mean, a pulse = an upward crossing of the threshold; its amplitude
# is the maximum before the trace falls back below threshold.
# Usage: python -I trace_pulses.py RUN OUTDIR
import sys, h5py, numpy as np

run, OUT = int(sys.argv[1]), sys.argv[2]
THR = {0: 300., 1: 60.}
NMAX = 8
f = h5py.File(f'/mnt/f/a1975/h5/run_{run:04d}.h5', 'r'); ev = f['frib/evt']
ids = sorted(int(k[3:-4]) for k in ev.keys() if k.endswith('_977'))
rows = []
for i in ids:
    w = float(np.asarray(ev[f'evt{i}_977'])[0])
    k = f'evt{i}_1903'
    if k not in ev:
        continue
    a = np.asarray(ev[k]).astype(np.float32)
    rec = [i, w]
    for c in (0, 1):
        s = np.convolve(a[:, c] - np.median(a[:, c]), np.ones(5) / 5, mode='same')
        up = s > THR[c]
        st = np.nonzero(up[1:] & ~up[:-1])[0] + 1
        if up[0]:
            st = np.r_[0, st]
        en = np.nonzero(~up[1:] & up[:-1])[0] + 1
        pos = np.full(NMAX, -1.); amp = np.zeros(NMAX)
        for j, b in enumerate(st[:NMAX]):
            e = en[en > b]; e = e[0] if e.size else len(s)
            pos[j] = b; amp[j] = s[b:e].max()
        rec += [len(st)] + list(pos) + list(amp)
    rows.append(rec)
cols = ['ev', 'word', 'icn'] + [f'icpos{j}' for j in range(NMAX)] + [f'icamp{j}' for j in range(NMAX)] + \
       ['mn'] + [f'mpos{j}' for j in range(NMAX)] + [f'mamp{j}' for j in range(NMAX)]
np.savez_compressed(f'{OUT}/pulses_{run:04d}.npz', cols=np.array(cols), d=np.array(rows, dtype=np.float32))
print('ok', run, len(rows), flush=True)
