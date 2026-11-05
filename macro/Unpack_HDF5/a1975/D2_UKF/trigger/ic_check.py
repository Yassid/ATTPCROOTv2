import numpy as np
D = '/mnt/f/a1975/sim_rebuild/trig/'
for r in (16, 22, 31, 40, 57, 102):
    a = np.loadtxt(D + f'ic_{r}.txt'); z = np.load(D + f'mesh/run_{r:04d}.npz'); names = list(z['names']); d = z['d']; ic = d[:, names.index('icamp')]; ev = d[:, 0].astype(int); m = dict(zip(ev, ic)); o = []
    for off in (-1, 0, 1, 2):
        x = [v for i, (_, v) in enumerate(a) if i + off in m and v > 0]; y = [m[i + off] for i, (_, v) in enumerate(a) if i + off in m and v > 0]; o.append(np.corrcoef(x, y)[0, 1] if len(x) > 10 else 0)
    print(f'run {r}: IC entries {len(a)}, entry 0 has evt {int(a[0,0])}, last evt {int(a[-1,0])}; frib first {ev[0]} n {len(ev)}; corr(entry i IC, waveform IC of evt i+off) off -1 0 +1 +2:', np.round(o, 3))
