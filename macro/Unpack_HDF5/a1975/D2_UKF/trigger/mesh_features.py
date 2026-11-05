# FRIB waveforms of one a1975 run -> per-event features of the IC (col 0) and the mesh (col 1). 12.5 MHz, 2048 bins. Baseline = median of bins 0-300.
import h5py, numpy as np, sys
run = int(sys.argv[1]); f = h5py.File(f'/mnt/f/a1975/h5/run_{run:04d}.h5', 'r'); ev = f['frib/evt']
n = sum(1 for k in f['get'].keys() if k.endswith('_data'))
names = ['ev','icbin','icamp','icn','mbase','mpk','mpkbin','minteg','m30','m60','m100','m200','m400','mlast100','msat','m0900','c2amp','c2bin']
out = []
for i in range(n + 2):
    k = f'evt{i}_1903'
    if k not in ev: continue
    a = np.asarray(ev[k]).astype(np.float32)
    ic = a[:, 0] - np.median(a[:300, 0]); m = a[:, 1]; mb = np.median(m[:300]); m = m - mb; c2 = a[:, 2] - np.median(a[:300, 2])
    win = ic[1100:1200]; icbin = 1100 + int(win.argmax()); icamp = float(win.max())
    above = ic > 300; icn = int(np.sum(above[1:] & ~above[:-1]))
    # 5-bin running mean to tame the digitiser noise before looking for crossings
    s = np.convolve(m, np.ones(5) / 5, mode='same')
    pk = float(s.max()); pkb = int(s.argmax()); integ = float(np.clip(s, 0, None).sum())
    def first(t):
        w = np.nonzero(s > t)[0]; return int(w[0]) if w.size else -1
    w = np.nonzero(s > 100)[0]; last = int(w[-1]) if w.size else -1
    plateau = int(np.sum(s > 0.98 * pk)) if pk > 50 else 0
    out.append([i, icbin, icamp, icn, mb, pk, pkb, integ, first(30), first(60), first(100), first(200), first(400), last, plateau, float(s[:900].max()), float(c2.max()), int(c2.argmax())])
np.savez_compressed(f'/mnt/f/a1975/sim_rebuild/trig/mesh/run_{run:04d}.npz', names=np.array(names), d=np.array(out, dtype=np.float32))
s = f['frib/scaler']; tot = np.zeros(11)
for k in s.keys():
    if k.endswith('_data'):
        a = np.asarray(s[k]).astype(float).ravel()
        if a.size >= 11: tot += a[:11]
np.save(f'/mnt/f/a1975/sim_rebuild/trig/mesh/scaler_{run:04d}.npy', tot); print('ok', run, len(out), flush=True)
