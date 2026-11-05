# Numbers behind the trigger report: mesh early-fire fraction vs residual radius, split at vz 900 because a vertex next to the pad plane drifts inside the 100-bin window whatever fires the mesh.
import numpy as np, glob, os, json
D = '/mnt/f/a1975/sim_rebuild/trig/'; rows = []; names = None
for f in sorted(glob.glob(D + 'mesh/run_*.npz')):
    z = np.load(f); names = list(z['names']); rows.append(np.c_[np.full(len(z['d']), int(os.path.basename(f)[4:8])), z['d']])
A = np.vstack(rows); col = {n: i + 1 for i, n in enumerate(names)}
J = np.load(D + 'joined.npz'); t, i = J['t'], J['i']
th = np.radians(t[:, 3]); ke = t[:, 5]; vz = t[:, 4]; p = np.sqrt(ke**2 + 2 * ke * 938.272); pT = p * np.sin(th); pz = p * np.cos(th)
Eb = 183.873 - 0.011099 * vz; pb = np.sqrt(Eb**2 + 2 * Eb * 16.0147013 * 931.494); rres = (1000 - vz) * np.tan(np.arctan2(pT, pb - pz))
m = A[i, col['m60']]; ic = A[i, col['icbin']]; rel = m - ic; early = rel < 100; ex = t[:, 6]
out = {'bins': [], 'zone0': []}
for a, b in [(0,5),(5,10),(10,14),(14,17),(17,20),(20,23),(23,26),(26,30),(30,35),(35,45)]:
    w = (rres >= a) & (rres < b); r = {'a': a, 'b': b}
    for nm, s in [('lt900', vz < 900), ('ge900', vz >= 900), ('all', vz > -1e9)]:
        k = w & s; n = int(k.sum()); f = float(early[k].mean()) if n else None; r[nm] = {'n': n, 'f': f, 'e': float(np.sqrt(f * (1 - f) / n)) if n else None}
    out['bins'].append(r)
for z0, z1 in zip(range(0, 1000, 100), range(100, 1100, 100)):
    k = (vz >= z0) & (vz < z1); out['zone0'].append({'z0': z0, 'z1': z1, 'n': int(k.sum()), 'f': float(early[k].mean()) if k.sum() else None, 'r': float(np.median(rres[k])) if k.sum() else None})
c = {}
for nm, w in [('bound', (ex > -1.5) & (ex < .729)), ('ex35', (ex > 3) & (ex < 5))]:
    k = w & (vz >= 900); c[nm + '_ge900'] = int(k.sum())
# candidates for the example traces: bound, gated, opposite ends of the vz range
bd = (ex > -1.5) & (ex < .729) & (t[:, 7] > 25) & (t[:, 7] < 40)
out['ex'] = {}
for nm, s in [('front', bd & (vz > 200) & (vz < 400) & (rres > 28) & early), ('far', bd & (vz > 700) & (vz < 850) & (rres < 12) & (rel > 300) & (rel < 800))]:
    k = np.where(s)[0][:6]; out['ex'][nm] = [{'run': int(t[j, 0]), 'recoev': int(t[j, 1]), 'frib': int(A[i[j], col['ev']]), 'vz': float(vz[j]), 'rres': float(rres[j]), 'tc': float(t[j, 7]), 'ex': float(ex[j]), 'rel': float(rel[j]), 'mpk': float(A[i[j], col['mpk']])} for j in k]
json.dump(out, open(D + 'report_data.json', 'w'), indent=0)
for r in out['bins']: print(r['a'], r['b'], r['lt900'], r['ge900'])
for r in out['zone0']: print(r)
print(json.dumps(out['ex'], indent=0)[:1800])
