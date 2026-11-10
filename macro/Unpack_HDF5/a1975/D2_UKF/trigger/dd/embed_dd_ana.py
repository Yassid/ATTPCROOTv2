# Survival of simulated elastic deuterons embedded in real D2 events (embed_dd.C outputs), 2026-11-10.
# Survival = mode-2 overlays that recover the mode-1 reference track / all mode-2 overlays, per theta_cm and vz of
# the reference, for all hosts and split by the host IC (carbon 1000-1350, oxygen > 1600).
# The scale test needs the CARBON-gated data, but the host species only changes the environment (beam-track size).
# Usage: python -I embed_dd_ana.py "EMBED_GLOB" OUT_PREFIX
import sys, glob
import numpy as np, uproot

files = sorted(glob.glob(sys.argv[1])); OUT = sys.argv[2]
A = [uproot.open(f)['e'].arrays(library='np') for f in files]
e = {k: np.concatenate([a[k] for a in A]) for k in A[0]}
out = []
P = lambda s='': (print(s), out.append(s))
m1 = e['mode'] == 1; m2 = e['mode'] == 2
P(f'Elastic-deuteron embedding (embed_dd_ana.py): {len(files)} files ' + ' '.join(f.rsplit("/", 1)[1] for f in files))
P(f'mode 1: {m1.sum()} sim deuterons, reference found for {np.mean(e["ref"][m1] > 0):.3f}; mode 2: {m2.sum()} overlays')
hosts = {'all': m2, 'C host': m2 & (e['icD'] >= 1000) & (e['icD'] <= 1350), 'O host': m2 & (e['icD'] > 1600),
         'low-IC host': m2 & (e['icD'] < 1000)}


def surv(m):
    n = m.sum()
    if n < 20:
        return '    -       '
    f = np.mean(e['found'][m] > 0)
    return f'{f:.3f}({np.sqrt(f*(1-f)/n)*1000:3.0f})'


P('\nSurvival by host: ' + '  '.join(f'{k} {surv(v)} n {v.sum()}' for k, v in hosts.items()))
P('\nBy theta_cm (reference), columns: ' + '  '.join(hosts))
for lo in range(12, 32, 2):
    s = (e['tc0'] >= lo) & (e['tc0'] < lo + 2)
    P(f'  {lo:2d}-{lo+2:<2d} ' + '  '.join(surv(v & s) for v in hosts.values()))
P('\nBy vz (reference), theta_cm 16-22, columns: ' + '  '.join(hosts))
for z0 in range(100, 1000, 100):
    s = (e['vz0'] >= z0) & (e['vz0'] < z0 + 100) & (e['tc0'] >= 16) & (e['tc0'] < 22)
    P(f'  {z0:3d}-{z0+100:<4d} ' + '  '.join(surv(v & s) for v in hosts.values()))
w = m2 & (e['tc0'] >= 16) & (e['tc0'] < 22) & (e['vz0'] >= 250) & (e['vz0'] < 900)
P(f'\nPLATEAU (theta_cm 16-22, vz 250-900): all hosts {surv(w)}, carbon hosts {surv(w & hosts["C host"])}')
P('By host size (data hits), plateau: ' + '  '.join(f'{a}-{b}: {surv(w & (e["nhD"] >= a) & (e["nhD"] < b))}' for a, b in [(0, 50), (50, 150), (150, 300), (300, 600), (600, 5000)]))
open(OUT + '.txt', 'w').write('\n'.join(out) + '\n')
