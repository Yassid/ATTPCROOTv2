#!/usr/bin/env python3
"""Parse PtolemyCpp text output into two-column theta_cm / dsigma-dOmega files.

Ptolemy prints TWO different tables and they must not be confused:

  elastic   "ANGLE ... SIGMA ..."  -> column 4 is dsigma/dOmega in the c.m.
  inelastic "COMPUTATION OF CROSS SECTIONS" then a header starting "ANGLE"
            -> column 2 is the reaction dsigma/dOmega in the c.m.

The inelastic table is PAGINATED into several blocks, and the last block in the file is the
analyzing power (negative numbers), not a cross section -- so the first value seen for each angle
is kept and later ones ignored.  Getting this wrong silently yields an analyzing power plotted as
a cross section.
"""
import sys, os

def elastic(path):
    d, on = {}, False
    for ln in open(path):
        if 'ANGLE' in ln and 'SIGMA' in ln:
            on = True; continue
        if on:
            p = ln.split()
            if len(p) >= 4:
                try: d[round(float(p[0]), 1)] = float(p[3])
                except ValueError: pass
    return d

def inelastic(path):
    d, on, seen = {}, False, 0
    for ln in open(path):
        if 'COMPUTATION OF CROSS SECTIONS' in ln:
            seen += 1; on = False; continue
        if seen and ln.strip().startswith('ANGLE'):
            on = True; continue
        if on:
            p = ln.split()
            if len(p) >= 6:
                try:
                    t = round(float(p[0]), 1)
                    if t not in d: d[t] = float(p[1])   # first block wins
                except ValueError: pass
    return d

here = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
n = 0
for f in sorted(os.listdir(os.path.join(here, 'outputs'))):
    if not f.endswith('.out'): continue
    src = os.path.join(here, 'outputs', f)
    txt = open(src).read()
    d = inelastic(src) if 'COMPUTATION OF CROSS SECTIONS' in txt else elastic(src)
    d = {t: v for t, v in d.items() if v > 0}
    if not d:
        print(f"  {f}: no cross section found"); continue
    with open(os.path.join(here, 'dat', f[:-4] + '.dat'), 'w') as o:
        for t in sorted(d): o.write(f"{t} {d[t]}\n")
    n += 1
print(f"  parsed {n} outputs into dat/")
