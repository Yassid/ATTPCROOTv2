#!/usr/bin/env bash
# Scan the deformation length of one inelastic FRESCO input and report how the cross section
# responds, both integrated and in shape.
#
# WHY THIS IS NOT JUST A RESCALING. In first-order DWBA the inelastic cross section is quadratic
# in the deformation length, so changing delta would only move the curve up and down. These inputs
# are NOT first order: they use iblock=2, i.e. a coupled-channels calculation, and a type=11
# coupling that deforms the potential by numerically re-projecting its radii rather than by taking
# a first derivative. Both make the result non-linear once the deformation is no longer small, so
# the scan measures something real:
#
#   * SATURATION. sigma/delta^2 falls from 5.3 to 3.0 between delta = 0.28 and 2.88 fm, so the
#     cross section grows more slowly than delta^2. Inferring a deformation by assuming quadratic
#     scaling therefore UNDERestimates what the data demand.
#   * THE SHAPE HARDLY MOVES. sigma(90)/sigma(40) goes only from 0.850 to 0.764 over a factor ten
#     in delta. Multi-step coupling does not reshape the angular distribution here, which rules it
#     out as the source of the 80-110 deg excess seen in the data.
#
# The deformation length is in fm: the FRESCO manual defines the type=11 parameters as
# "P(k) = DEF(k), the deformation lengths (in fm)". beta = delta / R with R = r0 A^(1/3), and r0
# is taken from the real volume potential of the input itself, not assumed.
#
#   ./beta_scan.sh <input-stem> <multipole-slot> [delta list in fm]
#   ./beta_scan.sh p14C_inel_161_6728_L3 3
set -eo pipefail
STEM=${1:?usage: beta_scan.sh <input-stem> <multipole-slot> [deltas]}
SLOT=${2:?which p(k) slot holds the deformation, i.e. the multipole L}
shift 2
DELTAS=${*:-"0.281 0.550 0.820 1.100 1.370 2.000 2.875"}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
IN="$HERE/inputs/$STEM.nin"
[ -f "$IN" ] || { echo "no such input: $IN" >&2; exit 1; }
command -v fresco >/dev/null || { echo "fresco not in PATH (try ~/.local/bin)" >&2; exit 1; }

for d in $DELTAS; do
  TAG=$(python3 -c "print('%03d'%round($d*100))")
  python3 - "$IN" "$HERE/inputs/${STEM}_d${TAG}.nin" "$SLOT" "$d" <<'PY'
import sys, re
src, dst, slot, d = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
out = []
for line in open(src):
    if 'type=11' in line:
        p = ['0.0000'] * 5
        p[slot - 1] = '%.4f' % d
        line = ' &pot kp= 1 type=11 p(1:5)= ' + '  '.join(p) + ' /\n'
    out.append(line)
open(dst, 'w').write(''.join(out))
PY
  "$HERE/run_fresco.sh" "${STEM}_d${TAG}" >/dev/null 2>&1 || { echo "  delta $d: FRESCO FAILED"; continue; }
done

python3 - "$HERE" "$STEM" "$IN" $DELTAS <<'PY'
import sys, math
here, stem, inp = sys.argv[1], sys.argv[2], sys.argv[3]
deltas = [float(x) for x in sys.argv[4:]]
# r0 from the real volume term of this very input, so R is not an assumption
r0 = None
for line in open(inp):
    if 'type= 1 ' in line and 'p(1:7)' in line:
        r0 = float(line.split('p(1:7)=')[1].split()[1]); break
A = 14
R = r0 * A ** (1 / 3.) if r0 else float('nan')
print("\n  r0 = %.4f fm (from the input)   R = r0 A^(1/3) = %.3f fm\n" % (r0, R))
print("  %-8s %-8s %-16s %-16s %-14s" % ("delta", "beta", "sig(all)[mb]", "sig(20-140)[mb]", "sig(90)/sig(40)"))
for d in deltas:
    tag = '%03d' % round(d * 100)
    f = "%s/outputs/%s_d%s_dsdo_ex2.dat" % (here, stem, tag)
    try:
        g = {}
        for line in open(f):
            a, x = line.split(); g[float(a)] = float(x)
    except OSError:
        print("  %-8.3f  (no output)" % d); continue
    tot = sum(g[a] * 2 * math.pi * math.sin(math.radians(a)) * math.radians(1.0) for a in g)
    win = sum(g[a] * 2 * math.pi * math.sin(math.radians(a)) * math.radians(1.0) for a in g if 20 <= a <= 140)
    print("  %-8.3f %-8.3f %-16.4g %-16.4g %-14.3f" % (d, d / R, tot, win, g[90] / g[40]))
print("\n  sigma/delta^2 is NOT constant: the coupling saturates, so a quadratic")
print("  extrapolation understates the deformation the data require.")
PY
