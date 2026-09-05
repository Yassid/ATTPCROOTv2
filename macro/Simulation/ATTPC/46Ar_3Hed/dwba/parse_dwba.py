#!/usr/bin/env python3
"""Parse a PtolemyCpp TRANSFER run into 'theta_cm  dsigma/dOmega[mb/sr]'.

Two traps, both of which produce a plausible wrong file rather than an error:

  * THE TABLE IS PAGINATED. Ptolemy reprints the header and repeats angles across blocks, so a
    naive dict-update keeps the LAST occurrence. For the inelastic tables the last block is the
    analyzing power (negative numbers); for transfer the repeats are partial-wave decompositions
    of the same angle. Either way the FIRST value seen for an angle is the total cross section,
    so that is what is kept.

  * COLUMN 2 IS THE CROSS SECTION, column 3 is sigma/Rutherford. Reading column 3 gives a smooth
    plausible curve of the wrong quantity.

Rows look like:  '   0.00   2.4602      0.000000   -1.35    0.00     0.0 ...'
The angle always carries two decimals; the cross section may be fixed or E-notation.
"""
import re
import sys

ROW = re.compile(r"^\s+(\d+\.\d{2})\s+([-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?)\s")
HDR = re.compile(r"\(3He,d\)")


def main(src: str, dst: str) -> int:
    seen: dict[float, float] = {}
    started = False
    for line in open(src, errors="replace"):
        if "COMPUTATION OF CROSS SECTIONS" in line:
            started = True
            continue
        if not started:
            continue
        m = ROW.match(line)
        if not m:
            continue
        ang = float(m.group(1))
        val = float(m.group(2).replace("D", "E").replace("d", "e"))
        if ang not in seen:          # first occurrence wins -- see the pagination note above
            seen[ang] = val
    if len(seen) < 10:
        print(f"  parse FAILED on {src}: only {len(seen)} angles", file=sys.stderr)
        return 1
    with open(dst, "w") as f:
        f.write("# theta_cm[deg]  dsigma/dOmega[mb/sr]   from " + src + "\n")
        for a in sorted(seen):
            f.write(f"{a:8.2f}  {seen[a]:14.6g}\n")
    print(f"  -> {dst}  ({len(seen)} angles, peak {max(seen.values()):.4g} mb/sr)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
