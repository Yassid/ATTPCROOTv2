#!/usr/bin/env python3
"""Digitise the CRC curves out of the report's vector PDFs.

These are vector plots, so the curve coordinates are IN the file -- there is no need to
pick pixels off a raster. pdftocairo -svg gives paths under a
matrix(0, 0.1, 0.1, 0, 0, 0) transform, i.e. a 90 degree rotation:
    canvas_X = 0.1 * path_y      (the theta_cm axis)
    canvas_Y = 0.1 * path_x      (the log cross-section axis)

Calibration is taken from the plot's own furniture, not assumed:
  * x   the frame edges carry the first and last labelled ticks, read from pdftotext -bbox
  * y   log-decade spacing is fitted to the MINOR ticks (2,3,...,9 within a decade), and the
        decade itself is anchored on the 10^n label position. Five independent minor ticks
        agree on the spacing to better than 0.5 %, which is the check that the fit is right.
"""
import re, subprocess, sys, math

def load_svg(pdf, svg):
    subprocess.run(["pdftocairo", "-svg", pdf, svg], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return open(svg).read()

def canvas_of(px, py):
    return (0.1 * py, 0.1 * px)

def segments(s):
    out = []
    for mm in re.finditer(r'<path[^>]*\sd="([^"]+)"[^>]*/?>', s):
        d = mm.group(1)
        n = [float(x) for x in re.findall(r'-?\d+\.?\d*', d)]
        if len(n) != 4 or 'C' in d:
            continue
        a = canvas_of(n[0], n[1]); b = canvas_of(n[2], n[3])
        out.append((a[0], a[1], b[0], b[1]))
    return out

def polylines(s, minpts=20):
    out = []
    for mm in re.finditer(r'<path([^>]*?)\sd="([^"]+)"([^>]*)/?>', s):
        attrs = mm.group(1) + mm.group(3); d = mm.group(2)
        if 'C' in d:
            continue
        n = [float(x) for x in re.findall(r'-?\d+\.?\d*', d)]
        if len(n) < 2 * minpts:
            continue
        pts = [canvas_of(n[i], n[i + 1]) for i in range(0, len(n) - 1, 2)]
        out.append((pts, 'dasharray' in attrs))
    return out

def words(pdf):
    txt = subprocess.run(["pdftotext", "-bbox", pdf, "-"], capture_output=True, text=True).stdout
    return [(float(a), float(b), float(c), float(d), e) for a, b, c, d, e in
            re.findall(r'<word xMin="([-\d.]+)" yMin="([-\d.]+)" xMax="([-\d.]+)" yMax="([-\d.]+)">([^<]*)</word>', txt)]

def xlabels(pdf, ybot, xleft):
    """(centre_x, value) for the numeric labels BELOW the frame -- the theta axis.

    Two filters, and BOTH are needed. Numeric-and-near-the-bottom alone lets the '10' of a
    10^-1 decade label through: it sits just above the frame bottom and would be read as an
    angle of 10 deg, silently stretching the whole x scale. That label is to the LEFT of the
    frame, where an angle label never is, so requiring x >= frame-left removes it."""
    cand = [(0.5 * (a + c), float(e)) for a, b, c, d, e in words(pdf)
            if re.fullmatch(r'\d+', e) and b > ybot - 5 and 0.5 * (a + c) > xleft - 5]
    return sorted(cand)

def decades(pdf):
    """(y_centre, exponent) for each 10^n label: a '10' with its exponent word up and to the right."""
    ws = words(pdf)
    tens = [(a, b, c, d) for a, b, c, d, e in ws if e == '10']
    exps = [(a, b, c, d, e) for a, b, c, d, e in ws if re.fullmatch(r'-?\d+', e)]
    out = []
    for a, b, c, d in tens:
        best, bd = None, 1e9
        for ea, eb, ec, ed, ee in exps:
            if ea < c - 1:            # exponent must start to the right of the '10'
                continue
            dist = abs(eb - b) + abs(ea - c)
            if dist < bd and abs(eb - b) < 25:
                bd, best = dist, int(ee)
        if best is not None:
            out.append((0.5 * (b + d), best))
    return sorted(out)

def ycal(segs, pdf):
    """Return (Y0, D) with value = 10**((Y0 - Y)/D).

    With two or more decade labels the spacing is measured directly from them, which needs no
    assumption about which minor ticks are present. With only one, fall back to the minor-tick
    ladder (0.9, 0.8, ... below the decade) and report its spread as the quality check."""
    decs = decades(pdf)
    if len(decs) >= 2:
        Ds = [(decs[i + 1][0] - decs[i][0]) / (decs[i][1] - decs[i + 1][1]) for i in range(len(decs) - 1)]
        D = sum(Ds) / len(Ds)
        spread = (max(Ds) - min(Ds)) / D if len(Ds) > 1 else 0.0
        # anchor on value 1 (exponent 0) if present, else extrapolate from the first label
        y0lab, n0 = next(((y, n) for y, n in decs if n == 0), decs[0])
        Y0 = y0lab + n0 * D
        print(f"    y: {len(decs)} decade labels {[n for _, n in decs]}, D={D:.1f} "
              f"(spread {100*spread:.2f} %), value 1 at Y={Y0:.1f}")
        return Y0, D
    xl = min(min(a, c) for a, b, c, d in segs)
    ticks = sorted(set(round(b, 2) for a, b, c, d in segs
                       if abs(b - d) < 0.5 and abs(a - xl) < 1.0 and abs(c - a) < 15))
    if not decs or len(ticks) < 6:
        raise SystemExit("  cannot calibrate the y axis")
    ylab, n = decs[0]
    Y0lab = min(ticks, key=lambda t: abs(t - ylab))
    below = [t for t in ticks if t > Y0lab + 1][:5]
    Ds = [(t - Y0lab) / (-math.log10(0.9 - 0.1 * i)) for i, t in enumerate(below)]
    D = sum(Ds) / len(Ds)
    spread = (max(Ds) - min(Ds)) / D
    Y0 = Y0lab + n * D
    print(f"    y: 1 decade label (10^{n}), D={D:.1f} from {len(Ds)} minor ticks "
          f"(spread {100*spread:.2f} %), value 1 at Y={Y0:.1f}")
    if spread > 0.02:
        print("    \033[1;31mWARNING: minor ticks disagree by >2 % -- calibration suspect\033[0m")
    return Y0, D

def run(stem, pdf):
    print(f"\n  {stem}")
    s = load_svg(pdf, f"{stem}.svg")
    segs = segments(s)
    ybot = max(max(b, d) for a, b, c, d in segs)
    xleft = min(min(a, c) for a, b, c, d in segs)
    xls = xlabels(pdf, ybot, xleft)
    if len(xls) < 2:
        raise SystemExit("  could not read the x-axis labels")
    (x1, v1), (x2, v2) = xls[0], xls[-1]
    print(f"    x: {v1:.0f} deg at X={x1:.1f}, {v2:.0f} deg at X={x2:.1f}")
    Y0, D = ycal(segs, pdf)
    todeg = lambda X: v1 + (X - x1) * (v2 - v1) / (x2 - x1)
    toval = lambda Y: 10 ** ((Y0 - Y) / D)

    pls = sorted(polylines(s), key=lambda t: -len(t[0]))[:3]
    solid = [p for p, dash in pls if not dash]
    dashed = [p for p, dash in pls if dash]
    for name, group in (("full", solid), ("reduced", dashed)):
        pts = []
        for p in group:
            pts += [(todeg(X), toval(Y)) for X, Y in p]
        pts = [(a, v) for a, v in pts if v > 0]
        pts.sort()
        if not pts:
            print(f"    {name}: no curve found")
            continue
        fn = f"{stem}_{name}.dat"
        with open(fn, "w") as f:
            for a, v in pts:
                f.write(f"{a:8.3f} {v:12.6g}\n")
        print(f"    {name:8s}: {len(pts):4d} points, theta {pts[0][0]:.1f}-{pts[-1][0]:.1f}, "
              f"sigma {min(v for _, v in pts):.3g}-{max(v for _, v in pts):.3g} -> {fn}")

for stem, pdf in (("gs", "/mnt/c/Users/Yassid/Downloads/gs.pdf"),
                  ("e074", "/mnt/c/Users/Yassid/Downloads/0.74.pdf"),
                  ("e310", "/mnt/c/Users/Yassid/Downloads/3.1.pdf")):
    run(stem, pdf)
