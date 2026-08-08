#!/usr/bin/env python3
"""Fit the p+14C optical potential to the a1954 g.s. angular distribution.

The comparison is made against the RESOLUTION-FOLDED prediction, not the raw DWBA:

    pred(r) = [ SUM_t P(r|t) A(t) dsigma/dOmega(t) sin(t) ] / ( A(r) sin(r) )

P(r|t) is the measured theta_cm migration matrix (truth-matched simulation, normalised per true
column, so it is conditional on reconstruction) and A is the simulated acceptance. Both are needed:
the data is acceptance-corrected, so A(r) is divided back out, while A(t) weights what can migrate
in the first place. Folding matters almost nowhere except the diffraction minimum, where it raises
the prediction by ~1.5x -- and that is exactly where the optical model is least constrained, so
fitting the RAW DWBA would quietly absorb the detector resolution into the imaginary depth.

Free parameters: V_V, r_V, a_V, W_D -- real depth/geometry plus the surface imaginary depth, which
between them control the position and depth of the minimum and the height of the secondary maximum.
Spin-orbit, Coulomb and the volume imaginary term stay at Koning-Delaroche 2003. The overall
normalisation is free and solved analytically at every evaluation (shape-only comparison: the data
carries no luminosity).

Optimiser is Nelder-Mead written out here rather than imported, so this has no dependency beyond
the standard library and the `fresco` binary.

    python3 fit_omp.py            # fit
    python3 fit_omp.py --check    # just evaluate the KD03 starting point
"""

import math, os, subprocess, sys, tempfile, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
FITDIR = os.path.join(HERE, "..", "pp", "plots", "ompfit")
SCRATCH = os.path.join(HERE, "outputs", "_fit")

# --- KD03 starting point for p + 14C at E_lab(p) = 11.581 MeV (see kd_params.py) ---
START = dict(VV=54.8068, rV=1.1357, aV=0.6757, WD=8.8795)
FIXED = dict(WV=0.9722, rD=1.3042, aD=0.5260, VSO=5.5100, rSO=0.9170, aSO=0.5900, rC=1.4778, A=14.0)
ELAB = 11.5810

TEMPLATE = """p + 14C elastic, OMP fit
NAMELIST
 &FRESCO  hcm= 0.05 rmatch= 30.000
     jtmin=  0.0 jtmax=   60.0 absend=  0.00001
     thmin=  1.00 thmax=180.00 thinc=  1.00
     iblock= 1
     chans= 1 smats= 2 xstabl= 1
     elab(1)= {elab:.4f} pel=1 exl=1 lab=1 lin=1 lex=1 /

 &PARTITION namep='Proton  ' massp= 1.00728 zp= 1 nex= 1
            namet=' 14C    ' masst=14.00324 zt= 6 qval= 0.0000 /
 &STATES jp= 0.5 ptyp= 1 ep= 0.0000 cpot= 1 jt= 0.0 ptyt= 1 et= 0.0000 /
 &partition /

 &pot kp= 1 type= 0 p(1:3)= {A:6.3f}  0.0000  {rC:6.4f} /
 &pot kp= 1 type= 1 p(1:7)= {VV:7.3f}  {rV:6.4f}  {aV:6.4f}  {WV:6.3f}  {rV:6.4f}  {aV:6.4f}  0.0000 /
 &pot kp= 1 type= 2 p(1:7)=  0.000  0.0000  0.0000  {WD:7.3f}  {rD:6.4f}  {aD:6.4f}  0.0000 /
 &pot kp= 1 type= 3 p(1:3)= {VSO:6.3f}  {rSO:6.4f}  {aSO:6.4f} /
 &pot /

 &overlap /
 &coupling /
"""


SYST = 0.0   # fractional point-to-point systematic added in quadrature (set by --syst)


def read_data():
    th, y, e = [], [], []
    for ln in open(os.path.join(FITDIR, "data.txt")):
        if ln.startswith("#"):
            continue
        p = ln.split()
        if len(p) >= 3:
            th.append(float(p[0])); y.append(float(p[1]))
            e.append(math.hypot(float(p[2]), SYST * float(p[1])))
    return th, y, e


def read_acc():
    a = {}
    for ln in open(os.path.join(FITDIR, "acceptance.txt")):
        p = ln.split()
        if len(p) >= 2:
            a[round(float(p[0]), 3)] = float(p[1])
    return a


def read_resp():
    rows = []
    for ln in open(os.path.join(FITDIR, "response.txt")):
        if ln.startswith("#"):
            continue
        p = ln.split()
        if len(p) == 3:
            rows.append((float(p[0]), float(p[1]), float(p[2])))
        elif len(p) == 2:
            pass  # the nx ny header line
    return rows


def run_fresco(par, tag="fit"):
    os.makedirs(SCRATCH, exist_ok=True)
    d = tempfile.mkdtemp(dir=SCRATCH)
    try:
        nin = os.path.join(d, "in.nin")
        allp = dict(FIXED); allp.update(par); allp["elab"] = ELAB
        open(nin, "w").write(TEMPLATE.format(**allp))
        with open(nin) as fi, open(os.path.join(d, "out"), "w") as fo:
            subprocess.run(["fresco"], stdin=fi, stdout=fo, stderr=subprocess.DEVNULL, cwd=d, timeout=120)
        xs = {}
        state = 0
        import re
        for ln in open(os.path.join(d, "out"), errors="ignore"):
            if "CROSS SECTIONS FOR OUTGOING" in ln:
                m = re.search(r"state\s*#\s*(\d+)", ln)
                state = int(m.group(1)) if m else 0
                continue
            if state == 1:
                m = re.match(r"\s*([0-9.]+)\s*deg.*X-S\s*=\s*([0-9.Ee+-]+)", ln)
                if m:
                    xs[float(m.group(1))] = float(m.group(2))
        return xs
    finally:
        shutil.rmtree(d, ignore_errors=True)


def interp(xs, a):
    ks = sorted(xs)
    if not ks:
        return None
    if a <= ks[0]:
        return xs[ks[0]]
    if a >= ks[-1]:
        return xs[ks[-1]]
    import bisect
    i = bisect.bisect_left(ks, a)
    x0, x1 = ks[i - 1], ks[i]
    f = (a - x0) / (x1 - x0)
    return xs[x0] * (1 - f) + xs[x1] * f


def fold(xs, resp, acc):
    """migrate the DWBA through the measured response, then undo A(r) as the data has been"""
    num = {}
    for t, r, p in resp:
        a_t = acc.get(round(t, 3), 0.0)
        if a_t <= 0:
            continue
        s_t = math.sin(math.radians(t))
        d = interp(xs, t)
        if d is None:
            continue
        num[r] = num.get(r, 0.0) + p * a_t * d * s_t
    out = {}
    for r, v in num.items():
        a_r = acc.get(round(r, 3), 0.0)
        s_r = math.sin(math.radians(r))
        if a_r > 0.05 and s_r > 1e-3:
            out[r] = v / a_r / s_r
    return out


def chi2(par, th, y, e, resp, acc, folded=True):
    xs = run_fresco(par)
    if not xs:
        return 1e12, None, None
    pred = fold(xs, resp, acc) if folded else xs
    p = [interp(pred, a) for a in th]
    if any(v is None or v <= 0 for v in p):
        return 1e12, None, None
    # analytic best normalisation k with k*p ~ y
    sn = sum(y[i] * p[i] / e[i] ** 2 for i in range(len(th)))
    sd = sum(p[i] * p[i] / e[i] ** 2 for i in range(len(th)))
    k = sn / sd if sd > 0 else 1.0
    c = sum(((y[i] - k * p[i]) / e[i]) ** 2 for i in range(len(th)))
    return c, k, p


def nelder_mead(f, x0, step, maxit=220, tol=1e-3):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0); p[i] += step[i]; pts.append(p)
    val = [f(p) for p in pts]
    for it in range(maxit):
        order = sorted(range(n + 1), key=lambda i: val[i])
        pts = [pts[i] for i in order]; val = [val[i] for i in order]
        if abs(val[-1] - val[0]) < tol * max(1.0, abs(val[0])):
            break
        cen = [sum(pts[i][j] for i in range(n)) / n for j in range(n)]
        xr = [cen[j] + 1.0 * (cen[j] - pts[-1][j]) for j in range(n)]
        fr = f(xr)
        if fr < val[0]:
            xe = [cen[j] + 2.0 * (cen[j] - pts[-1][j]) for j in range(n)]
            fe = f(xe)
            pts[-1], val[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < val[-2]:
            pts[-1], val[-1] = xr, fr
        else:
            xc = [cen[j] + 0.5 * (pts[-1][j] - cen[j]) for j in range(n)]
            fc = f(xc)
            if fc < val[-1]:
                pts[-1], val[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    pts[i] = [pts[0][j] + 0.5 * (pts[i][j] - pts[0][j]) for j in range(n)]
                    val[i] = f(pts[i])
    order = sorted(range(n + 1), key=lambda i: val[i])
    return pts[order[0]], val[order[0]]


def main():
    global SYST
    for a in sys.argv:
        if a.startswith("--syst="):
            SYST = float(a.split("=")[1])
    if SYST:
        print("adding a %.0f %% point-to-point systematic in quadrature" % (100 * SYST))
    th, y, e = read_data()
    acc = read_acc()
    resp = read_resp()
    print("data: %d points, theta_cm %.1f - %.1f deg" % (len(th), th[0], th[-1]))

    names = ["VV", "rV", "aV", "WD"]
    x0 = [START[k] for k in names]

    c0, k0, p0 = chi2(START, th, y, e, resp, acc)
    print("\nKD03 start : VV %.3f rV %.4f aV %.4f WD %.3f" % tuple(x0))
    print("             chi2 = %.1f  (chi2/N = %.2f, N = %d)" % (c0, c0 / len(th), len(th)))
    if "--check" in sys.argv:
        return

    nev = [0]
    def f(v):
        nev[0] += 1
        if v[0] < 20 or v[0] > 90 or v[1] < 0.9 or v[1] > 1.5 or v[2] < 0.35 or v[2] > 1.4 or v[3] < 0.5 or v[3] > 30:
            return 1e12
        c, _, _ = chi2(dict(zip(names, v)), th, y, e, resp, acc)
        return c

    best, cb = nelder_mead(f, x0, [3.0, 0.05, 0.04, 1.5])
    par = dict(zip(names, best))
    cB, kB, pB = chi2(par, th, y, e, resp, acc)
    print("\n%d FRESCO evaluations" % nev[0])
    print("best fit   : VV %.3f rV %.4f aV %.4f WD %.3f" % tuple(best))
    print("             chi2 = %.1f  (chi2/N = %.2f)" % (cB, cB / len(th)))
    print("\n  change from KD03: " + "  ".join(
        "%s %+.1f%%" % (n, 100 * (best[i] / x0[i] - 1)) for i, n in enumerate(names)))

    print("\n  theta_cm |     data      err |   KD03 folded  ratio |   fit folded  ratio")
    for i, a in enumerate(th):
        print("  %8.1f | %9.4g %8.4g | %12.4g %6.2f | %11.4g %6.2f"
              % (a, y[i], e[i], k0 * p0[i], y[i] / (k0 * p0[i]), kB * pB[i], y[i] / (kB * pB[i])))

    with open(os.path.join(FITDIR, "bestfit.txt"), "w") as fo:
        fo.write("# best-fit p+14C OMP, fitted to the folded prediction over theta_cm %.0f-%.0f\n" % (th[0], th[-1]))
        for n in names:
            fo.write("%s %.5f\n" % (n, par[n]))
        for n, v in FIXED.items():
            fo.write("%s %.5f\n" % (n, v))
        fo.write("chi2 %.3f\nndata %d\n" % (cB, len(th)))
    print("\nwrote %s/bestfit.txt" % FITDIR)


if __name__ == "__main__":
    main()
