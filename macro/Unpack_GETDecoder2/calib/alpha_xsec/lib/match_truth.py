#!/usr/bin/env python3
"""Truth-matched tracking efficiency and angular resolution.

Track multiplicity alone is a poor proxy for efficiency: RANSAC can return three segments
that are two real arms plus a fragment, or two segments that are the two arms with the beam
missed.  What the analysis actually needs is BOTH ARMS, because the vertex is built from the
arms and the beam is only a refinement.  So match every found track to a true direction and
ask which of the three the finder actually recovered.

Matching is on direction only (|cos| -- a line fit has no sense), which is also what makes
the angular resolution fall out for free.

  usage: match_truth.py <truth2.csv> <sim_tracks.csv> [out.png]
"""
import sys, csv, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import TB_NS, V_DRIFT

truth_csv, tracks_csv = sys.argv[1], sys.argv[2]
out_png = sys.argv[3] if len(sys.argv) > 3 else "match.png"
MM_PER_TB = TB_NS * 1e-3 * V_DRIFT * 10.0      # ns -> us -> cm/us -> mm
TOL_DEG   = 10.0
MIN_HITS  = 8

truth = {}
for r in csv.DictReader(open(truth_csv)):
    g = lambda *k: np.array([float(r[x]) for x in k])
    truth[int(r["event"])] = dict(Ecm=float(r["Ecm"]), tcm=float(r["theta_cm"]),
                                  th1=float(r["theta1"]), th2=float(r["theta2"]),
                                  b=g("bdx","bdy","bdz"), d1=g("d1x","d1y","d1z"),
                                  d2=g("d2x","d2y","d2z"), v=g("vx","vy","vz")*10.0)

hits = collections.defaultdict(list)
for r in csv.DictReader(open(tracks_csv)):
    hits[(int(r["event"]), int(r["track"]))].append(
        (float(r["x"]), float(r["y"]), float(r["tb"])))

def fit_dir(p):
    """principal direction of a hit cloud"""
    c = p.mean(0)
    u, s, vt = np.linalg.svd(p - c, full_matrices=False)
    return vt[0], c, s[0] / max(s[1], 1e-9)

def build(zsign):
    tr = collections.defaultdict(list)
    for (e, k), h in hits.items():
        if len(h) < MIN_HITS:
            continue
        p = np.array([(x, y, zsign * tb * MM_PER_TB) for x, y, tb in h])
        d, c, el = fit_dir(p)
        tr[e].append(dict(d=d, c=c, n=len(h), el=el, k=k))
    return tr

def score(tr):
    """mean best |cos| between truth arms and found tracks -- picks the z sign"""
    tot, m = 0.0, 0
    for e, t in truth.items():
        for key in ("d1", "d2"):
            best = max((abs(float(np.dot(f["d"], t[key]))) for f in tr.get(e, ())), default=0.0)
            tot += best; m += 1
    return tot / max(m, 1)

trp, trm = build(+1), build(-1)
sp, sm = score(trp), score(trm)
zsign = +1 if sp >= sm else -1
tr = trp if zsign > 0 else trm
print(f"  drift-axis sign: z = {'+' if zsign>0 else '-'}tb*{MM_PER_TB:.4f} mm"
      f"   (mean |cos| {max(sp,sm):.4f} vs {min(sp,sm):.4f} for the other sign)")

rows, dang = [], {"beam": [], "arm": []}
for e in sorted(truth):
    t = truth[e]
    cand = list(tr.get(e, ()))
    got, used = {}, set()
    for key in ("d1", "d2", "b"):                      # arms first, they are the backbone
        best, bi = None, None
        for i, f in enumerate(cand):
            if i in used:
                continue
            a = np.degrees(np.arccos(min(1.0, abs(float(np.dot(f["d"], t[key]))))))
            if best is None or a < best:
                best, bi = a, i
        if bi is not None and best < TOL_DEG:
            got[key] = best; used.add(bi)
            dang["beam" if key == "b" else "arm"].append(best)
    rows.append((e, t["Ecm"], t["tcm"], len(cand),
                 "d1" in got, "d2" in got, "b" in got))

E   = np.array([r[1] for r in rows]);  T = np.array([r[2] for r in rows])
a1  = np.array([r[4] for r in rows]);  a2 = np.array([r[5] for r in rows])
bm  = np.array([r[6] for r in rows])
both = a1 & a2
print(f"\n  {len(rows)} simulated reactions, matching within {TOL_DEG:.0f} deg\n")
print(f"  forward arm found      : {100*a1.mean():5.1f} %")
print(f"  backward arm found     : {100*a2.mean():5.1f} %")
print(f"  BOTH arms (the vertex) : {100*both.mean():5.1f} %   <- eps_reco")
print(f"  beam found             : {100*bm.mean():5.1f} %")
print(f"  both arms + beam       : {100*(both&bm).mean():5.1f} %")
for k in ("arm", "beam"):
    v = np.array(dang[k])
    if len(v):
        print(f"  {k:>5} angular residual: median {np.median(v):.2f} deg, "
              f"68th pct {np.percentile(v,68):.2f} deg")

def curve(x, sel, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    num, _ = np.histogram(x[sel], bins=edges); den, _ = np.histogram(x, bins=edges)
    eff = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-9) / np.maximum(den, 1))
    return c, eff, np.where(den > 0, err, np.nan), den

Eed, Ted = np.arange(0.9, 6.01, 0.3), np.arange(0, 185, 10)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for sel, lab, col in ((both, "both arms  ($\\varepsilon_{reco}$)", "#C44E52"),
                      (both & bm, "both arms + beam", "#4C72B0"),
                      (bm, "beam", "#55A868")):
    c, eff, err, den = curve(E, sel, Eed)
    ax[0].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color=col, label=lab)
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\varepsilon_{reco}$")
ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=8); ax[0].grid(alpha=.25)
ax[0].set_title("truth-matched tracking efficiency")
c, eff, err, den = curve(T, both, Ted)
ax[1].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color="#C44E52")
ax[1].set_xlabel(r"$\theta_{cm}$ [deg]"); ax[1].set_ylabel(r"$\varepsilon_{reco}$")
ax[1].set_ylim(0, 1.05); ax[1].grid(alpha=.25); ax[1].set_title("vs angle")
v = np.array(dang["arm"])
ax[2].hist(v, bins=np.arange(0, TOL_DEG + .25, .25), color="#C44E52", alpha=.85)
ax[2].set_xlabel("arm angular residual [deg]"); ax[2].set_ylabel("tracks")
ax[2].set_title(f"angular resolution (median {np.median(v):.2f} deg)")
fig.tight_layout(); fig.savefig(out_png, dpi=140)
print(f"\n  wrote {out_png}")

with open(out_png.replace(".png", ".csv"), "w") as fh:
    fh.write("Ecm,eff,err,N\n")
    c, eff, err, den = curve(E, both, Eed)
    for i in range(len(c)):
        fh.write(f"{c[i]:.3f},{eff[i]:.4f},{err[i]:.4f},{int(den[i])}\n")
        if den[i] > 20:
            print(f"  E_cm {c[i]:5.2f}  N {int(den[i]):4d}  eps {eff[i]:.3f} +- {err[i]:.3f}")
