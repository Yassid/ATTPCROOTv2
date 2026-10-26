#!/usr/bin/env python3
"""Reconstruction efficiency, measured through the student's own pipeline.

Matching reconstructed track DIRECTIONS against MC truth turns out to be a trap: the
digitised hit z is a DRIFT DISTANCE, so it runs opposite to the Geant4 z, and chasing that
sign through an absolute frame invites exactly the kind of spurious "frame rotation" a
Kabsch fit will happily invent.

The way out is to measure what the analysis actually produces.  theta_cm comes from the
arms' angle to the RECONSTRUCTED beam, and E_cm from the beam path length -- both are
invariant under any rotation, reflection or z offset between the simulation frame and the
pad frame.  So: run steps 2-4 unchanged on the simulation, and ask what fraction of the
reactions that happened come out the far end with the right answer.

That is the number the cross section must be divided by, and it folds in everything --
track finding, the beam/arm labelling rule, the vertex fit and the 90-degree cut -- instead
of track finding alone.

  usage: efficiency2.py <truth2.csv> <sim_tracks.csv> [out.png]
"""
import sys, os, csv, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
from lib.geometry import reconstruct
from lib.eloss import RangeTable

truth_csv, tracks_csv = sys.argv[1], sys.argv[2]
out_png = sys.argv[3] if len(sys.argv) > 3 else "efficiency2.png"
BEAM_ALIGN = 0.85
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

truth = {}
for r in csv.DictReader(open(truth_csv)):
    truth[int(r["event"])] = (float(r["Ecm"]), float(r["theta_cm"]),
                              float(r["E1"]), float(r["E2"]))

ev = {}
for r in csv.DictReader(open(tracks_csv)):
    ev.setdefault(int(r["event"]), {}).setdefault(int(r["track"]), []).append(
        (float(r["tb"]), float(r["x"]), float(r["y"]), float(r["q"])))

def fit_dir(P, w):
    w = np.asarray(w, float); w = w / w.sum()
    c = (P * w[:, None]).sum(axis=0); D = P - c
    val, vec = np.linalg.eigh((D * w[:, None]).T @ D)
    d = vec[:, np.argmax(val)]
    return c, d / np.linalg.norm(d)

MM = config.TB_NS * 1e-3 * config.V_DRIFT * 10.0
rec_of = {}
n3 = 0
n2 = 0
for e, tracks in ev.items():
    if len(tracks) < 2:
        continue
    n2 += 1
    if len(tracks) >= 3:
        n3 += 1
    segs = []
    for k, rows in tracks.items():
        a = np.asarray(rows)
        P = np.column_stack([a[:, 1], a[:, 2], (a[:, 0] - config.T0_BUCKETS) * MM])
        c, d = fit_dir(P, a[:, 3])
        segs.append(dict(rows=a, d=d, length=float(np.ptp(P @ d)), qmean=float(a[:, 3].mean())))
    # The beam is the near-axial track with the LOWEST charge per pad (a fast beam alpha
    # ionises less than the slower scattered ones).  It is allowed to be absent: with only
    # two tracks the arms alone fix the vertex and the tilt supplies the entrance.
    beamish = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
    beam = min(beamish, key=lambda s: s["qmean"]) if beamish else None
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
    if len(arms) < 2:
        continue
    parts, bnd, off = [], [], 0
    empty = np.zeros((0, 4))
    for s in (beam, arms[0], arms[1]):
        if s is None:
            parts.append(empty); bnd.append((off, off)); continue
        parts.append(s["rows"]); bnd.append((off, off + len(s["rows"]))); off += len(s["rows"])
    A = np.vstack(parts)
    def mk(i):
        m = np.zeros(len(A), bool); m[bnd[i][0]:bnd[i][1]] = True; return m
    cand = dict(tb=A[:, 0], x=A[:, 1], y=A[:, 2], q=A[:, 3],
                beam=mk(0), a1=mk(1), a2=mk(2))
    r = reconstruct(cand, config.V_DRIFT, config.T0_BUCKETS, config.TB_NS)
    if r is None:
        continue
    L1, L2 = r["length"]
    E1, E2 = float(tab.energy_of_range(L1)), float(tab.energy_of_range(L2))
    th = min(r["theta"])                       # the more forward arm
    rec_of[e] = dict(hasbeam=beam is not None, tcm=2.0 * th, sum=r["theta_sum"], path=r["path"],
                     Ecm_arms=(E1 + E2) / 2.0,
                     Ecm_vtx=float(tab.energy_after(r["path"])) / 2.0)

evs = sorted(truth)
Et = np.array([truth[e][0] for e in evs]); Tt = np.array([truth[e][1] for e in evs])
has = np.array([e in rec_of for e in evs])
sel = np.array([e in rec_of and abs(rec_of[e]["sum"] - 90.0) < config.SUM_WINDOW
                for e in evs])
print(f"  {len(evs)} simulated reactions")
print(f"  RANSAC gave >=2 tracks          : {100*n2/len(evs):5.1f} %")
print(f"  RANSAC gave >=3 tracks          : {100*n3/len(evs):5.1f} %")
nb = sum(1 for e in rec_of if rec_of[e]["hasbeam"])
print(f"    of the reconstructed, with a beam track: {100*nb/max(len(rec_of),1):5.1f} %")
print(f"  labelled + vertex reconstructed : {100*has.mean():5.1f} %")
print(f"  passing |theta1+theta2-90|<{config.SUM_WINDOW:.0f}    : {100*sel.mean():5.1f} %   <- eps_reco")

m = np.array([e for e in evs if e in rec_of])
if len(m):
    dt = np.array([rec_of[e]["tcm"] - truth[e][1] for e in m])
    de = np.array([rec_of[e]["Ecm_vtx"] - truth[e][0] for e in m])
    da = np.array([rec_of[e]["Ecm_arms"] - truth[e][0] for e in m])
    q = lambda v: (np.median(v), 0.5 * (np.percentile(v, 84) - np.percentile(v, 16)))
    print(f"\n  resolution (reconstructed - truth), median +- half 68% width:")
    print(f"    theta_cm         : {q(dt)[0]:+6.2f} +- {q(dt)[1]:5.2f} deg")
    print(f"    E_cm from vertex : {q(de)[0]:+6.3f} +- {q(de)[1]:5.3f} MeV")
    print(f"    E_cm from arms   : {q(da)[0]:+6.3f} +- {q(da)[1]:5.3f} MeV")

def curve(x, s, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    num, _ = np.histogram(x[s], bins=edges); den, _ = np.histogram(x, bins=edges)
    eff = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-9) / np.maximum(den, 1))
    return c, eff, np.where(den > 0, err, np.nan), den

Eed, Ted = np.arange(0.9, 6.01, 0.3), np.arange(0, 185, 10)
fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
for s, lab, col in ((has, "reconstructed", "#4C72B0"),
                    (sel, r"+ $90^\circ$ cut  ($\varepsilon_{reco}$)", "#C44E52")):
    c, eff, err, den = curve(Et, s, Eed)
    ax[0].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color=col, label=lab)
ax[0].set_xlabel(r"true $E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\varepsilon_{reco}$")
ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=8); ax[0].grid(alpha=.25)
ax[0].set_title("efficiency through the full pipeline")
c, eff, err, den = curve(Tt, sel, Ted)
ax[1].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color="#C44E52")
ax[1].set_xlabel(r"true $\theta_{cm}$ [deg]"); ax[1].set_ylabel(r"$\varepsilon_{reco}$")
ax[1].set_ylim(0, 1.05); ax[1].grid(alpha=.25); ax[1].set_title("vs angle")
if len(m):
    ax[2].hist(dt, bins=np.arange(-40, 40.5, 2), color="#C44E52", alpha=.85)
    ax[2].axvline(0, color="k", lw=1)
    ax[2].set_xlabel(r"$\theta_{cm}$ reconstructed $-$ true [deg]")
    ax[2].set_ylabel("events"); ax[2].set_title("angular resolution")
fig.tight_layout(); fig.savefig(out_png, dpi=140)
print(f"\n  wrote {out_png}")

# what step 7 consumes: truth E_cm / theta_cm per simulated reaction plus the pass flag
np.savez_compressed(os.path.join(os.path.dirname(out_png) or ".", "efficiency.npz"),
                    E=Et, T=Tt, sel=sel, has=has)
with open(out_png.replace(".png", ".csv"), "w") as fh:
    fh.write("Ecm,eff,err,N\n")
    c, eff, err, den = curve(Et, sel, Eed)
    for i in range(len(c)):
        fh.write(f"{c[i]:.3f},{eff[i]:.4f},{err[i]:.4f},{int(den[i])}\n")
        if den[i] > 20:
            print(f"  E_cm {c[i]:5.2f}  N {int(den[i]):4d}  eps {eff[i]:.3f} +- {err[i]:.3f}")
