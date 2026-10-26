#!/usr/bin/env python3
"""Convert a RANSAC track CSV into the pickle format steps 3-6 expect.

AtSampleConsensus returns an unordered list of tracks.  The downstream analysis wants them
labelled: which is the beam, which are the two scattered arms.  Same rule as the Python
finder uses, so the two are compared on equal terms -- the beam is the track pointing
closest to the drift axis with the LOWEST mean charge per pad (a fast beam alpha ionises
less than the slower scattered ones), and the arms are the two longest of the rest.

Two tracks are enough.  The beam is a refinement, not a requirement.

  usage: ransac_to_pickle.py <ransac_XXXX.csv> <tracks_XXXX.pkl>
"""
import sys, os, csv, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config

BEAM_ALIGN = 0.85

def fit_dir(P, w):
    w = np.asarray(w, float); w = w / w.sum()
    c = (P * w[:, None]).sum(axis=0)
    D = P - c
    C = (D * w[:, None]).T @ D
    val, vec = np.linalg.eigh(C)
    d = vec[:, np.argmax(val)]
    return c, d / np.linalg.norm(d)

src, dst = sys.argv[1], sys.argv[2]
ev = {}
with open(src) as fh:
    for r in csv.DictReader(fh):
        e = int(r["event"]); k = int(r["track"])
        ev.setdefault(e, {}).setdefault(k, []).append(
            (float(r["tb"]), float(r["x"]), float(r["y"]), float(r["q"])))

cands = []
for e, tracks in ev.items():
    if len(tracks) < 2:
        continue
    segs = []
    for k, rows in tracks.items():
        a = np.asarray(rows)
        tb, x, y, q = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
        z = (tb - config.T0_BUCKETS) * config.TB_NS * 1e-3 * config.V_DRIFT * 10.0
        P = np.column_stack([x, y, z])
        c, d = fit_dir(P, q)
        segs.append(dict(rows=a, d=d, length=float(np.ptp(P @ d)), qmean=float(q.mean())))
    # The beam may legitimately be missing: when the reaction happens close to the window
    # there is almost no beam track to find.  The two arms alone fix the vertex, and the
    # nominal tilt then supplies the entrance (see geometry.reconstruct).  Measured on
    # simulation this recovers the whole E_cm > 4.7 MeV half of the acceptance, which is
    # otherwise reconstructed with ZERO efficiency.
    beamish = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
    beam = min(beamish, key=lambda s: s["qmean"]) if beamish else None
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
    if len(arms) < 2:
        continue
    # flatten into the one-array-plus-masks layout the rest of the chain uses
    parts, masks, off = [], [], 0
    empty = np.zeros((0, 4))
    for s in (beam, arms[0], arms[1]):
        rows = empty if s is None else s["rows"]
        parts.append(rows); masks.append((off, off + len(rows))); off += len(rows)
    allrows = np.vstack(parts)
    def mask(i):
        m = np.zeros(len(allrows), bool); m[masks[i][0]:masks[i][1]] = True; return m
    cands.append(dict(event=e, tb=allrows[:, 0], x=allrows[:, 1], y=allrows[:, 2],
                      q=allrows[:, 3], beam=mask(0), a1=mask(1), a2=mask(2)))

pickle.dump({"nevents": len(ev), "cands": cands}, open(dst, "wb"))
print(f"  {os.path.basename(src)}: {len(ev)} tracked events -> {len(cands)} candidates "
      f"-> {os.path.basename(dst)}")
