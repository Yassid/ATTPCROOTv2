#!/usr/bin/env python3
"""alpha+alpha elastic excitation function from the 68 magnet-OFF runs.

Thick-target inverse kinematics: the beam enters at 7.80 MeV and slows continuously, so the
vertex position IS an energy measurement and one run scans E_cm from 3.9 MeV down to 0.

B = 0 makes this the cleanest possible sample: omega*tau = 0, so there is no E x B lateral
drift and the tracks are straight, meaning a full-arm line fit is the tangent at the vertex.
Both the chord-vs-tangent bias and the Lorentz correction -- the two model-dependent steps --
are simply absent.

The tabulated energy loss is for 150 torr; these runs are at 299.5 torr. Stopping power is
proportional to electron density, hence to pressure at fixed temperature, so
E(d, P) = E_150(d * P / 150) exactly, with no refitting.

Selection follows the campaign fit: opening angle within mu +- 2 sigma of 90 deg, where the
window is 99% pure, rather than the old 55-125 deg window which is only 90% pure.
"""
import glob, os, pickle, sys, math, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from calibrate_run import geometry

VD, B, E, TILT0 = 2.2373, 0.0, 30000.0, 7.0
P_RUN, P_TABLE = 299.5, 150.0
Z_WINDOW = 1000.0
MU, SIG = 89.987, 1.56          # measured on this very sample

def eloss():
    d, e = [], []
    for line in open("/home/yassid/dec2014_calib/eloss_alpha_heco2.txt"):
        if line.startswith("#"):
            continue
        a, b = line.split(); d.append(float(a)); e.append(float(b))
    return np.array(d), np.array(e)

DIST, ENER = eloss()

def run_one(path):
    run = int(os.path.basename(path).replace("cands_", "").replace(".pkl", ""))
    d = pickle.load(open(path, "rb"))
    ecm, zv, ops_out = [], [], []
    for rec in d["cands"]:
        try:
            o, bt, vtx, _ = geometry(rec, VD, B, E, TILT0)
        except Exception:
            continue
        if abs(o - MU) > 2 * SIG:            # 99% pure elastic window
            continue
        z = vtx[2]
        if not (0 < z < Z_WINDOW):
            continue
        # gas traversed before the reaction, converted to the 150 torr table's scale
        path_len = (Z_WINDOW - z) * P_RUN / P_TABLE
        eb = float(np.interp(path_len, DIST, ENER))
        if eb <= 0.05:
            continue
        ecm.append(eb / 2.0); zv.append(z); ops_out.append(o)
    return run, np.array(ecm), np.array(zv), np.array(ops_out)

if __name__ == "__main__":
    files = sorted(glob.glob("/home/yassid/dec2014_calib/campaign/P300_B0/cands_*.pkl"))
    with Pool(12) as p:
        res = p.map(run_one, files)
    ecm = np.concatenate([r[1] for r in res]); zv = np.concatenate([r[2] for r in res])
    print(f"elastic events with usable vertex: {len(ecm)}")
    print(f"  vertex z  : {zv.min():.0f} - {zv.max():.0f} mm, median {np.median(zv):.0f}")
    print(f"  E_cm      : {ecm.min():.2f} - {ecm.max():.2f} MeV, median {np.median(ecm):.2f}")
    np.save("/home/yassid/dec2014_calib/campaign/ecm_b0.npy", ecm)
    np.save("/home/yassid/dec2014_calib/campaign/zvtx_b0.npy", zv)
