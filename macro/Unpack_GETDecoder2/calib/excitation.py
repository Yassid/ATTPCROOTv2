#!/usr/bin/env python3
"""alpha+alpha elastic excitation function from the Dec 2014 thick-target data.

The beam enters the gas at 7.8 MeV and slows continuously, so the vertex position
*is* an energy measurement: a single run scans E_cm from 3.9 MeV down to 0. That is
the thick-target inverse-kinematics trick, and it makes the excitation function fall
out of the vertex distribution we already reconstruct.

  E_beam(vertex) = energy after traversing (z_window - z_vertex) of He:CO2
  E_cm           = E_beam / 2          (equal masses, non-relativistic)
  E_x(8Be)       = E_cm + Q            with Q = 0 for alpha+alpha elastic

The energy-loss table comes from CATIMA via the port in 926468f8 (dumped to text by
eloss_table.cxx so this script has no ROOT dependency).

Why this is a real test of the calibration and not just an analysis: the resonance
energy depends on the vertex position, which depends on v_D, the tilt and the Lorentz
correction. A 12% error in v_D would move the 8Be 2+ by several hundred keV. No
geometric check can catch an error common to the whole chain; this can.
"""
import math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")

Z_WINDOW = 1000.0          # entrance window, mm (detector frame, pad plane at z=0)
E_BEAM0 = 7.80             # MeV, 1.95 MeV/u * 4
M_ALPHA = 3727.379         # MeV/c^2


def load_eloss(path="/home/yassid/dec2014_calib/eloss_alpha_heco2.txt"):
    """Energy remaining after a given distance, from the CATIMA dump."""
    d, e = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        a, b = line.split()
        d.append(float(a)); e.append(float(b))
    return np.array(d), np.array(e)


def main():
    dist, ener = load_eloss(sys.argv[1] if len(sys.argv) > 1 else
                            "/home/yassid/dec2014_calib/eloss_alpha_heco2.txt")
    res = pickle.load(open("/home/yassid/dec2014_calib/prod_results.pkl", "rb"))

    ecm, runs = [], []
    for r in res:
        v = r["vtx"]
        ops = r["ops"]
        n = min(len(v), len(ops))
        for i in range(n):
            if not (55 < ops[i] < 125):          # alpha+alpha elastic selection
                continue
            z = v[i][2]
            if not (0 < z < Z_WINDOW):
                continue
            path_len = Z_WINDOW - z              # mm of gas traversed before reacting
            eb = float(np.interp(path_len, dist, ener))
            if eb <= 0:
                continue
            ecm.append(eb / 2.0)
            runs.append(r["run"])
    ecm = np.array(ecm)
    print(f"{len(ecm)} elastic events with a usable vertex")
    print(f"E_cm range {ecm.min():.2f} - {ecm.max():.2f} MeV, median {np.median(ecm):.2f}")
    np.save("/home/yassid/dec2014_calib/ecm.npy", ecm)


if __name__ == "__main__":
    main()
