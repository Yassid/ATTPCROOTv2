"""Langevin drift velocity for an arbitrary tilt azimuth.

The AT-TPC commissioning paper (eq. 4) writes the drift velocity as

    v_D = v/(1+wt^2) [ Ehat + wt (Ehat x Bhat) + wt^2 (Ehat.Bhat) Bhat ]

and then specialises it with B = B[sin(theta) yhat + cos(theta) zhat], i.e. it
*assumes* the tilt lies in the y-z plane. That is an assumption about the tilt's
azimuth, not just its magnitude, and it is what AtPSA::CalcLorentzVector encodes.

For the Dec 2014 data the beam (hence B) is measured at azimuth ~ -162 deg, almost
entirely along -x, so the specialised form applies the shear in the wrong direction.
Keeping phi general:

    Bhat = (sin(t)cos(p), sin(t)sin(p), cos(t)),  Ehat = zhat
    Ehat x Bhat = (-sin(t)sin(p), sin(t)cos(p), 0)
    Ehat . Bhat = cos(t)
"""
import math


def drift_vector(v_d, B, E, tilt_deg, azim_deg=90.0):
    """Langevin drift components (same units as v_d). azim_deg=90 -> paper's form."""
    if B == 0.0 or E is None:
        return 0.0, 0.0, v_d
    t, p = math.radians(tilt_deg), math.radians(azim_deg)
    wt = (B / E) * (v_d * 1e4)          # v_d in cm/us -> m/s
    f = v_d / (1.0 + wt * wt)
    st, ct, cp, sp = math.sin(t), math.cos(t), math.cos(p), math.sin(p)
    vx = f * (wt * (-st * sp) + wt * wt * ct * st * cp)
    vy = f * (wt * (st * cp) + wt * wt * ct * st * sp)
    vz = f * (1.0 + wt * wt * ct * ct)
    return vx, vy, vz
