"""Helix geometry, for the magnet-on settings.

With B along the drift axis a track is a HELIX: a circle in the pad plane (x, y) advancing
linearly in z.  Everything the B = 0 analysis does with straight lines has to be redone,
because the sagitta is 5 - 59 mm against a 2.8 mm pad pitch -- a line fit is wrong by ten
to twenty times the granularity, not by a little.

What the curvature BUYS, and the reason this is worth doing: the radius gives the momentum
directly,

    p_T [MeV/c] = 0.29979 * B[T] * q * R[mm]           (q = 2 for an alpha)

so the magnet-on data carry an energy measurement that needs neither the range table nor
the reconstructed path length.  That is a third, independent handle on the energy scale
that the B = 0 data simply do not have.
"""
import numpy as np

M_ALPHA = 3727.379     # MeV
Q_ALPHA = 2


def fit_circle(x, y, w=None):
    """Charge-weighted algebraic circle fit (Kasa).  -> (xc, yc, R, rms residual)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    w = w / w.sum()
    # solve  2 x xc + 2 y yc + c = x^2 + y^2   in a weighted least-squares sense
    A = np.column_stack([2 * x, 2 * y, np.ones(len(x))])
    b = x * x + y * y
    W = w[:, None]
    sol, *_ = np.linalg.lstsq(A * W, b * w, rcond=None)
    xc, yc, c = sol
    R2 = c + xc * xc + yc * yc
    if R2 <= 0:
        return xc, yc, np.inf, np.inf
    R = np.sqrt(R2)
    r = np.hypot(x - xc, y - yc) - R
    return xc, yc, R, float(np.sqrt(np.sum(w * r * r)))


def fit_line_xy(x, y, w=None):
    """Weighted total-least-squares line in the pad plane -> rms perpendicular residual.

    Only used to demonstrate that a line does NOT describe these tracks.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    w = w / w.sum()
    cx, cy = np.sum(w * x), np.sum(w * y)
    dx, dy = x - cx, y - cy
    C = np.array([[np.sum(w * dx * dx), np.sum(w * dx * dy)],
                  [np.sum(w * dx * dy), np.sum(w * dy * dy)]])
    val, vec = np.linalg.eigh(C)
    d = vec[:, np.argmax(val)]
    n = np.array([-d[1], d[0]])
    r = dx * n[0] + dy * n[1]
    return float(np.sqrt(np.sum(w * r * r)))


def fit_helix(P, w=None):
    """Fit a helix to 3D points.  P is (N, 3) with z along the field.

    Returns a dict with the circle, the dip angle, the arc length and the momentum.
    """
    P = np.asarray(P, float)
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    xc, yc, R, rms = fit_circle(x, y, w)
    if not np.isfinite(R):
        return None
    # unwrapped turning angle about the circle centre, then z against arc length
    phi = np.unwrap(np.arctan2(y - yc, x - xc))
    s = R * phi                                    # signed arc length in the pad plane
    A = np.column_stack([s, np.ones(len(s))])
    (dzds, z0), *_ = np.linalg.lstsq(A, z, rcond=None)
    lam = np.arctan(dzds)                          # dip angle out of the pad plane
    arc = float(np.ptp(s) / max(np.cos(lam), 1e-6))  # true 3D path along the helix
    return dict(xc=xc, yc=yc, R=R, rms=rms, phi=phi, s=s, lam=lam,
                dzds=dzds, z0=z0, arc=arc)


def tangent_at(fit, x, y, sense=+1):
    """Unit tangent to the helix at a point on it -- the DIRECTION at the vertex.

    For a curved track the direction is a local quantity: a global fit direction is
    meaningless, and this is what replaces it when computing theta.
    """
    ux, uy = x - fit["xc"], y - fit["yc"]
    n = np.hypot(ux, uy)
    if n < 1e-9:
        return None
    # tangent to the circle, times the dip
    tx, ty = -uy / n, ux / n
    tz = fit["dzds"]
    v = np.array([tx, ty, tz], float) * sense
    return v / np.linalg.norm(v)


def momentum_MeV(R_mm, lam, B_T, q=Q_ALPHA):
    """Momentum from the curvature.  R in mm, B in tesla."""
    # p[GeV/c] = 0.29979 * B[T] * R[m] * q, so with R in mm the 1/1000 and the GeV->MeV
    # 1000 cancel exactly and the constant is used as-is.  Getting this wrong by the
    # obvious factor of ten puts an 11 MeV alpha at 0.1 MeV, which is how it was caught.
    pT = 0.29979 * B_T * q * R_mm                # MeV/c, R in mm
    return pT / max(np.cos(lam), 1e-6)


def energy_from_radius(R_mm, lam, B_T, q=Q_ALPHA, m=M_ALPHA):
    """Kinetic energy implied by the curvature -- independent of any range table."""
    p = momentum_MeV(R_mm, lam, B_T, q)
    return float(np.sqrt(p * p + m * m) - m)
