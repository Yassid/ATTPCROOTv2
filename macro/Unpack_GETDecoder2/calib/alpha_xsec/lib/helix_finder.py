#!/usr/bin/env python3
"""Pattern recognition for the MAGNET-ON settings: find helices, not lines.

Why neither existing finder works here.

  * lib/clustering.py clusters in [position, ORIENTATION] space.  That is what makes it
    good at separating two straight arms meeting at a vertex, and it is exactly wrong for a
    curved track: the local direction rotates continuously along a helix, so hits at the two
    ends of one track look like different objects and the track is shattered into short
    arcs.  Measured on real 0.57 T data it returned fitted radii of ~23 mm where 400-600 mm
    is expected.

  * ATTPCROOT's AtSampleConsensus with kCircle2D fits the XY projection and IGNORES z in its
    distance metric, so it happily links hits from different tracks that happen to overlap
    in projection.  On the same data it gave ~47 mm radii and "arc lengths" of 873 mm, i.e.
    nearly three turns, which no alpha here can do.

So the model has to be a genuine 3D helix: a circle in the pad plane AND z advancing
linearly with arc length, with both entering the inlier test.

The sampling is the standard RANSAC one -- three hits define a circle -- plus the dip from
a fit of z against arc length.  A candidate is kept only if BOTH residuals are small, which
is what stops it linking unrelated tracks.
"""
import numpy as np

__all__ = ["find_helices"]


def _circle_from_three(p1, p2, p3):
    """Centre and radius of the circle through three 2D points, or None if collinear."""
    ax, ay = p1; bx, by = p2; cx, cy = p3
    d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-9:
        return None
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    return ux, uy, float(np.hypot(ax - ux, ay - uy))


def _dip_residual(P, xc, yc, R):
    """Fit z against arc length about (xc, yc); return per-hit residual and the slope.

    This is the half that kCircle2D leaves out, and it is what prevents hits from a
    different track at a different z being counted as inliers.
    """
    phi = np.arctan2(P[:, 1] - yc, P[:, 0] - xc)
    # order by angle before unwrapping, otherwise the unwrap is meaningless for a point set
    o = np.argsort(phi)
    s = np.empty(len(P)); s[o] = R * np.unwrap(phi[o])
    A = np.column_stack([s, np.ones(len(s))])
    sol, *_ = np.linalg.lstsq(A, P[:, 2], rcond=None)
    return P[:, 2] - A @ sol, sol[0], s


def find_helices(P, q=None, n_iter=400, tol_xy=6.0, tol_z=12.0,
                 min_hits=12, max_tracks=4, r_min=50.0, r_max=5000.0, rng=None):
    """Greedy RANSAC over 3D helices.

    P        (N, 3) hit positions, z along the field
    tol_xy   inlier tolerance on the circle, mm
    tol_z    inlier tolerance on z against the arc-length fit, mm
    returns  list of boolean masks into P, largest track first
    """
    P = np.asarray(P, float)
    rng = np.random.default_rng(0 if rng is None else rng)
    q = np.ones(len(P)) if q is None else np.asarray(q, float)
    remaining = np.ones(len(P), bool)
    out = []

    for _ in range(max_tracks):
        idx = np.flatnonzero(remaining)
        if len(idx) < min_hits:
            break
        Q = P[idx]
        best = None
        for _ in range(n_iter):
            pick = rng.choice(len(Q), size=3, replace=False)
            c = _circle_from_three(Q[pick[0], :2], Q[pick[1], :2], Q[pick[2], :2])
            if c is None:
                continue
            xc, yc, R = c
            if not (r_min < R < r_max):
                continue
            dxy = np.abs(np.hypot(Q[:, 0] - xc, Q[:, 1] - yc) - R)
            near = dxy < tol_xy
            if near.sum() < min_hits:
                continue
            # the dip is fitted on the circle's inliers only, then applied to all of them
            dz, _, _ = _dip_residual(Q[near], xc, yc, R)
            keep = np.zeros(len(Q), bool)
            sub = np.flatnonzero(near)[np.abs(dz) < tol_z]
            keep[sub] = True
            n = keep.sum()
            if n >= min_hits and (best is None or n > best[0]):
                best = (n, keep, xc, yc, R)
        if best is None:
            break
        n, keep, xc, yc, R = best
        mask = np.zeros(len(P), bool)
        mask[idx[keep]] = True
        out.append(mask)
        remaining[idx[keep]] = False
    return out
