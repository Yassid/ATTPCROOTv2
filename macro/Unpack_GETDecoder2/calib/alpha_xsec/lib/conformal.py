"""Conformal mapping: turn helices back into straight lines.

The problem this solves.  With B on, every track from the vertex is an arc, and a short arc
does not determine a circle -- over the ~90 mm a found fragment spans, every radius from
200 to 800 mm sits within 3.8 mm of sagitta.  Fitting circles directly therefore gives a
radius biased low by about a factor two and a tangent at the vertex wrong by 20-40 degrees,
measured against simulated truth.

The way out is a change of variables.  Shift the origin to the VERTEX, so every track is a
circle THROUGH the origin, and apply

    u = X / (X^2 + Y^2),     v = Y / (X^2 + Y^2)

A circle through the origin, X^2 + Y^2 = 2aX + 2bY, divides through to 1 = 2au + 2bv --
a STRAIGHT LINE.  So in conformal space the three tracks are three straight lines, which is
the problem the existing machinery was built for, and the line's coefficients give the
circle back exactly:  centre (a, b), radius sqrt(a^2 + b^2).

The radius then comes from a line fit over the whole track rather than from the sagitta of a
fragment, which is what breaks the degeneracy: a line's slope is well determined even when
the arc is short, because the map has spread the near-vertex hits far apart.
"""
import numpy as np

__all__ = ["to_conformal", "circle_from_line", "find_tracks_conformal"]


def to_conformal(P, vertex, r_min=12.0, scale=1.0e4):
    """(N,3) hits -> conformal (u, v) plus the mask of hits that were usable.

    Hits closer than r_min to the vertex are dropped: the map sends them to infinity, and
    they carry no information about the circle anyway.
    """
    X = P[:, 0] - vertex[0]
    Y = P[:, 1] - vertex[1]
    r2 = X * X + Y * Y
    keep = r2 > r_min * r_min
    u = scale * X[keep] / r2[keep]
    v = scale * Y[keep] / r2[keep]
    return np.column_stack([u, v]), keep


def circle_from_line(nrm, rhs, scale=1.0e4):
    """Invert the map.

    A circle through the origin with centre (a, b) satisfies X^2 + Y^2 = 2aX + 2bY, so
    dividing by r^2 gives 2a*u + 2b*v = 1 in UNSCALED conformal coordinates.  With
    u' = s*u the same line reads

        2a*u' + 2b*v' = s

    so if the fitted line is written with a UNIT normal as  nrm . (u', v') = rhs, then
    (2a, 2b) = nrm * s / rhs: the centre is nrm * s / (2 rhs) and the radius is its length.
    Dropping the factor s here is what made the first version come back 270 times too large.
    """
    if abs(rhs) < 1e-12:
        return None
    a = float(nrm[0]) * scale / (2.0 * rhs)
    b = float(nrm[1]) * scale / (2.0 * rhs)
    return a, b, float(np.hypot(a, b))


def _line_ransac(UV, n_iter=300, tol=None, min_hits=10, rng=None):
    """Greedy RANSAC for straight lines in the conformal plane."""
    rng = np.random.default_rng(0 if rng is None else rng)
    if tol is None:
        tol = 0.05 * np.median(np.hypot(UV[:, 0], UV[:, 1]))
    remaining = np.ones(len(UV), bool)
    out = []
    for _ in range(4):
        idx = np.flatnonzero(remaining)
        if len(idx) < min_hits:
            break
        Q = UV[idx]
        best = None
        for _ in range(n_iter):
            i, j = rng.choice(len(Q), size=2, replace=False)
            d = Q[j] - Q[i]
            n = np.hypot(*d)
            if n < 1e-12:
                continue
            nrm = np.array([-d[1] / n, d[0] / n])
            dist = np.abs((Q - Q[i]) @ nrm)
            inl = dist < tol
            if inl.sum() >= min_hits and (best is None or inl.sum() > best[0]):
                best = (int(inl.sum()), inl)
        if best is None:
            break
        m = np.zeros(len(UV), bool)
        m[idx[best[1]]] = True
        out.append(m)
        remaining[idx[best[1]]] = False
    return out


def find_tracks_conformal(P, vertex, q=None, r_min=12.0, min_hits=10, scale=1.0e4):
    """Find tracks as straight lines in conformal space.

    Returns a list of dicts: mask into P, and the circle (xc, yc, R) recovered by inverting
    the map.  The mask is in the ORIGINAL hit indexing.
    """
    UV, keep = to_conformal(P, vertex, r_min=r_min, scale=scale)
    if len(UV) < min_hits:
        return []
    kidx = np.flatnonzero(keep)
    out = []
    for m in _line_ransac(UV, min_hits=min_hits):
        if m.sum() < min_hits:
            continue
        # total least squares line through the conformal inliers, then invert
        Q = UV[m]
        c = Q.mean(axis=0)
        _, _, vt = np.linalg.svd(Q - c, full_matrices=False)
        d = vt[0]
        nrm = np.array([-d[1], d[0]])
        # line is  nrm . (x - c) = 0, i.e. nrm.x = nrm.c; put in the form 2(a u + b v) = 1
        rhs = float(nrm @ c)
        circ = circle_from_line(nrm, rhs, scale=scale)
        if circ is None:
            continue
        full = np.zeros(len(P), bool)
        full[kidx[m]] = True
        out.append(dict(mask=full, xc=circ[0] + vertex[0], yc=circ[1] + vertex[1], R=circ[2]))
    return out


def straightness(P, vertex, r_min=12.0, scale=1.0e4, n_dir=36):
    """How line-like is the conformal picture about this candidate vertex?

    In conformal coordinates about the TRUE vertex every track is a straight line; about a
    wrong vertex they are not.  So the vertex can be found with no seed, as the point that
    maximises straightness.

    The first version of this scored badly (median vertex error 145 mm) for two reasons
    worth recording, because both are easy to repeat:
      * it histogrammed with FIXED bins over a 1/r dynamic range, so almost every point
        landed in the central bins whatever the vertex, and the score barely moved;
      * it took the max over directions, rewarding ONE line, where the signal is that
        THREE lines exist at once.
    This version normalises the projection by its own spread before binning, and sums the
    best few directions instead of taking the single best.
    """
    UV, keep = to_conformal(P, vertex, r_min=r_min, scale=scale)
    if len(UV) < 12:
        return -1.0
    scores = []
    for th in np.linspace(0, np.pi, n_dir, endpoint=False):
        nrm = np.array([-np.sin(th), np.cos(th)])
        proj = UV @ nrm
        sd = np.std(proj)
        if sd <= 0:
            continue
        # scale-free: bin in units of the projection's own spread, so the score measures
        # concentration rather than the absolute size of the conformal cloud
        h, _ = np.histogram(proj / sd, bins=np.linspace(-3, 3, 31))
        scores.append(float(h.max()) / len(UV))
    if not scores:
        return -1.0
    scores.sort(reverse=True)
    return float(np.mean(scores[:3]))   # three tracks -> three good directions


def find_vertex(P, q=None, r_scan=None, n_coarse=13, n_fine=9, r_min=12.0):
    """Locate the vertex as the point maximising conformal straightness.

    Two passes: a coarse grid over the hit cloud, then a fine grid around the best cell.
    """
    P = np.asarray(P, float)
    if r_scan is None:
        r_scan = 0.5 * max(np.ptp(P[:, 0]), np.ptp(P[:, 1]))
    c = np.array([np.median(P[:, 0]), np.median(P[:, 1])])
    best = (-1.0, c)
    for npts, rad in ((n_coarse, r_scan), (n_fine, 2.0 * r_scan / max(n_coarse - 1, 1))):
        gx = np.linspace(best[1][0] - rad, best[1][0] + rad, npts)
        gy = np.linspace(best[1][1] - rad, best[1][1] + rad, npts)
        for x in gx:
            for y in gy:
                v = np.array([x, y])
                sc = straightness(P, v, r_min=r_min)
                if sc > best[0]:
                    best = (sc, v)
    return best[1], best[0]
