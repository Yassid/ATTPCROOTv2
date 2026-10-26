"""Turning a stored candidate event into physics quantities.

A candidate pickle holds, per event, every hit of the event (tb, x, y, q) plus three
boolean masks naming the beam segment and the two arms.  Nothing is pre-computed in 3D,
because z depends on the drift velocity, which is still being pinned down -- so every
geometric quantity is rebuilt here from (x, y, tb) for whatever V_DRIFT is in force.
"""
import numpy as np


def fit_line(P, w=None):
    """Charge-weighted total-least-squares 3D line -> (centroid, unit direction)."""
    w = np.ones(len(P)) if w is None else np.asarray(w, dtype=float)
    w = w / w.sum()
    c = (P * w[:, None]).sum(axis=0)
    D = P - c
    C = (D * w[:, None]).T @ D
    val, vec = np.linalg.eigh(C)
    d = vec[:, np.argmax(val)]
    return c, d / np.linalg.norm(d)


def vertex_from_tracks(tracks):
    """Least-squares point closest to all the given lines.  tracks = [(centroid, dir), ...]"""
    A = np.zeros((3, 3)); b = np.zeros(3)
    for c, d in tracks:
        M = np.eye(3) - np.outer(d, d)
        A += M; b += M @ c
    return np.linalg.solve(A, b)


def positions(rec, v_drift, t0_buckets, tb_ns=160.0):
    """3D hit positions in mm.  x, y are pad coordinates; z comes from the drift time."""
    z = v_drift * (rec["tb"] - t0_buckets) * tb_ns * 1e-3 * 10.0
    return np.column_stack([rec["x"], rec["y"], z])


def reconstruct(rec, v_drift, t0_buckets, tb_ns=160.0, tilt_deg=6.47):
    """
    Rebuild one event.  Returns None if it cannot be reconstructed, else a dict with:

      vertex     3D reaction point, from the two arms
      entrance   where the beam entered, see below
      beam_dir   unit vector along the beam's direction of travel
      theta      (theta1, theta2) each arm's LAB angle to the beam, degrees
      theta_sum  theta1 + theta2 -- must be 90 deg for equal-mass elastic
      length     (L1, L2) each arm's length in mm
      path       distance the beam travelled from entrance to vertex, mm
      opening    angle between the two arms, degrees

    The beam segment is OPTIONAL.  If rec["beam"] names fewer than two hits the entrance
    is rebuilt from the nominal tilt instead -- see below.  Only the two arms are required,
    because they are what fixes the vertex.

    ENTRANCE.  The trigger is not synchronised to the beam, so each event carries its own
    unknown t0 and absolute z is meaningless event by event.  What is reliable: the tilt
    pivots at the entrance window, so the beam enters ON the detector axis and only
    diverges as it goes deeper.  The entrance is therefore the point where the fitted beam
    line comes closest to x = y = 0 -- per event, with no timing assumption.
    """
    P = positions(rec, v_drift, t0_buckets, tb_ns)
    seg = {}
    for k in ("beam", "a1", "a2"):
        m = rec[k]
        if m.sum() < 2:
            if k == "beam":
                seg["beam"] = None          # no beam track: fall back to the nominal tilt
                continue
            return None
        seg[k] = fit_line(P[m], rec["q"][m])

    vertex = vertex_from_tracks([seg["a1"], seg["a2"]])

    if seg["beam"] is not None:
        cb, db = seg["beam"]
        den = db[0] ** 2 + db[1] ** 2
        if den < 1e-12:                   # beam exactly along z: axis crossing undefined
            return None
        t = -(cb[0] * db[0] + cb[1] * db[1]) / den
        entrance = cb + t * db
    else:
        # THE BEAM WAS NOT FOUND.  The tilt pivots at the window, so the beam enters on the
        # detector axis and leaves it only with depth: the transverse distance of the vertex
        # from the axis already fixes how far the beam travelled, with no timing at all.
        #   r = |(Vx, Vy)|,  path = r / sin(tilt),  dz = r / tan(tilt)
        # Two entrances satisfy that, one on each side in z.  The right one is the one the
        # arms agree with: for equal-mass elastic scattering both arms go forward of 90 deg,
        # so both must have a POSITIVE projection on the beam direction.  That picks the
        # side without needing to know which way z runs in this frame.
        r = float(np.hypot(vertex[0], vertex[1]))
        if r < 1e-6:
            return None
        dz = r / np.tan(np.radians(tilt_deg))
        best, entrance = None, None
        for sgn in (-1.0, +1.0):
            e = np.array([0.0, 0.0, vertex[2] + sgn * dz])
            w = vertex - e
            nw = np.linalg.norm(w)
            if nw < 1e-9:
                continue
            u = w / nw
            proj = []
            for k in ("a1", "a2"):
                pts = P[rec[k]]
                c, d = seg[k]
                tv = float(vertex @ d)
                far = pts[np.argmax(np.abs(pts @ d - tv))]
                v = far - vertex
                nv = np.linalg.norm(v)
                if nv < 1e-9:
                    proj = None; break
                proj.append(float((v / nv) @ u))
            if proj and (best is None or min(proj) > best):
                best, entrance = min(proj), e
        if entrance is None or best is None or best <= 0:
            return None

    bd = vertex - entrance
    n = np.linalg.norm(bd)
    if n < 1e-9:
        return None
    bd = bd / n                           # direction of travel, entrance -> vertex

    theta, length, udir = [], [], []
    for k in ("a1", "a2"):
        m = rec[k]
        pts = P[m]
        c, d = seg[k]
        proj = pts @ d
        tv = float(vertex @ d)
        far = pts[np.argmax(np.abs(proj - tv))]     # the end furthest from the vertex
        w = far - vertex
        nn = np.linalg.norm(w)
        if nn < 1e-9:
            return None
        u = w / nn
        udir.append(u)
        theta.append(np.degrees(np.arccos(np.clip(u @ bd, -1, 1))))
        length.append(float(np.ptp(proj)))

    opening = np.degrees(np.arccos(np.clip(udir[0] @ udir[1], -1, 1)))
    return dict(vertex=vertex, entrance=entrance, beam_dir=bd,
                theta=tuple(theta), theta_sum=theta[0] + theta[1],
                length=tuple(length), path=float(n), opening=float(opening))
