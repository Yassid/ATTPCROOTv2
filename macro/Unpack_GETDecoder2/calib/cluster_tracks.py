#!/usr/bin/env python3
"""Track clustering for the Dec2014 alpha events.

Plain HDBSCAN on (x,y,z) merges the two outgoing arms at the vertex, which is
exactly where we least want it to: the vertex is a real density maximum where
both tracks overlap. What separates them there is not position but *direction*.

So each hit gets a local direction from a PCA over its k nearest neighbours, and
we cluster in [position, orientation] space. A line direction is sign-ambiguous
(d and -d describe the same track), so the orientation is fed in as the
sign-invariant tensor d (x) d rather than as d itself.

The vertex is then recovered as the intersection of the fitted arms, not as a
cluster of its own -- which is what makes it well defined.
"""
import math
import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import HDBSCAN


def local_directions(P, k=10):
    """Unit local track direction at each hit, from PCA over k nearest neighbours."""
    k = max(3, min(k, len(P)))
    tree = cKDTree(P)
    _, idx = tree.query(P, k=k)
    dirs = np.empty_like(P, dtype=float)
    lin = np.empty(len(P))
    for i in range(len(P)):
        nb = P[idx[i]]
        nb = nb - nb.mean(axis=0)
        # singular values give how line-like the neighbourhood is
        _, s, vt = np.linalg.svd(nb, full_matrices=False)
        dirs[i] = vt[0]
        lin[i] = s[0] / max(s[1], 1e-9)
    return dirs, lin


def orientation_tensor(dirs):
    """Sign-invariant encoding of a direction: the 6 independent parts of d (x) d."""
    x, y, z = dirs[:, 0], dirs[:, 1], dirs[:, 2]
    return np.column_stack([x * x, y * y, z * z,
                            math.sqrt(2) * x * y,
                            math.sqrt(2) * x * z,
                            math.sqrt(2) * y * z])


def cluster(P, q=None, k=10, dir_weight=120.0, min_cluster_size=6, pos_scale=40.0):
    """Cluster hits into tracks. Returns labels (-1 = noise) and the local dirs.

    dir_weight sets how much the orientation counts relative to position: it is
    the distance in mm that a fully orthogonal direction difference is worth.
    """
    dirs, lin = local_directions(P, k)
    feat = np.hstack([P / pos_scale, orientation_tensor(dirs) * (dir_weight / pos_scale)])
    lab = HDBSCAN(min_cluster_size=min_cluster_size,
                  min_samples=3,
                  cluster_selection_method="eom").fit_predict(feat)
    return lab, dirs, lin


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


def closest_approach(c1, d1, c2, d2):
    """Midpoint of the common perpendicular of two 3D lines, and their distance."""
    w0 = c1 - c2
    a, b, c = d1 @ d1, d1 @ d2, d2 @ d2
    d_, e = d1 @ w0, d2 @ w0
    den = a * c - b * b
    if abs(den) < 1e-12:
        return 0.5 * (c1 + c2), np.linalg.norm(np.cross(d1, w0))
    s = (b * e - c * d_) / den
    t = (a * e - b * d_) / den
    p1, p2 = c1 + s * d1, c2 + t * d2
    return 0.5 * (p1 + p2), np.linalg.norm(p1 - p2)


def vertex_from_tracks(tracks):
    """Least-squares point closest to all track lines. tracks = [(c,d),...]."""
    A = np.zeros((3, 3))
    b = np.zeros(3)
    for c, d in tracks:
        M = np.eye(3) - np.outer(d, d)      # projector orthogonal to the line
        A += M
        b += M @ c
    return np.linalg.solve(A, b)
