#!/usr/bin/env python
"""Greedy KF-style spiral track follower.

Seed = PRA's largest track arc. Then extend from BOTH ends by local-curvature
following (predict -> search window -> add -> refit local circle), over the
unassigned hits. The local circle refit each step makes the radius shrink with
the spiral (energy loss handled implicitly), and bidirectional extension recovers
the dropped START of the track (the ev57 failure). No global model assumed.

Reusable API: follow_event(hits_df_for_one_event, params) -> dict with the
extended track hit-index set and diagnostics.
"""
import numpy as np

# ---------- geometry helpers ----------
def fit_circle(x, y):
    """Algebraic (Kasa) circle fit -> (cx, cy, R). Returns None if degenerate."""
    A = np.c_[2*x, 2*y, np.ones(len(x))]
    b = x*x + y*y
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    cx, cy, c = sol
    r2 = c + cx*cx + cy*cy
    if r2 <= 0:
        return None
    return cx, cy, np.sqrt(r2)

def order_along_arc(P, cx, cy):
    """Order points sequentially along the arc around (cx,cy). Sort by angle, then
    start the sequence right after the largest angular gap so a partial arc (even
    one crossing +-pi) comes out contiguous."""
    ang = np.arctan2(P[:,1]-cy, P[:,0]-cx)
    idx = np.argsort(ang)
    sa = ang[idx]
    gaps = np.diff(np.r_[sa, sa[0] + 2*np.pi])   # gap to next, wrapping
    cut = int(np.argmax(gaps))                    # largest gap -> arc break
    order = np.r_[idx[cut+1:], idx[:cut+1]]       # roll so sequence starts after the gap
    return order, ang[order]

# ---------- the follower ----------
DEFAULTS = dict(
    K=14,          # hits in the local fit window
    step=4.0,      # arc-length advance per prediction [mm]
    win=11.0,      # search radius around the prediction [mm]
    maxmiss=6,     # stop after this many consecutive empty predictions
    maxsteps=4000, # hard cap on prediction steps
    min_seed=12,   # need at least this many seed hits to follow
    rmin=8.0,      # ignore predictions whose local radius collapses below this
)

def _tangent(win_pts):
    """Local 3D tangent via PCA of the window, oriented start->end."""
    c = win_pts.mean(axis=0)
    Q = win_pts - c
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    d = Vt[0]
    if np.dot(d, win_pts[-1] - win_pts[0]) < 0:
        d = -d
    return d / (np.linalg.norm(d) + 1e-9)

def _extend(seed_xyz, pool_xyz, pool_idx, used, direction, par):
    """Walk from one end of the seed using a local PCA tangent (robust on thick,
    curving point clouds). direction=+1 leading end, -1 trailing end."""
    K = par['K']
    grabbed = []
    win_pts = seed_xyz[-K:].copy() if direction > 0 else seed_xyz[:K][::-1].copy()
    miss = 0
    for _ in range(par['maxsteps']):
        d = _tangent(win_pts)
        P = win_pts[-1]
        pred = P + par['step'] * d
        free = ~used
        if not free.any():
            break
        cand = pool_xyz[free]
        # gate: forward of the current point (don't walk backward into the seed)
        fwd = (cand - P) @ d
        dist = np.linalg.norm(cand - pred, axis=1)
        ok = fwd > -par['step']                       # allow slight lateral
        if not ok.any():
            miss += 1
            if miss >= par['maxmiss']: break
            win_pts = np.vstack([win_pts, pred])[-K:]
            continue
        dist_ok = np.where(ok, dist, np.inf)
        j = int(np.argmin(dist_ok))
        if dist_ok[j] > par['win']:
            miss += 1
            if miss >= par['maxmiss']:
                break
            win_pts = np.vstack([win_pts, pred])[-K:]   # coast over a gap
            continue
        miss = 0
        gi = np.where(free)[0][j]
        used[gi] = True
        grabbed.append(pool_idx[gi])
        win_pts = np.vstack([win_pts, pool_xyz[gi]])[-K:]
    return grabbed

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

GROW = dict(
    dmerge=10.0,    # connectivity radius for region growing [mm]
    dbridge=45.0,   # max gap a directional bridge will jump [mm]
    Kend=14,        # hits used to estimate the local tangent at an endpoint
    min_dense=3,    # a hit is absorbable only if it has >=this many neighbors (rejects noise)
    min_seed=12,
    max_iter=40,
)

def _endpoints(track_xyz):
    """Two extreme points of the track along its 1st PCA axis -> (idx_lo, idx_hi)."""
    c = track_xyz.mean(axis=0)
    Q = track_xyz - c
    _, _, Vt = np.linalg.svd(Q, full_matrices=False)
    proj = Q @ Vt[0]
    return int(np.argmin(proj)), int(np.argmax(proj))

def _bridge_from(end_xyz, track_xyz, pool_xyz, free, dense_ok, par):
    """From a track endpoint, jump to the nearest FORWARD *dense* unassigned hit
    within dbridge (outward continuation of the local tangent). Bridging only to a
    dense hit makes the jump land in the real track body, not on stray noise."""
    dt = np.linalg.norm(track_xyz - end_xyz, axis=1)
    near = track_xyz[np.argsort(dt)[:par['Kend']]]
    tang = _tangent(near)
    if np.dot(tang, end_xyz - near.mean(axis=0)) < 0:   # orient outward
        tang = -tang
    cand_mask = free & dense_ok
    if not cand_mask.any():
        return None
    cand = pool_xyz[cand_mask]
    fwd = (cand - end_xyz) @ tang
    dist = np.linalg.norm(cand - end_xyz, axis=1)
    ok = (fwd > 0) & (dist < par['dbridge'])
    if not ok.any():
        return None
    di = np.where(ok, dist, np.inf)
    j = int(np.argmin(di))
    return int(np.where(cand_mask)[0][j])

def follow_event(ev_hits, par=None):
    """Grow + bridge + grow recovery of a dropped track around the PRA seed.
    ev_hits: dict with arrays x,y,z,praTrk for ONE event."""
    p = dict(GROW)
    if par: p.update(par)
    P = np.c_[ev_hits['x'], ev_hits['y'], ev_hits['z']]
    pra = ev_hits['praTrk']
    seed_ids = [t for t in set(pra) if t >= 0]
    if not seed_ids:
        return dict(ok=False, reason='no PRA seed', seed=0, extended=0, grabbed=0)
    best = max(seed_ids, key=lambda t: (pra == t).sum())
    intrk = (pra == best).copy()
    n_seed = int(intrk.sum())
    if n_seed < p['min_seed']:
        return dict(ok=False, reason='seed too small', seed=n_seed, extended=n_seed, grabbed=0)
    tree = cKDTree(P)
    # density per hit: reject isolated noise from region growing
    dense_ok = tree.query_ball_point(P, r=p['dmerge'], return_length=True) >= p['min_dense']

    def grow():
        frontier = np.where(intrk)[0].tolist()
        while frontier:
            nb = tree.query_ball_point(P[frontier], r=p['dmerge'])
            nxt = []
            for lst in nb:
                for j in lst:
                    if not intrk[j] and dense_ok[j]:
                        intrk[j] = True; nxt.append(j)
            frontier = nxt

    for _ in range(p['max_iter']):
        grow()
        track_xyz = P[intrk]
        free = ~intrk
        if not free.any():
            break
        lo, hi = _endpoints(track_xyz)
        added = False
        for e in (track_xyz[lo], track_xyz[hi]):
            gi = _bridge_from(e, track_xyz, P, free, dense_ok, p)
            if gi is not None and not intrk[gi]:
                intrk[gi] = True; added = True
        if not added:
            break

    grabbed_idx = [i for i in np.where(intrk)[0] if not (pra[i] == best)]
    return dict(ok=True, seed=n_seed, extended=int(intrk.sum()),
                grabbed=len(grabbed_idx), grabbed_idx=grabbed_idx,
                seed_idx=np.where(pra == best)[0].tolist())
