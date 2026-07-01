#!/usr/bin/env python
"""Manual point-cloud labeler for real a1975 16C(d,p) events (GNN gold set).

Lasso hits on either 2D projection and assign a class. Semantic labels match the
overlay convention:  0 = proton , 1 = 17C recoil , 2 = beam/noise (default).
Optional HDBSCAN "guide" coloring reveals track structure before you lasso.

Run (GUI, needs WSLg display + PySide6 in gnn_env):
    ~/gnn_env/bin/python label_tool.py
Smoke test (no GUI, renders event 0 to a PNG and exits):
    ~/gnn_env/bin/python label_tool.py --smoke

Workflow per event: lasso the proton track -> press P ; (optional) lasso 17C -> R ;
everything else stays beam/noise. Press N for next (autosaves). Output: data/labels.parquet
(columns event,x,y,z,q,label,reviewed). Only reviewed=True events are meant for training.

Keys:
  P proton | R 17C | B beam/noise   (class to paint the next lasso / last selection)
  G toggle HDBSCAN guide colors      C reset event to all beam/noise
  D toggle 'reviewed' flag           N / -> next event     M / <- prev event
  S save now                         Q quit (autosaves)
Mouse: drag to lasso on either panel. Point size scales with charge.
"""
import argparse, sys, numpy as np, pandas as pd
from matplotlib.path import Path

CLASS = {0: ("proton", "#2f7bff"), 1: ("17C", "#e0301e"), 2: ("beam/noise", "#b7b7b7")}
IN_PARQUET  = "data/label_input.parquet"
OUT_PARQUET = "data/labels.parquet"


def load():
    df = pd.read_parquet(IN_PARQUET)
    events = []
    for ev, g in df.groupby("event", sort=True):
        g = g.reset_index(drop=True)
        xyz = g[["x", "y", "z"]].to_numpy(float)
        q = g["q"].to_numpy(float)
        events.append({"ev": int(ev), "xyz": xyz, "q": q,
                       "label": np.full(len(g), 2, np.int8), "reviewed": False})
    # resume from prior labels if present and consistent
    try:
        prev = pd.read_parquet(OUT_PARQUET)
        for e in events:
            p = prev[prev.event == e["ev"]]
            if len(p) == len(e["label"]):
                e["label"] = p["label"].to_numpy(np.int8)
                e["reviewed"] = bool(p["reviewed"].iloc[0])
        print(f"resumed: {int(prev.groupby('event').reviewed.first().sum())} events already reviewed")
    except FileNotFoundError:
        pass
    return events


def save(events):
    frames = []
    for e in events:
        frames.append(pd.DataFrame({
            "event": e["ev"], "x": e["xyz"][:, 0], "y": e["xyz"][:, 1], "z": e["xyz"][:, 2],
            "q": e["q"], "label": e["label"], "reviewed": e["reviewed"]}))
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(OUT_PARQUET, index=False)
    nrev = int(out.groupby("event").reviewed.first().sum())
    return nrev


def sizes(q):
    s = np.sqrt(np.clip(q, 1, None))
    return 6 + 34 * (s - s.min()) / (np.ptp(s) + 1e-9)


def guide_colors(xyz):
    from sklearn.cluster import HDBSCAN
    import matplotlib.cm as cm
    X = xyz.copy()
    X[:, :2] *= 1.0            # keep transverse & drift comparable-ish
    lab = HDBSCAN(min_cluster_size=15, min_samples=6).fit_predict(X)
    pal = cm.tab20(np.linspace(0, 1, 20))
    cols = np.array([pal[l % 20] if l >= 0 else (0.85, 0.85, 0.85, 1.0) for l in lab])
    return cols


def render(e, ax_xy, ax_zy, guide=False):
    xyz, q = e["xyz"], e["q"]
    if guide:
        cols = guide_colors(xyz)
    else:
        cols = np.array([CLASS[int(l)][1] for l in e["label"]])
    sz = sizes(q)
    ax_xy.clear(); ax_zy.clear()
    ax_xy.scatter(xyz[:, 0], xyz[:, 1], s=sz, c=cols, edgecolors="none")
    ax_zy.scatter(xyz[:, 2], xyz[:, 1], s=sz, c=cols, edgecolors="none")
    ax_xy.set_title("pad plane  X vs Y"); ax_xy.set_xlabel("x [mm]"); ax_xy.set_ylabel("y [mm]")
    ax_xy.set_aspect("equal", "box")
    ax_zy.set_title("drift  Z vs Y"); ax_zy.set_xlabel("z [mm]"); ax_zy.set_ylabel("y [mm]")


def counts(e):
    u, c = np.unique(e["label"], return_counts=True)
    d = dict(zip(u.tolist(), c.tolist()))
    return f"p={d.get(0,0)} 17C={d.get(1,0)} bg={d.get(2,0)}"


def smoke():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    events = load()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.5))
    # fake a proton lasso on event 0 to exercise the label path
    e = events[0]
    off = np.hypot(e["xyz"][:, 0], e["xyz"][:, 1]) >= 60
    e["label"][off] = 0
    render(e, a1, a2)
    fig.suptitle(f"SMOKE event {e['ev']}  n={len(e['q'])}  {counts(e)}")
    plt.tight_layout(); out = "data/smoke_event0.png"; plt.savefig(out, dpi=110)
    print(f"smoke OK: {len(events)} events loaded, rendered {out}")


def gui():
    import matplotlib
    matplotlib.use("QtAgg")
    import matplotlib.pyplot as plt
    from matplotlib.widgets import LassoSelector
    events = load()
    st = {"i": 0, "cls": 0, "guide": False}
    fig, (ax_xy, ax_zy) = plt.subplots(1, 2, figsize=(14, 6.5))
    fig.subplots_adjust(bottom=0.12, top=0.9)

    def title():
        e = events[st["i"]]
        nrev = sum(x["reviewed"] for x in events)
        fig.suptitle(f"event {st['i']+1}/{len(events)} (src {e['ev']})   "
                     f"paint class = [{CLASS[st['cls']][0]}]   {counts(e)}   "
                     f"reviewed {nrev}/{len(events)}   {'GUIDE' if st['guide'] else ''}")

    def draw():
        render(events[st["i"]], ax_xy, ax_zy, st["guide"])
        title(); fig.canvas.draw_idle()

    def apply(mask):
        e = events[st["i"]]
        e["label"][mask] = st["cls"]
        if st["cls"] in (0, 1):
            e["reviewed"] = True
        draw()

    def on_xy(verts):
        e = events[st["i"]]
        m = Path(verts).contains_points(e["xyz"][:, :2])
        if m.any(): apply(m)

    def on_zy(verts):
        e = events[st["i"]]
        m = Path(verts).contains_points(np.c_[e["xyz"][:, 2], e["xyz"][:, 1]])
        if m.any(): apply(m)

    ls1 = LassoSelector(ax_xy, on_xy)          # keep refs alive
    ls2 = LassoSelector(ax_zy, on_zy)

    def key(evt):
        k = (evt.key or "").lower()
        if k == "p": st["cls"] = 0; title(); fig.canvas.draw_idle()
        elif k == "r": st["cls"] = 1; title(); fig.canvas.draw_idle()
        elif k == "b": st["cls"] = 2; title(); fig.canvas.draw_idle()
        elif k == "g": st["guide"] = not st["guide"]; draw()
        elif k == "c":
            events[st["i"]]["label"][:] = 2; draw()
        elif k == "d":
            events[st["i"]]["reviewed"] = not events[st["i"]]["reviewed"]; title(); fig.canvas.draw_idle()
        elif k in ("n", "right"):
            save(events); st["i"] = min(st["i"] + 1, len(events) - 1); st["guide"] = False; draw()
        elif k in ("m", "left"):
            save(events); st["i"] = max(st["i"] - 1, 0); st["guide"] = False; draw()
        elif k == "s":
            print(f"saved, reviewed {save(events)}/{len(events)}")
        elif k == "q":
            print(f"quit, reviewed {save(events)}/{len(events)}"); plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", key)
    draw()
    print("labeler ready. Lasso a track, press P/R/B to paint. N next, S save, Q quit.")
    plt.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    smoke() if a.smoke else gui()
