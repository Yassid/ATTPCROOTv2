#!/usr/bin/env python
"""Manual point-cloud labeler for real a1975 16C(d,p) events (GNN gold set / training labels).

Label individual particle TRACKS as instances; beam + noise haze = background. Trains track
SEPARATION (species ID left to genfit/PID). Labels: 0 = background, 1,2,3,... = tracks.

** KEY FEATURE: each unreviewed event is PRE-SEEDED with dircluster (direction+dE/dx clustering) **
so you CORRECT rather than label from scratch. dircluster rarely merges two tracks (so you almost
never have to split) and only over-segments -> your job is mostly to MERGE fragments of the same
track and confirm. Noise is pre-marked as background.

Run (GUI, needs WSLg display + PySide6 in gnn_env):
    ~/gnn_env/bin/python label_tool.py [--input data/label_input.parquet] [--limit N]
Smoke test (no GUI, pre-seeds event 0 with dircluster and renders it):
    ~/gnn_env/bin/python label_tool.py --smoke

Workflow per event: it opens PRE-SEEDED (colored tracks + gray background). To merge fragments of
one physical track: press its number (1..9), then lasso the other fragments -> they join it.
To fix noise/beam grabbed into a track: press B, lasso it. Press D to confirm (reviewed), N next.

Keys:
  (lasso)  paint hits with the CURRENT track id      T / Enter  start a NEW track (next id)
  1..9     select an existing track id (merge/fix)   B / 0      paint selection = background
  G        re-seed this event with dircluster        C          reset event to all background
  Z        undo last lasso                            D          toggle 'reviewed' flag
  S        save now      N/-> next (autosave)         M/<-       previous event      Q quit
Mouse: drag to lasso on either panel. Point size scales with charge.
Output: data/labels.parquet (event,x,y,z,q,label,reviewed) -- only reviewed=True events train.
"""
import argparse, os, sys, numpy as np, pandas as pd
from matplotlib.path import Path

BG_COLOR = "#c9c9c9"
IN_PARQUET  = "data/label_input.parquet"
OUT_PARQUET = "data/labels.parquet"

# import dircluster from the parent dir (16C_dp_gnn/) for the pre-seed
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
try:
    from dircluster import cluster as _dircluster
except Exception:
    _dircluster = None


def track_palette():
    import matplotlib.cm as cm
    return cm.tab20(np.linspace(0, 1, 20))


def color_for(lbls, pal):
    out = np.empty((len(lbls), 4))
    for i, l in enumerate(lbls):
        out[i] = (0.79, 0.79, 0.79, 1.0) if l <= 0 else pal[(int(l) - 1) % 20]
    return out


def preseed(xyz, q):
    """dircluster -> initial track labels (fragments 1,2,...; noise -> 0=background)."""
    if _dircluster is None:
        return np.zeros(len(q), np.int16)
    lab = _dircluster(xyz, q, qratio=0.65, min_hits=4)   # noise=-1
    out = np.zeros(len(lab), np.int16)
    for newid, c in enumerate(sorted(set(lab[lab >= 0])), start=1):
        out[lab == c] = newid
    return out


def compact(lbl):
    """renumber positive track ids to 1..K contiguous (clean output)."""
    out = np.zeros_like(lbl)
    for newid, t in enumerate(sorted(set(lbl[lbl > 0].tolist())), start=1):
        out[lbl == t] = newid
    return out


def load(in_parquet, limit=None):
    df = pd.read_parquet(in_parquet)
    events = []
    for ev, g in df.groupby("event", sort=True):
        g = g.reset_index(drop=True)
        events.append({"ev": int(ev), "xyz": g[["x", "y", "z"]].to_numpy(float),
                       "q": g["q"].to_numpy(float),
                       "label": np.zeros(len(g), np.int16), "reviewed": False, "seeded": False})
        if limit and len(events) >= limit:
            break
    try:
        prev = pd.read_parquet(OUT_PARQUET)
        for e in events:
            p = prev[prev.event == e["ev"]]
            if len(p) == len(e["label"]):
                e["label"] = p["label"].to_numpy(np.int16); e["reviewed"] = bool(p["reviewed"].iloc[0])
                # keep any event that already has work (reviewed OR any track label); only re-seed
                # truly-untouched (all-background) events -> never wipes un-confirmed progress
                e["seeded"] = e["reviewed"] or bool((e["label"] > 0).any())
        print(f"resumed: {int(prev.groupby('event').reviewed.first().sum())} events already reviewed")
    except FileNotFoundError:
        pass
    return events


def save(events):
    frames = [pd.DataFrame({"event": e["ev"], "x": e["xyz"][:, 0], "y": e["xyz"][:, 1],
                            "z": e["xyz"][:, 2], "q": e["q"], "label": compact(e["label"]),
                            "reviewed": e["reviewed"]}) for e in events]
    out = pd.concat(frames, ignore_index=True)
    out.to_parquet(OUT_PARQUET, index=False)
    return int(out.groupby("event").reviewed.first().sum())


def sizes(q):
    s = np.sqrt(np.clip(q, 1, None))
    return 6 + 34 * (s - s.min()) / (np.ptp(s) + 1e-9)


def render(e, ax_xy, ax_zy, pal):
    xyz, q = e["xyz"], e["q"]; cols = color_for(e["label"], pal); sz = sizes(q)
    ax_xy.clear(); ax_zy.clear()
    ax_xy.scatter(xyz[:, 0], xyz[:, 1], s=sz, c=cols, edgecolors="none")
    ax_zy.scatter(xyz[:, 2], xyz[:, 1], s=sz, c=cols, edgecolors="none")
    ax_xy.set_title("pad plane  X vs Y"); ax_xy.set_xlabel("x [mm]"); ax_xy.set_ylabel("y [mm]")
    ax_xy.set_aspect("equal", "box")
    ax_zy.set_title("drift  Z vs Y"); ax_zy.set_xlabel("z [mm]"); ax_zy.set_ylabel("y [mm]")


def nspirals(e):
    return int((np.unique(e["label"]) > 0).sum())


def free_id(e):
    return int(e["label"].max()) + 1 if e["label"].max() > 0 else 1


def ensure_seeded(e):
    if not e["seeded"]:
        e["label"] = preseed(e["xyz"], e["q"]); e["seeded"] = True


def smoke(in_parquet):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    events = load(in_parquet); pal = track_palette()
    e = next((ev for ev in events if not ev["reviewed"]), events[0])   # test pre-seed on an UNreviewed event
    ensure_seeded(e)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.5)); render(e, a1, a2, pal)
    fig.suptitle(f"SMOKE (dircluster pre-seed) event {e['ev']}  n={len(e['q'])}  tracks={nspirals(e)}")
    plt.tight_layout(); plt.savefig("data/smoke_event0.png", dpi=110)
    print(f"smoke OK: {len(events)} events, dircluster pre-seed -> {nspirals(e)} tracks, rendered data/smoke_event0.png")


def gui(in_parquet, limit):
    import matplotlib; matplotlib.use("QtAgg")
    import matplotlib.pyplot as plt
    from matplotlib.widgets import LassoSelector
    events = load(in_parquet, limit); pal = track_palette()
    st = {"i": 0, "cur": 1, "undo": None}
    fig, (ax_xy, ax_zy) = plt.subplots(1, 2, figsize=(14, 6.5)); fig.subplots_adjust(bottom=0.1, top=0.9)

    def title():
        e = events[st["i"]]; nrev = sum(x["reviewed"] for x in events)
        cur = "background" if st["cur"] == 0 else f"track #{st['cur']}"
        fig.suptitle(f"event {st['i']+1}/{len(events)} (src {e['ev']})   painting = [{cur}]   "
                     f"tracks here: {nspirals(e)}   reviewed {nrev}/{len(events)}"
                     f"   {'[reviewed]' if e['reviewed'] else ''}")

    def draw():
        render(events[st["i"]], ax_xy, ax_zy, pal); title(); fig.canvas.draw_idle()

    def apply(mask):
        e = events[st["i"]]; st["undo"] = (st["i"], e["label"].copy())
        e["label"][mask] = st["cur"]
        if st["cur"] >= 1: e["reviewed"] = True
        draw()

    def on_xy(v):
        e = events[st["i"]]; m = Path(v).contains_points(e["xyz"][:, :2])
        if m.any(): apply(m)

    def on_zy(v):
        e = events[st["i"]]; m = Path(v).contains_points(np.c_[e["xyz"][:, 2], e["xyz"][:, 1]])
        if m.any(): apply(m)

    ls1 = LassoSelector(ax_xy, on_xy); ls2 = LassoSelector(ax_zy, on_zy)  # keep refs

    def goto(j):
        save(events); st["i"] = min(max(j, 0), len(events)-1); ensure_seeded(events[st["i"]])
        st["cur"] = free_id(events[st["i"]]); st["undo"] = None; draw()

    def key(evt):
        k = (evt.key or "").lower(); e = events[st["i"]]
        if k in ("t", "enter"): st["cur"] = free_id(e); title(); fig.canvas.draw_idle()
        elif k in ("b", "0"): st["cur"] = 0; title(); fig.canvas.draw_idle()
        elif k in list("123456789"): st["cur"] = int(k); title(); fig.canvas.draw_idle()
        elif k == "g": e["seeded"] = False; ensure_seeded(e); e["reviewed"] = False; draw()
        elif k == "z" and st["undo"] and st["undo"][0] == st["i"]: e["label"] = st["undo"][1]; st["undo"] = None; draw()
        elif k == "c": e["label"][:] = 0; e["seeded"] = True; st["cur"] = 1; draw()
        elif k == "d": e["reviewed"] = not e["reviewed"]; title(); fig.canvas.draw_idle()
        elif k in ("n", "right"): goto(st["i"]+1)
        elif k in ("m", "left"): goto(st["i"]-1)
        elif k == "s": print(f"saved, reviewed {save(events)}/{len(events)}")
        elif k == "q": print(f"quit, reviewed {save(events)}/{len(events)}"); plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", key)
    ensure_seeded(events[0]); st["cur"] = free_id(events[0]); draw()
    print("labeler ready (dircluster pre-seed). Merge fragments: press track#, lasso them. "
          "T=new, B=bg, Z=undo, G=re-seed, D=confirm, N=next, Q=quit.")
    plt.show()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--input", default=IN_PARQUET)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    smoke(a.input) if a.smoke else gui(a.input, a.limit)
