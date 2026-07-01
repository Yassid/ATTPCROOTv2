#!/usr/bin/env python
"""Interactive side-by-side point-cloud viewer: Spyral (left) vs ATTPCROOT (right).
Self-contained HTML with an event slider — opens in any browser (no OpenGL needed,
works under WSL). Spyral clusters colored by HDBSCAN label; ATTPCROOT hits colored
by PRA track id. Faint gray = unclustered/background cloud each framework saw.

Run: ~/Spyral/venv/bin/python scripts/viewer.py
Out: plots/d2_pointcloud_viewer.html
"""
import h5py, json, glob, numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

WS = "/home/yassid/spyral_d2/workspace"
ATTPC = json.load(open("/home/yassid/spyral_d2/attpc_hits.json"))
SEL = json.load(open("/home/yassid/spyral_d2/viewer_events.json"))
OUT = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/spyral_compare/plots/d2_pointcloud_viewer.html"
CHUNK = 6935
PALETTE = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628","#f781bf","#00ced1"]

cat_of = {}
for c, evs in SEL["categories"].items():
    for e in evs: cat_of[e] = c

def spyral_event(ev):
    """Return (background_cloud Nx4, list of (label, cluster Nx4))."""
    run = 300 + ev // CHUNK
    bg = None; clusters = []
    for r in (run, run-1, run+1):
        try:
            fp = h5py.File(f"{WS}/Pointcloud/run_{r:04d}.h5", "r")
            if f"cloud_{ev}" in fp["cloud"]:
                bg = fp["cloud"][f"cloud_{ev}"][:][:, :4]
            fc = h5py.File(f"{WS}/Cluster/run_{r:04d}.h5", "r")
            if f"event_{ev}" in fc["cluster"]:
                g = fc["cluster"][f"event_{ev}"]
                for ck in g.keys():
                    cl = g[ck]["cloud"][:][:, :4]
                    clusters.append((int(g[ck].attrs["label"]), cl))
                break
        except (OSError, KeyError):
            continue
    return bg, clusters

events = SEL["events"]
fig = make_subplots(rows=1, cols=2, specs=[[{"type":"scene"},{"type":"scene"}]],
                    subplot_titles=("Spyral (HDBSCAN clusters)", "ATTPCROOT (PRA tracks)"),
                    horizontal_spacing=0.02)

trace_event = []   # which event each trace belongs to
titles = []
for ev in events:
    bg, clusters = spyral_event(ev)
    a = ATTPC.get(str(ev), {"hits": [], "tracks": []})
    # ---- Spyral (col 1) ----
    if bg is not None and len(bg):
        fig.add_trace(go.Scatter3d(x=bg[:,0], y=bg[:,1], z=bg[:,2], mode="markers",
            marker=dict(size=1.5, color="lightgray"), name="cloud", visible=False,
            hoverinfo="skip"), row=1, col=1); trace_event.append(ev)
    sp_th = []
    for i,(lab, cl) in enumerate(clusters):
        col = PALETTE[i % len(PALETTE)]
        fig.add_trace(go.Scatter3d(x=cl[:,0], y=cl[:,1], z=cl[:,2], mode="markers",
            marker=dict(size=2.5, color=col), name=f"sp_clu{lab}", visible=False,
            hovertemplate="x%{x:.0f} y%{y:.0f} z%{z:.0f}"), row=1, col=1); trace_event.append(ev)
    # ---- ATTPCROOT (col 2) ----
    hits = np.array(a["hits"]) if a["hits"] else np.zeros((0,4))
    if len(hits):
        fig.add_trace(go.Scatter3d(x=hits[:,0], y=hits[:,1], z=hits[:,2], mode="markers",
            marker=dict(size=1.5, color="lightgray"), name="hits", visible=False,
            hoverinfo="skip"), row=1, col=2); trace_event.append(ev)
    at_th=[]
    for i,tr in enumerate(a["tracks"]):
        pts = np.array(tr["pts"]) if tr["pts"] else np.zeros((0,4))
        if not len(pts): continue
        col = PALETTE[i % len(PALETTE)]
        at_th.append(round(tr["theta"],1))
        fig.add_trace(go.Scatter3d(x=pts[:,0], y=pts[:,1], z=pts[:,2], mode="markers",
            marker=dict(size=2.5, color=col), name=f"trk{tr['id']} th{tr['theta']:.0f}",
            visible=False, hovertemplate="x%{x:.0f} y%{y:.0f} z%{z:.0f}"), row=1, col=2); trace_event.append(ev)
    sp_th = [f"{np.degrees(np.arctan2(np.hypot(cl[:,0].mean(),cl[:,1].mean()), cl[:,2].mean())):.0f}" for _,cl in clusters]
    titles.append(f"event {ev}  [{cat_of.get(ev,'?')}]   Spyral: {len(clusters)} clusters   |   "
                  f"ATTPCROOT: {len([1 for t in a['tracks'] if t['pts']])} tracks  theta={at_th}")

# slider steps toggle visibility per event
trace_event = np.array(trace_event)
steps = []
for k, ev in enumerate(events):
    vis = (trace_event == ev).tolist()
    steps.append(dict(method="update", label=str(ev),
                      args=[{"visible": vis}, {"title.text": titles[k]}]))
# show first event
for j in range(len(fig.data)):
    fig.data[j].visible = (trace_event[j] == events[0])

fig.update_layout(
    sliders=[dict(active=0, currentvalue={"prefix":"event "}, steps=steps, pad={"t":40})],
    title=titles[0], height=720, showlegend=True, legend=dict(font=dict(size=9)),
    scene=dict(xaxis_title="x[mm]", yaxis_title="y[mm]", zaxis_title="z[mm] (drift)", aspectmode="data"),
    scene2=dict(xaxis_title="x[mm]", yaxis_title="y[mm]", zaxis_title="z[mm] (drift, native)", aspectmode="data"),
)
fig.write_html(OUT, include_plotlyjs="cdn")
print(f"saved {OUT}  ({len(fig.data)} traces, {len(events)} events)")
print("NOTE: z is reflected between the two frameworks (handedness) — compare shapes, not absolute z.")
