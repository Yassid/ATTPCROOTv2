#!/usr/bin/env python3
"""Vertex with the arms as the backbone and the beam as a refinement, not an equal partner.

The old vertex_from_tracks() gave beam and arms equal weight with no residual check, so a
beam line pointing somewhere the arms disagree with dragged the common point along the beam
and inflated the entrance->vertex path (the "pathological" events: beam line missing the
arms' intersection by ~150 mm).

Here the two arms define the vertex; the beam line is then tested against it. Under 10 mm --
the minimum of the clearly bimodal residual distribution -- the beam is consistent and is
folded into a 3-line fit, which constrains the vertex much better ALONG the beam, the one
direction two arms opening at 90 deg leave soft. Over 10 mm the beam is dropped and flagged.

Entrance: the tilt pivots at the window, so the beam enters on the detector axis and the
entrance is where the beam line comes closest to x=y=0. With no usable beam, the same line
is constructed through the vertex along the nominal tilt direction.
"""
import numpy as np, glob, pickle, warnings, sys
warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from cluster_tracks import fit_line, vertex_from_tracks

VD, T0 = 2.142, -34.0
BEAM_RES_CUT = 10.0
TILT, AZIM = np.radians(6.76), np.radians(-160.3)
NOMINAL = np.array([np.sin(TILT)*np.cos(AZIM), np.sin(TILT)*np.sin(AZIM), np.cos(TILT)])

def dist_line(p, c, d):
    w = p - c
    return float(np.linalg.norm(w - (w @ d) * d))

def axis_crossing(c, d):
    den = d[0]**2 + d[1]**2
    if den < 1e-12:
        return None
    return c - ((c[0]*d[0] + c[1]*d[1]) / den) * (-d) * -1.0 if False else c + (-(c[0]*d[0]+c[1]*d[1])/den) * d

def reconstruct(rec):
    T = (rec["tb"] - T0) * 160.0 * 1e-3
    P3 = np.column_stack([rec["x"], rec["y"], VD * T * 10])
    segs = {k: fit_line(P3[rec[k]], rec["q"][rec[k]]) for k in ("a1", "a2", "beam")}
    d1, d2 = segs["a1"][1], segs["a2"][1]
    op = float(np.degrees(np.arccos(abs(np.clip(d1 @ d2, -1, 1)))))
    v = vertex_from_tracks([segs["a1"], segs["a2"]])          # arms are the backbone
    res = dist_line(v, *segs["beam"])
    used_beam = res < BEAM_RES_CUT
    if used_beam:
        v = vertex_from_tracks([segs["a1"], segs["a2"], segs["beam"]])   # beam refines it
        ent = axis_crossing(*segs["beam"])
    else:
        ent = axis_crossing(v, NOMINAL)                        # nominal tilt through the vertex
    return dict(op=op, vtx=v, ent=ent, res=res, used_beam=used_beam, P3=P3, segs=segs)
