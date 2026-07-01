#!/usr/bin/env python
"""Convert Spyral Pointcloud + Cluster HDF5 (run_0300..0305) into a FLAT ROOT TTree
(one row per point) so the ROOT macro view_d2.C can draw Spyral clouds via TTree::Draw.

Branches (all scalar -> trivially readable in old FairSoft ROOT, no RNTuple):
  event/I  global event number (matches ATTPCROOT entry)
  kind/I   0 = clustered point, 1 = background (full point cloud)
  lab/I    HDBSCAN cluster label (kind==0); -1 for background
  x,y,z,q/F
Out: /home/yassid/spyral_d2/spyral_clouds.root
"""
import h5py, glob, numpy as np, uproot

WS = "/home/yassid/spyral_d2/workspace"
OUT = "/home/yassid/spyral_d2/spyral_clouds.root"

EV, KIND, LAB, X, Y, Z, Q = [], [], [], [], [], [], []

def add(ev, kind, lab, arr, xi=0, yi=1, zi=2, qi=3):
    n = arr.shape[0]
    EV.extend([ev]*n); KIND.extend([kind]*n); LAB.extend([lab]*n)
    X.extend(arr[:, xi]); Y.extend(arr[:, yi]); Z.extend(arr[:, zi]); Q.extend(arr[:, qi])

for cf in sorted(glob.glob(f"{WS}/Cluster/run_030*.h5")):
    run = cf.split("run_")[1][:4]
    fc = h5py.File(cf, "r")
    cloudp = h5py.File(f"{WS}/Pointcloud/run_{run}.h5", "r")["cloud"]
    for ekey in fc["cluster"].keys():
        ev = int(ekey.split("_")[1])
        for ck in fc["cluster"][ekey].keys():
            g = fc["cluster"][ekey][ck]
            add(ev, 0, int(g.attrs["label"]), g["cloud"][:])      # (N,5)
        pk = f"cloud_{ev}"
        if pk in cloudp:
            add(ev, 1, -1, cloudp[pk][:])                          # (N,8)

with uproot.recreate(OUT) as f:
    f.mktree("spyral", {"event": "int32", "kind": "int32", "lab": "int32",
                        "x": "float32", "y": "float32", "z": "float32", "q": "float32"})
    f["spyral"].extend({"event": np.array(EV, np.int32), "kind": np.array(KIND, np.int32),
                        "lab": np.array(LAB, np.int32), "x": np.array(X, np.float32),
                        "y": np.array(Y, np.float32), "z": np.array(Z, np.float32),
                        "q": np.array(Q, np.float32)})
print(f"wrote {len(EV)} points ({len(set(EV))} events) -> {OUT}")
