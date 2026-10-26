"""Loop over the candidate files once and turn every event into a flat table.

Reconstruction is the slow part (a few minutes over all 68 runs), so the result is cached
to a .npz.  Delete the cache, or pass force=True, to rebuild it after changing V_DRIFT.
"""
import os, glob, pickle
import numpy as np
from .geometry import reconstruct
from .eloss import RangeTable


def build(config, max_files=None, force=False, quiet=False):
    tag = f"events_vd{config.V_DRIFT:.3f}".replace(".", "p")
    if max_files:
        tag += f"_n{max_files}"
    cache = os.path.join(config.OUT_DIR, tag + ".npz")
    if os.path.exists(cache) and not force:
        if not quiet:
            print(f"  using cached {os.path.basename(cache)}  (delete it to rebuild)")
        return dict(np.load(cache))

    tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)
    files = sorted(glob.glob(os.path.join(config.CANDS_DIR, "tracks_*.pkl"))) \
            or sorted(glob.glob(os.path.join(config.CANDS_DIR, "cands_*.pkl")))
    if max_files:
        files = files[:max_files]
    if not files:
        raise SystemExit(f"no track files in {config.CANDS_DIR} -- run ./01_dump.sh then 02_track.py")

    cols = {k: [] for k in ("run", "th1", "th2", "thsum", "opening",
                            "L1", "L2", "E1", "E2", "path", "vz", "nhit1", "nhit2",
                            "hasbeam", "event")}
    for i, path in enumerate(files):
        run = int(os.path.basename(path).split("_")[-1].split(".")[0])
        if not quiet:
            print(f"\r  reconstructing {i+1}/{len(files)}  run {run}   ", end="", flush=True)
        for rec in pickle.load(open(path, "rb"))["cands"]:
            r = reconstruct(rec, config.V_DRIFT, config.T0_BUCKETS, config.TB_NS)
            if r is None:
                continue
            L1, L2 = r["length"]
            cols["run"].append(run)
            cols["th1"].append(r["theta"][0]); cols["th2"].append(r["theta"][1])
            cols["thsum"].append(r["theta_sum"]); cols["opening"].append(r["opening"])
            cols["L1"].append(L1); cols["L2"].append(L2)
            cols["E1"].append(float(tab.energy_of_range(L1)))
            cols["E2"].append(float(tab.energy_of_range(L2)))
            cols["path"].append(r["path"]); cols["vz"].append(r["vertex"][2])
            cols["nhit1"].append(int(rec["a1"].sum())); cols["nhit2"].append(int(rec["a2"].sum()))
            # whether a beam track was actually found.  The vertex is better defined when
            # it was, so this is the handle for splitting the sample as a systematic.
            cols["hasbeam"].append(int(rec["beam"].sum() >= 2))
            cols["event"].append(int(rec.get("event", -1)))
    if not quiet:
        print()
    out = {k: np.asarray(v) for k, v in cols.items()}
    os.makedirs(config.OUT_DIR, exist_ok=True)
    np.savez_compressed(cache, **out)
    if not quiet:
        print(f"  {len(out['run'])} events -> {os.path.basename(cache)}")
    return out
