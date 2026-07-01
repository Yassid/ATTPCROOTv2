#!/usr/bin/env python3
"""Split one a1975 legacy-merger HDF5 run into N chunk-"runs" for parallel Spyral.

Spyral parallelizes ACROSS runs (n_processes), not within a run, so a single big
run is single-process. This splits run_XXXX.h5 into N valid chunk files, each a
self-contained "run" with its own meta/meta event range, so Spyral can process
them n-way parallel and we merge the per-chunk parquet afterwards.

Legacy layout copied verbatim per event i in the chunk's range:
  get/evt{i}_*           (evt{i}_data, evt{i}_header)
  frib/evt/evt{i}_*      (evt{i}_1903, evt{i}_header, [evt{i}_977])
Run-level groups copied to every chunk:
  meta/*  (cobo*asad* datasets; meta/meta rewritten to the chunk's [min,ts,max,ts])
  frib/{scaler,runinfo,title}

Usage:
  python chunk_run.py --src /mnt/f/a1975/h5/run_0016.h5 \
                      --outdir /home/yassid/spyral_d2/h5 \
                      --out-base 300 --nchunks 6
  -> writes run_0300.h5 .. run_0305.h5  (set Spyral run_min=300 run_max=305)
"""
import argparse
import h5py
import numpy as np
from pathlib import Path


def copy_event(src, dst, group_path, evt):
    """Copy all datasets named evt{evt}_* from src[group_path] to dst[group_path]."""
    sg = src[group_path]
    dg = dst.require_group(group_path)
    prefix = f"evt{evt}_"
    for name in sg.keys():
        if name.startswith(prefix):
            sg.copy(name, dg, name=name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--out-base", type=int, required=True,
                    help="first chunk run number, e.g. 300 -> run_0300.h5 ...")
    ap.add_argument("--nchunks", type=int, default=6)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.src, "r") as src:
        meta = src["meta"]["meta"][:]
        ev_min, ev_max = int(meta[0]), int(meta[2])
        ts0, ts1 = meta[1], meta[3]
        n_events = ev_max - ev_min + 1
        edges = np.linspace(ev_min, ev_max + 1, args.nchunks + 1).astype(int)
        print(f"src {args.src}: events {ev_min}..{ev_max} ({n_events}); "
              f"{args.nchunks} chunks, edges={edges.tolist()}")

        for c in range(args.nchunks):
            lo, hi = int(edges[c]), int(edges[c + 1]) - 1  # inclusive
            run_no = args.out_base + c
            out = outdir / f"run_{run_no:04d}.h5"
            with h5py.File(out, "w") as dst:
                # run-level meta group: copy cobo*asad* datasets, rewrite meta/meta
                mg = dst.require_group("meta")
                for name in src["meta"].keys():
                    if name == "meta":
                        continue
                    src["meta"].copy(name, mg, name=name)
                mg.create_dataset("meta", data=np.array([lo, ts0, hi, ts1],
                                                        dtype=meta.dtype))
                # run-level frib non-event datasets/groups
                fr = dst.require_group("frib")
                for name in ("scaler", "runinfo", "title"):
                    if name in src["frib"]:
                        src["frib"].copy(name, fr, name=name)
                dst.require_group("get")
                dst.require_group("frib/evt")
                # per-event payload
                for ev in range(lo, hi + 1):
                    copy_event(src, dst, "get", ev)
                    copy_event(src, dst, "frib/evt", ev)
            print(f"  wrote {out.name}: events {lo}..{hi} ({hi - lo + 1})")


if __name__ == "__main__":
    main()
