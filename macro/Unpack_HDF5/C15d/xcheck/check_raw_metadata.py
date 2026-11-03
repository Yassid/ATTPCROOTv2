import h5py, sys
for run in sys.argv[1:]:
    p=f"/media/yassid/Seagate Hub/ATTPC/Data/a2091/run_{run}.h5"
    try:
        with h5py.File(p,'r') as f:
            ok_get=True; n=None
            try: n=len(f['get'].keys())
            except Exception as e: ok_get=False; err=f"{type(e).__name__}: {e}"
            meta=list(f['meta']['meta'][:]) if 'meta' in f and 'meta' in f['meta'] else None
            print(f"run_{run}: /get {'OK '+str(n)+' datasets' if ok_get else 'DAMAGED -> '+err}; meta={meta}", flush=True)
    except Exception as e:
        print(f"run_{run}: CANNOT OPEN {type(e).__name__}: {e}", flush=True)
