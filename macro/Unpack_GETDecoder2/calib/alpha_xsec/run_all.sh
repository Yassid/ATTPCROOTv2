#!/bin/bash
# Re-track every P300_B0 run with the BEAM-OPTIONAL selection, so that the selection
# applied to the data matches the simulation eps_reco was measured on.
# Written once and NEVER edited while running -- bash re-reads a running script by byte
# offset, which silently truncates it mid-loop.
# NOTE: source the environment BEFORE set -u.  setup_fr19port.sh expands variables that
# may be unset, and with -u in force that kills the shell -- silently, because its stderr
# is redirected.  That is exactly how this script failed with a bare "exit 1" and no output.
cd /home/yassid/dec2014_package
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
set -u
W=/home/yassid/dec2014_retrack
echo "[$(date +%H:%M:%S)] start"
n=0
for f in "$W"/P300_B0/alpha_run_*_hits.root; do
  run=$(basename "$f" | sed 's/alpha_run_\([0-9]*\)_hits.root/\1/')
  csv="$W/csv/hits_$run.csv"
  pkl="$W/tracks/tracks_$run.pkl"
  if [ -s "$pkl" ]; then echo "  skip $run"; continue; fi
  [ -s "$csv" ] || root -l -b -q "lib/dump_hits.C(\"$f\",\"$csv\")" 2>&1 | grep -E "hits ->|ERROR"
  [ -s "$csv" ] || { echo "  run $run: NO CSV, skipping"; continue; }
  CSV="$csv" PKL="$pkl" RUN="$run" python3 -c '
import os, sys, pickle, importlib.util
sys.path.insert(0, "/home/yassid/dec2014_package")
spec = importlib.util.spec_from_file_location("tk", "/home/yassid/dec2014_package/02_track.py")
tk = importlib.util.module_from_spec(spec); spec.loader.exec_module(tk)
csv_path, pkl_path, run = os.environ["CSV"], os.environ["PKL"], os.environ["RUN"]
evs, groups = tk.load_csv(csv_path)
found = []
for e, h in zip(evs, groups):
    try:
        r = tk.find_tracks(h)
    except Exception:
        r = None
    if r:
        r["event"] = int(e)
        found.append(r)
pickle.dump({"nevents": len(evs), "cands": found}, open(pkl_path, "wb"))
pct = 100.0 * len(found) / max(len(evs), 1)
print("  run %s: %d events -> %d candidates (%.0f%%)" % (run, len(evs), len(found), pct))
'
  n=$((n+1))
  echo "[$(date +%H:%M:%S)] done $n"
done
echo "[$(date +%H:%M:%S)] ALL DONE: $(ls "$W"/tracks/*.pkl 2>/dev/null | wc -l) runs tracked"
