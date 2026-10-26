#!/usr/bin/env python3
"""Step 0 -- check the setup before anything else.  Run this first.

    python3 00_check.py
"""
import sys, os, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

print("Python", sys.version.split()[0])
missing = []
for mod in ("numpy", "scipy", "matplotlib"):
    try:
        __import__(mod)
        print(f"  OK   {mod}")
    except ImportError:
        print(f"  MISSING {mod}")
        missing.append(mod)
if missing:
    print(f"\ninstall with:  pip3 install --user {' '.join(missing)}")
    sys.exit(1)

import config
print(f"\nDATA_ROOT = {config.DATA_ROOT}")
if not os.path.isdir(config.DATA_ROOT):
    print("  NOT FOUND.  Edit DATA_ROOT in config.py, or set the environment variable:")
    print("     export DEC2014_DATA=/path/to/dec2014_data")
    sys.exit(1)

rootf = sorted(glob.glob(os.path.join(config.ROOT_DIR, "alpha_run_*_hits.root")))
print(f"ROOT_DIR  = {config.ROOT_DIR}")
print(f"  {len(rootf)} unpacked ROOT files")
import shutil
print(f"  ROOT on PATH: {'yes, ' + os.popen('root-config --version 2>/dev/null').read().strip() if shutil.which('root') else 'NO -- source thisroot.sh before step 1'}")
files = sorted(glob.glob(os.path.join(config.TRACKS_DIR, "tracks_*.pkl")))
print(f"TRACKS_DIR= {config.TRACKS_DIR}")
print(f"  {len(files)} tracked runs found" + ("" if files else "   <-- nothing here, check the path"))
if not files:
    print("     (expected -- step 2 has not run yet)")
    for mod in ("sklearn",):
        try:
            __import__(mod); print(f"  OK   {mod}  (needed for tracking)")
        except ImportError:
            print(f"  MISSING {mod}  -- pip3 install --user scikit-learn==1.3.2")
    print("\nSetup is good.  Next:  ./01_dump.sh   then   python3 02_track.py")
    sys.exit(0)

if not os.path.isfile(config.ELOSS):
    print(f"  range table MISSING: {config.ELOSS}"); sys.exit(1)

import pickle
d = pickle.load(open(files[0], "rb"))
rec = d["cands"][0]
print(f"  first file holds {len(d['cands'])} candidates from {d['nevents']} events")
print(f"  a candidate has keys: {sorted(rec.keys())}")
print(f"  event 0 has {len(rec['tb'])} hits; "
      f"beam={rec['beam'].sum()}  a1={rec['a1'].sum()}  a2={rec['a2'].sum()}")

from lib.eloss import RangeTable
t = RangeTable(config.ELOSS, config.PRESSURE)
print(f"\nrange table OK: an {config.E_BEAM} MeV alpha travels "
      f"{t.total_mm:.0f} mm at {config.PRESSURE} torr")
os.makedirs(config.OUT_DIR, exist_ok=True)
print(f"output will go to {config.OUT_DIR}")
print("\nSetup is good.  Next:  python3 03_kinematics.py")
