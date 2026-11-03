# ATTPCROOT vs Spyral cross-check: the azimuth/polar decomposition

Scripts behind the second half of `../CROSSCHECK_vs_SPYRAL.md`. Run order:

```bash
# 1. ATTPCROOT side: the merged dd tree has no phi, the per-run kin trees do
root -b -q dump_phi_to_csv.C                 # -> attpc_phi.csv  (410,723 tracks, 102 runs)

# 2. the geometry question: are the two pad maps the same object?
root -b -q dump_padmap.C                     # -> attpc_padxy.csv (10,240 pads)
python compare_padmaps.py                    # answer: yes, to 0.7 mm and 0.0000 deg

# 3. the physics questions (need the Spyral venv for polars)
V=/home/yassid/Data/spyral/venv_1.1.2/bin/python
$V theta_vs_azimuth.py    # the phi convention, and the tilt fit (tilt is excluded)
$V theta_harmonics.py     # harmonics n=1..8; n=6 is the SMALLEST, so not pad geometry
$V phi_vs_curvature.py    # dphi vs 1/pT, correlation +0.977 -> arc-length effect
$V theta_dependences.py   # dtheta vs theta, KE, chi2
```

The paths are absolute on purpose -- these read two analysis trees that live outside the repo
(`~/C15d_fit/run_*_kin_d.root` and
`~/Data/spyral/C15d_a2091/Workspace/InterpSolver/*.parquet`).

⚠ Spyral's parquets store the vertex in METRES; `spyral_kin.root` stores mm.

## The fitter decomposition (the part that resolved the theta offset)

```bash
# ATTPCROOT's SIMPLE estimator (AtSpyralPID), straight out of the PID caches
root -b -q dump_pid_polar.C            # -> attpc_pidpolar.csv (2,020,162 valid tracks)

V=/home/yassid/Data/spyral/venv_1.1.2/bin/python
$V compare_simple_estimators.py        # simple vs simple: agree to +-0.25 deg
$V theta_vs_vz.py                      # the offset is FLAT in vz: -0.091 +- 0.020 deg/m
```

`compare_simple_estimators.py` is the one that matters: AtSpyralPID is a port of Spyral's own
estimator, so agreement between the two codes' simple estimators while their FITTED thetas differ
by 0.45 deg places the offset in the fitters, not in the hits or the geometry.

### The Garfield on/off test

Not scripted here because it is a Spyral pipeline re-run, kept at
`/home/yassid/Data/spyral/gtest/` (`Run_on.py` / `Run_off.py`, differing only in
`do_garfield_correction`, phases `[True,True,True,False,False]`, runs 20-21, raw staged in
`gtest/raw`). Result: dTheta = -0.083 +- 0.005 deg.

⚠ Pair tracks between two Spyral passes ONLY via events where both found exactly ONE cluster.
`cluster_index` is a position, not an identity; the passes find different cluster counts per event,
and matching on it silently pairs different tracks (-0.249 deg, rms 1.08, and a 155 mm vertex jump
that gave it away -- against the correct -0.083, rms 0.055).
