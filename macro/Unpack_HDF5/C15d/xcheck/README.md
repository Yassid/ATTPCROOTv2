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
