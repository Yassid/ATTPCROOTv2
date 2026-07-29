# a2091 15C(p,p') — UKF / GENFIT workspace

Modern staged reconstruction of the **a2091 15C proton-target runs**, ported from the
validated a1954 14C `UKF/` pipeline. Reaction: **15C(p,p')15C\*** (inverse kinematics,
15C beam on H2 active target at **300 torr**; detected ejectile = recoil PROTON).

Built by Claude on 2026-07-27 at the user's request ("analyze 15C+p, adapt everything to
the latest 12Be / 14C macros, gas pressure 300 torr; 15C+d comes later"). Source the env
first (see **Builds**), then run macros from this folder.

## Run set & data

- **Proton-target (15C+p) runs = 0138–0182** — 41 present on disk (see `runs_pp.txt`).
  Missing in that range: 0154 0155 0156 0168.
- **Deuteron-target (15C+d) runs = 0013–0133** — for the LATER (p,d) analysis, not touched yet.
- Raw HDF5: `/media/yassid/Seagate Hub/ATTPC/Data/a2091/run_XXXX.h5`
  (external Seagate drive; proton runs 0.35–37 GB, smallest = run_0179 @ 0.35 GB).
- Parameters: `parameters/ATTPC.a2091_C15.par` (B=2.85 T, H2 **300 torr**, drift 1.30 cm/µs,
  EField 45000 V/m, TBEntrance 300 — from the experimental ATTPC.a2091.par, with the full
  pad-plane geometry filled in from the a1954 template).
- Reco/fit outputs: `/home/yassid/a2091_C15_reco/` (+ `logs/`); slim cache + FRIB in
  `/home/yassid/a2091_C15_reco_slim/`; fits in `/home/yassid/a2091_C15_fit/`.

## Fitting plan (what the user asked for)

Fit 15C+p with **both** fitters, GENFIT **with and without** material effects — three passes
over the SAME gated inputs so they compare track-by-track:

| Pass | Driver | Fitter | matEffects | Notes |
|------|--------|--------|:----------:|-------|
| UKF          | `fitpipe_C15.sh`     | AtFitterUKFMulti | — (CATIMA ρ=3.308e-5) | fast, reference |
| GENFIT (noMat) | `fitpipe_C15.sh`   | AtGenfitter | **OFF** | helix, geometry = navigation only |
| GENFIT (matFX) | `genfit_mat_C15.sh` | AtGenfitter | **ON**  | dE/dx in the track model; uses `ATTPC_H300torr_RT` |

UKF (CATIMA) and GENFIT-matFX (geometry medium `H_300torr_RT`) are pinned to the **same**
ρ = 3.308e-5 g/cm³ so the only difference is the fitter, not the material.

## Pipeline

| Macro | Does | Output |
|-------|------|--------|
| `pipeline/unpackReco_C15.C(run, nEv, ...)` | unpack + **AtPSAMultiFit** + **AtDirDeDxCleaner** + PRA | `<run>_reco.root` (AtEventH + AtPatternEvent) |
| `pipeline/unpackFRIB_C15.C(run, ...)` | unpack the FRIB/IC group only | `<run>_FRIB.root` (small, for IC/PID gating) |
| `pipeline/slim_cache_C15.C` | strip reco.root to the AtPatternEvent cache | `<run>_slim.root` (fast local refits) |
| `pipeline/gate_events_C15.C` | IC-window + PID-polygon event gate | gated `<run>_reco.root` in `in/` |
| `pipeline/fitUKF_C15.C(run, nEv, "proton", -1)` | **UKF** (AtFitterUKFMulti) | `<run>_ukf.root` |
| `pipeline/fitGenfit_C15.C(run, ..., matEffects, ..., geoName)` | **GENFIT** (AtGenfitter) | `<run>_genfit.root` |
| `pp/ex_C15.C(runsCSV, inDir, Ebeam)` | two-body kinematics → 15C excitation spectrum | `pp/plots/ex_C15*.png` + kin cache |
| `pp/explore_C15.C` / `pp/make_explorer_html.C` | interactive / browser Ex explorer | X11 GUI / self-contained HTML |

Batch drivers (resumable, skip existing outputs; take a space-separated run list + parallelism):
- `./reco_hdb_C15.sh "$(cat runs_pp.txt)" 2`  — HDBSCAN reco (I/O-capped on the external drive)
- `./frib_C15_batch.sh "$(cat runs_pp.txt)" 4` — FRIB/IC unpack
- `./slim_cache_batch.sh "$(cat runs_pp.txt)" 4` — build slim caches
- `./fitpipe_C15.sh "$(cat runs_pp.txt)" 4`   — gate → UKF → GENFIT(noMat)
- `./genfit_mat_C15.sh "$(cat runs_pp.txt)" 4` — GENFIT(matFX) pass
- `./follow_C15.sh "$(cat runs_pp.txt)"`      — chase reco: slim+fit each run as reco finishes

## Builds (two, on purpose)

- **`build/`** — default build. Use for reco + UKF. `source build/config.sh`
- **`build_genfit/`** — GENFIT-enabled (`-DGENFIT=~/fair_install/GenFitInst`, Yassid/GenFit
  fork). Use for `fitGenfit_C15.C` / `genfit_mat_C15.sh`. `source build_genfit/config.sh`

Repo root is `/home/yassid/fair_install/ATTPCROOTv2` (the batch scripts set `REPO` to it).

## Handedness

Experimental data → **Bz = −2.85 T** (`bFieldSign = -1`, the default), same as 12Be/14C/a1975.

## Geometry

GENFIT material effects use **`ATTPC_H300torr_RT`** (medium `H_300torr_RT`, ρ = 3.308e-5 g/cm³,
room temperature). Generated on 2026-07-27 →
`geometry/ATTPC_H300torr_RT_geomanager.root`. Reco itself uses `ATTPC_H1bar.root` (UKF energy
loss is CATIMA/density-based, so reco geometry is navigation-only, as in the 14C pipeline).

## Flagged unknowns — CONFIRM before trusting absolute physics

1. **Beam energy** — `ex_C15.C` defaults `Ebeam = 195 MeV`, a PLACEHOLDER taken from the old
   a2091 `C15_pp_anaFit.C` (`Ebeam_buff = 195`; ~13 MeV/u for 15C). The absolute Ex scale
   depends on it — calibrate by putting the (p,p') elastic peak at Ex≈0.
2. **Gas density / temperature** — 3.308e-5 g/cm³ assumes room temperature. If the g.s. sits
   off zero after the beam-energy calibration, revisit the temperature (14C used 3.553e-5).
3. **IC / PID gates** — the IC window (500–900) and the deuteron/proton PID polygons are
   14C-era placeholders. Build a2091 gates from the sqrt(dEdx)-vs-Brho plane before physics.
4. **Space-charge / pad-time corrections** — OFF. The `resources/corrections/a1954` LUTs are
   referenced but not applied (`doSC=false`); no a2091 correction set exists yet.

## Typical commands

```bash
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source ../../../../build/config.sh          # from this folder: the UKF build

# smoke test: reco 200 events of the smallest proton run
root -b -q 'pipeline/unpackReco_C15.C("run_0179", 200)'
# UKF fit
root -b -q 'pipeline/fitUKF_C15.C("run_0179", -1, "proton", -1)'
# 15C excitation spectrum (SET Ebeam!)
root -b -q 'pp/ex_C15.C("run_0179", "/home/yassid/a2091_C15_reco/", 195.0)'

# GENFIT (needs the genfit build); matEffects OFF then ON
source ../../../../build_genfit/config.sh
root -b -q 'pipeline/fitGenfit_C15.C("run_0179")'                                   # noMat
root -b -q 'pipeline/fitGenfit_C15.C("run_0179",-1,"/home/yassid/a2091_C15_reco/","","",-2.85,2,5,"",4.0,10.0,170.0,kTRUE,kFALSE,"proton","ATTPC_H300torr_RT")'  # matFX
```
