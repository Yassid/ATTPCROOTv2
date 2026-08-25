# a1954 12Be(p,p') — UKF / GENFIT workspace

Modern staged reconstruction of the **a1954 12Be proton-target runs**, ported from the
a1975 16C `UKF/` pipeline. Reaction: **12Be(p,p')12Be*** (inverse kinematics, 12Be beam
on H2 active target; detected ejectile = recoil PROTON).

Built by Claude on 2026-07-21. Source the env first (see **Builds** below), then run
macros from this folder.

## Run set & data

- **Proton-target physics runs = 0142–0215** (67 present; missing 0160, 0177–0181, 0183).
- Raw HDF5: `/media/yassid/NSCL_e15250/data/a1954_remerged/h5/run_XXXX.h5`
  (0142 is 0.6 GB; 0143–0215 are ~5–9 GB each).
- Parameters: `parameters/ATTPC.a1954_Be12.par` (B=2.85 T, H2 600 torr, drift 1.00 cm/µs).
- Reco/fit outputs: `/home/yassid/a1954_Be12_reco/` (+ `logs/`).

## Pipeline

| Macro | Does | Output |
|-------|------|--------|
| `pipeline/unpackReco_Be12.C(run, nEv, ...)` | unpack + **AtPSAMultiFit** + **AtDirDeDxCleaner** + PRA(triplclust) | `<run>_reco.root` (AtEventH + AtPatternEvent) |
| `pipeline/fitUKF_Be12.C(run, nEv, particles, bFieldSign)` | **new UKF** (AtFitterUKFMulti) — fast, reads `_reco.root` | `<run>_ukf.root` (AtTrackingEvent) |
| `pipeline/fitGenfit_Be12.C(run, nEv, ioDir, ...)` | **updated GENFIT** (AtGenfitter) — needs the genfit build | `<run>_genfit.root` (AtTrackingEvent + AtPIDEvent) |
| `pp/ex_Be12.C(runsCSV, inDir, Ebeam)` | two-body kinematics → 12Be excitation spectrum + KE-vs-θ | `pp/plots/ex_Be12.png` + `pp/plots/proton_kin<tag>.root` (cache) |
| `pp/explore_Be12.C(cache, ...)` | **interactive GUI** on that cache: Ebeam / cuts / all binnings, Ex recomputed live, `Zero g.s.` solves Ebeam | `pp/plots/explore_Be12.png` |

| `pp/make_explorer_html.C(cache, out, ...)` | bakes a cache into a standalone **browser** explorer (same controls, no X11) | `~/a1954_Be12_explorer.html` |

Browser explorer — for when WSLg/X11 misbehaves. Self-contained HTML, recomputes
E<sub>x</sub>/θ<sub>cm</sub> in JavaScript with the same two-body expressions:
```bash
root -b -q '<...>/pp/make_explorer_html.C()'      # bake cache -> ~/a1954_Be12_explorer.html
explorer.exe ~/a1954_Be12_explorer.html           # opens in the Windows browser
```
The g.s. fit is a Gauss-Newton χ² fit (Poisson errors) — verified equal to ROOT's
`TH1::Fit("gaus")`: μ +0.0195 vs +0.019, σ 0.2337 vs 0.234.
Note it defaults to θ<sub>lab</sub> ≥ 0, which drops 7 tracks with a negative fitted angle
that `ex_Be12.C` kept (26,597 vs 26,604).

Interactive explorer — plain ROOT is enough (no `config.sh`, no `VMCWORKDIR`; it finds
`plots/proton_kin_clean155.root` next to itself). Must run **without** `-b`:
```bash
root -l ~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954_Be12/UKF/pp/explore_Be12.C
```

Batch drivers (resumable, skip existing outputs):
- `./reco_batch.sh "run_0142 run_0147 ..." 2`   — reco (2-parallel; external drive I/O cap)
- `./fit_batch.sh  "run_0142 run_0147 ..." 4`    — UKF fit (CPU-bound, 4-parallel ok)

## Builds (two, on purpose)

- **`build/`** — the default build. GENFIT2 NOT FOUND. Use for UKF + reco.
  `source build/config.sh`
- **`build_genfit/`** — configured with `-DGENFIT=~/fair_install/GenFitInst` so
  `AtGenfitter` is compiled against the updated GenFit fork (Yassid/GenFit, SRIM
  energy-loss tables + `2584bfe` dEdxParam guard). Use for `fitGenfit_Be12.C`.
  `source build_genfit/config.sh`

Both build from the same source; kept separate so a GENFIT rebuild never disturbs a
running reco/UKF batch.

## Handedness — VALIDATED

Experimental data has the opposite helix handedness to simulation, so **Bz = −2.85 T**
(`bFieldSign = -1`, the default). Confirmed on run_0142 (proton hypothesis):

| Bz | median χ²/ndf | median KE | median θ |
|----|-----|-----|-----|
| **−2.85 (exp)** | **0.03** | **5.96 MeV** | 63° |

Clean elastic-recoil locus — matches the a1975 result.

## Flagged unknowns (confirm before trusting absolute physics)

1. **Beam energy** — `ex_Be12.C` defaults `Ebeam = 100 MeV` as a PLACEHOLDER. The
   absolute Ex scale depends on it. Set the real a1954 12Be beam energy.
2. ~~**Gas density for CATIMA**~~ — **RESOLVED 2026-08-25: a1954 ran H2 at 300 torr,
   rho = 3.308e-5 g/cm³, geometry `ATTPC_H300torr_RT`** (`media.geo` `H_300torr_RT`).
   `fitUKF_Be12.C`, `fitGenfit_Be12.C` and all drivers now default to it, and
   `ATTPC.a1954_Be12.par` reads `GasPressure 300 / Density 0.0331`.
   Everything produced before that date used 6.5e-5 / `ATTPC_H600torr` — twice the
   material. Live for the UKF (its eloss correction uses the density) and for the
   (p,d) genfit set (`matEffects=kTRUE`); inert only for the (p,p') genfit set, which
   ran `matEffects=kFALSE`.
3. **Pad-time correction** — OFF (`applyTimeCorr=false`); the a1975 CSV is experiment-
   specific. Measure an a1954 pad-time table if needed.
4. **No IC beam gate / no proton PID gate yet** — `ex_Be12.C` uses a fit-quality cut
   only. Add: (a) FRIB IC gate via the `_FRIBDAQ` unpacker → `<run>_FRIB.root`;
   (b) a 12Be proton gate from the sqrt(dEdx)-vs-Brho plane (AtSpyralPID).

## Typical commands

```bash
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source ../../../../build/config.sh          # from this folder: the UKF build

# reco one run (slow, once)
root -b -q 'pipeline/unpackReco_Be12.C("run_0142")'
# UKF fit (fast, iterate)
root -b -q 'pipeline/fitUKF_Be12.C("run_0142", -1, "proton", -1)'
# 12Be excitation spectrum (SET Ebeam!)
root -b -q 'pp/ex_Be12.C("run_0142", "/home/yassid/a1954_Be12_reco/", 100.0)'

# GENFIT fit (needs the genfit build)
source ../../../../build_genfit/config.sh
root -b -q 'pipeline/fitGenfit_Be12.C("run_0142")'
```
