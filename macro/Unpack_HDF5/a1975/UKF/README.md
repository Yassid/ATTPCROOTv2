# a1975 16C+p — new-UKF workspace

Incremental pipeline for the a1975 16C+p (proton-target) experimental data, using
the new UKF (`EventFit::AtFitterUKF` + OpenKF). Raw data on the F: drive
(`/mnt/f/a1975/h5/`); 10-run reco/PID production on `/mnt/f/a1975/reco/` (runs
0106–0115). Source the env first: `source build/config.sh` from the repo root.

**Run all macros from this folder** (the UKF root). Macros read bare-name data
(`run_0116_reco.root`, …) from here and write plots into each folder's `plots/`
subfolder, so paths resolve relative to this directory.

## Directory layout (organized by reaction / role)

```
UKF/
├── README.md
├── <data>.root              unique LOCAL data — run_0116/0117/0118 reco/FRIB/pid/ukf/…
│                            (NOT on F:; kept here, read by bare name) + small caches
│                            (proton_kin.root, combined_pidobs.root, *_pidobs.root)
│
├── pipeline/                reaction-agnostic unpack → fit → PID-run → display infra
│   ├── unpack_a1975_UKF.C, unpackPSA_, unpackFRIB_, unpackReco_, unpackNFit_
│   ├── fitUKF_a1975.C, runPID_a1975.C, discover_aux_a1975.C
│   ├── event_png_UKF.C, event_png_batch_UKF.C, run_eve_UKF.C
│   ├── ukf_scan.sh, ukf_scan_mom.sh
│   └── plots/               event displays
│
├── pid/                     PID port, planes, gates, leakage (shared p/d separation)
│   ├── pid_spyral_a1975.C   faithful Spyral PID plane (the trustworthy estimator)
│   ├── pid_stats_a1975.C    high-stats old-vs-faithful side-by-side
│   ├── pid_clean/_plane/_ic/_gates/_compare_methods/_a1975_UKF.C
│   ├── make_gates_, make_proton_gate_, build_proton_gate_, gate_proton_, leak_a1975.C
│   ├── proton_band.json     active proton gate (sqrtdEdx vs brho)
│   └── plots/
│
├── pp/                      16C(p,p′) — proton elastic / inelastic
│   ├── ex_a1975.C           16C excitation-energy spectrum (two-body kinematics)
│   ├── ex_eval.C, ex_overlay.C
│   ├── ukf_clean_a1975.C, ukf_results_a1975.C, vertex_correct_a1975.C
│   └── plots/
│
└── pd/                      16C(p,d)15C — neutron pickup
    ├── ex_pd_a1975.C        15C excitation-energy spectrum from DEUTERON-hyp fits
    └── plots/
```

Gate defaults point to `pid/proton_band.json`; data products live at the UKF root
and are read by bare name. New plots go to `<folder>/plots/`.

## Pipeline stages (`pipeline/`)

| Macro | Does | Output |
|-------|------|--------|
| `pipeline/unpack_a1975_UKF.C(run, nEv)` | unpack only | `<run>_unpack.root` (AtRawEvent) |
| `pipeline/unpackPSA_a1975_UKF.C(run, nEv, persistRaw)` | unpack + PSA | `<run>_psa.root` / `<run>_disp.root` |
| `pipeline/unpackReco_a1975_UKF.C(run, nEv, persistRaw)` | unpack + PSA + **SC + PRA** | `<run>_reco.root` (AtPatternEvent) |
| `pipeline/fitUKF_a1975.C(run, nEv, particles, bFieldSign, ...)` | **UKF only** (reads `<run>_reco.root`) — fast | `<run>_ukf<suffix>.root` (AtTrackingEvent) |
| `pipeline/unpackNFit_a1975_UKF.C(run, nEv, particles, bFieldSign, ...)` | **full end-to-end** | `<run>_ukf_full.root` |

Event display (WSL has no working OpenGL, so use the static one):
- `pipeline/event_png_UKF.C(file, evt)` / `pipeline/event_png_batch_UKF.C(file, evStart, evEnd)` — PNG XY/XZ/YZ + 3D → `pipeline/plots/`.
- `pipeline/run_eve_UKF.C` — interactive Eve viewer (does NOT work on this WSL; OpenGL).

## Particle selection (macro level)

`particles` is a comma-separated list, e.g. `"proton"` or `"proton,deuteron"`.
Supported: proton, deuteron, triton, He3, alpha. Each becomes one UKF hypothesis
(`AtFitterUKFMulti` keeps the best χ²/ndf per track). The `outSuffix` arg of
`fitUKF_a1975.C` writes to a separate `<run>_ukf<suffix>.root` (e.g. `_d` for a
deuteron-hypothesis refit) so different hypotheses don't clobber each other.

## Reactions

- **16C(p,p′)** (`pp/`): proton ejectile → 16C excitation spectrum. `pp/ex_a1975.C`
  with two-body kinematics m1=m4=16C, m2=m3=p. Elastic → Ex≈0.
- **16C(p,d)15C** (`pd/`): deuteron ejectile (neutron pickup) → 15C excitation
  spectrum. `pd/ex_pd_a1975.C` with m1=16C, m2=p, m3=d, m4=15C; needs a
  **deuteron-hypothesis** UKF fit (`fitUKF_a1975.C(...,"deuteron",...,"_d",...)`).
  The dEdx PID alone does NOT isolate (p,d) — the decisive discriminator is the
  deuteron-hypothesis kinematics closing onto a sharp 15C ground-state peak.

## Orientation / handedness (IMPORTANT)

Experimental data has the **opposite helix handedness** to simulation (sim
reverses drift z in digitization, `AtClusterize.cxx:128`; experiment does not).
So the UKF uses **Bz = −2.85 T** for exp data. Controlled by `bFieldSign`:
`-1` = experimental (default), `+1` = simulation convention.

Empirically confirmed on run_0116:

| Bz | χ²/ndf (med) | KE (med) | physical KE |
|----|-----|-----|-----|
| **−2.85 (exp)** | **0.04** | **7.9 MeV** | 75% |
| +2.85 (sim) | 3.09 | 152 MeV | 34% |

KE-vs-θ shows a clean elastic recoil locus at θ≈70–90°, KE 2–7 MeV.

## Resolution-tuned UKF parameters

From the Ex-resolution scans (elastic-peak FWHM metric), the best parameters are
**`measSigma = 0.5 mm`, `momSigmaFrac = 0.1`** (~26% better FWHM than the 2.0/0.3
default). PID (brho/dEdx) is charge-sign agnostic, so it is unaffected by the
handedness flip.

Remaining notes:
- `gasDensity = 9.0e-5 g/cm³` (H₂ at 1 bar) is an estimate — verify against run conditions.
- A few diverged fits give huge KE outliers — filter on `AtFitMetadata` convergence / a χ²/ndf cut.
- `chargeSign` filtering is OFF in the multi-fitter (PRA sign is sim-calibrated, unreliable for exp until validated).

## Typical commands (run from this folder)

```bash
# stage to track candidates (slow, once per run)
root -b -q 'pipeline/unpackReco_a1975_UKF.C("run_0116", 1000, true)'

# iterate the UKF — proton hypothesis (fast)
root -b -q 'pipeline/fitUKF_a1975.C("run_0116", -1, "proton", -1)'

# deuteron-hypothesis refit for the (p,d) channel (separate output file)
root -b -q 'pipeline/fitUKF_a1975.C("run_0106", -1, "deuteron", -1, 2.85, 9.0e-5, "_d", "/mnt/f/a1975/reco/", 0.5, 0.1, 1, 10)'

# PID quality plane (faithful Spyral estimator) → pid/plots/
root -b -q 'pid/pid_spyral_a1975.C("run_0116")'

# 16C(p,p') excitation spectrum → pp/plots/
root -b -q 'pp/ex_a1975.C("run_0106,run_0107,...", "/mnt/f/a1975/reco/")'

# 16C(p,d) 15C excitation spectrum (needs <run>_ukf_d.root) → pd/plots/
root -b -q 'pd/ex_pd_a1975.C("run_0106", "/mnt/f/a1975/reco/")'
```
