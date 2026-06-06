# a1975 16C+p — new-UKF workspace

Incremental pipeline for the a1975 16C+p (proton-target) experimental data, using
the new UKF (`EventFit::AtFitterUKF` + OpenKF). Raw data on the F: drive
(`/mnt/f/a1975/h5/`). Source the env first: `source build/config.sh` from the repo root.

## Pipeline stages (macros)

| Macro | Does | Output |
|-------|------|--------|
| `unpack_a1975_UKF.C(run, nEv)` | unpack only | `<run>_unpack.root` (AtRawEvent) |
| `unpackPSA_a1975_UKF.C(run, nEv, persistRaw)` | unpack + PSA | `<run>_psa.root` / `<run>_disp.root` (AtEventH [+AtRawEvent]) |
| `unpackReco_a1975_UKF.C(run, nEv, persistRaw)` | unpack + PSA + **SC + PRA** | `<run>_reco.root` (AtPatternEvent) |
| `fitUKF_a1975.C(run, nEv, particles, bFieldSign, ...)` | **UKF only** (reads `<run>_reco.root`) — fast | `<run>_ukf.root` (AtTrackingEvent) |
| `unpackNFit_a1975_UKF.C(run, nEv, particles, bFieldSign, ...)` | **full end-to-end** | `<run>_ukf_full.root` |

Event display (WSL has no working OpenGL, so use the static one):
- `event_png_UKF.C(file, evt)` / `event_png_batch_UKF.C(file, evStart, evEnd)` — PNG XY/XZ/YZ + 3D.
- `run_eve_UKF.C` — interactive Eve viewer (does NOT work on this WSL; OpenGL).

## Particle selection (macro level)

`particles` is a comma-separated list, e.g. `"proton"` or `"proton,deuteron"`.
Supported: proton, deuteron, triton, He3, alpha. Each becomes one UKF hypothesis
(`AtFitterUKFMulti` keeps the best χ²/ndf per track).

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

## Known tuning TODO

- χ²/ndf piles near 0 ⇒ `SetMeasurementSigma(2.0)` is too loose; tighten to ~0.5–1 mm
  for χ²/ndf≈1 and meaningful uncertainties.
- A few diverged fits give huge KE outliers (~1e7 MeV) — filter on `AtFitMetadata`
  convergence in analysis.
- `gasDensity = 9.0e-5 g/cm³` (H₂ at 1 bar) is an estimate — verify against run conditions.
- `chargeSign` filtering is OFF in the multi-fitter (PRA sign is sim-calibrated,
  unreliable for exp until validated).

## Typical commands

```bash
# stage to track candidates (slow, once per run)
root -b -q 'unpackReco_a1975_UKF.C("run_0116", 1000, true)'
# iterate the UKF (fast)
root -b -q 'fitUKF_a1975.C("run_0116", 1000, "proton", -1)'
# or full end-to-end
root -b -q 'unpackNFit_a1975_UKF.C("run_0116", 1000, "proton", -1)'
```
