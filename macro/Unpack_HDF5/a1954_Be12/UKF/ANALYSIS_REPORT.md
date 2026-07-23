# a1954 ¹²Be — full analysis report (UKF + GENFIT)

Autonomous end-to-end analysis, 2026-07-21/22. 54 runs (0142–0215, usable), HDBSCAN reco,
IC + PID gated, two-body kinematics, both fitters. All plots in `pp/plots/`, `pid/plots/`.

## Pipeline
raw HDF5 → **HDBSCAN reco** (multipeak PSA + AtDirDeDxCleaner + AtTrackFinderHDBSCAN) →
**UKF** (`AtFitterUKF`) and **GENFIT** (`AtGenfitter`, Yassid/GenFit fork 2584bfe) fits →
`dump_kine_Be12.C` (KE/θ + PID + IC per track, event-ID matched to FRIB) →
`compute_ex_Be12.py` (IC(12Be 625–750) + PID-gate + χ² cut → Ex).

## Particle ID
- PID plane = `AtSpyralPID` (√dE/dx vs Bρ, pulse **maximum** = repo default). HDBSCAN gives
  the tightest band (beats Triplclust; ≈ Spyral 1.0.0). arclen>200 sharpens it.
- Beam = fragment cocktail; **¹²Be = IC-amplitude peak 625–750** (not the dominant ~1900).
- Gates: `pid/proton_band_Be12.json` (auto), `pid/deuteron_12Be.json` (user-drawn).

## Beam energy — CALIBRATED FROM DATA
The (p,p′) elastic peak must sit at Ex=0. Scanning Ebeam → **Ebeam ≈ 155 MeV (12.9 MeV/u)**
(UKF). Self-consistency check PASSES: the ¹²Be 2⁺₁ lands at 2.10 MeV (its known value).
⚠ Confirm against the real a1954 beam tune. (`calibrate_ebeam.py`, `ebeam_calib.png`.)

## ¹²Be(p,p′) — 13,786 gated protons (UKF)
- Sharp **elastic peak at Ex=0**, FWHM **0.47 MeV**.
- **2⁺₁ inelastic peak at 2.10 MeV** (~5% of elastic) — the known first excited state.
- Resolution: flat vs vertex-z (good calib); degrades for far vertices; vz<600 → FWHM 0.43
  but costs the weak 2⁺ stats → use full stats.

## ¹²Be(p,d)¹¹Be — 1,723 gated deuterons (UKF)
- ¹¹Be g.s. region + excited-state structure at ~2.5, ~4.5, ~7.5 MeV (Ex(11Be)).
- ⚠ **~0.5 MeV systematic offset** on the g.s. (both fitters, Ebeam-insensitive) → a deuteron
  KE/energy-loss systematic to calibrate out. Relative peak spacing is meaningful.

## UKF vs GENFIT (key comparison result)
| | UKF | GENFIT |
|---|---|---|
| median χ²/ndf | 0.10 (loose meas errors) | 6.31 (proper errors) |
| (p,p′) elastic FWHM | 0.47 MeV | 0.35 MeV (sharper) |
| calibrated Ebeam | **155 MeV** | 93–98 MeV (spurious) |
| 2⁺ resolved at 2.10? | **YES** | no (smeared) |

**GENFIT underestimates backward-angle KE** (UKF/GENFIT KE ratio 1.1 at θ 40–60° → **1.6 at
60–80°**), giving a wrong beam-energy calibration and smearing the inelastic states. **UKF is
the reliable fitter** for the full spectrum here (2⁺ at the correct energy). GENFIT's backward-
angle KE bias is worth investigating (energy-loss / backward-seed handling; a1975 hit similar).

## Key files
- Kinematics CSVs (persisted): `pid/kine_data/kine_{ukf,genfit}_{pp,pd}.csv`
- Macros: `pid/dump_kine_Be12.C`, `pid/draw_gate_Be12.C`, `pid/apply_gate_Be12.C`
- Analysis: `~/attpc_spyral_1.0.0/{compute_ex,calibrate_ebeam,compare_fitters_ex}.py`
- Plots: `pp/plots/ex_FINAL_summary.png`, `ex_UKF_vs_GENFIT.png`, `ebeam_calib.png`;
  `pid/plots/pid_Be12_full.png`

## Open items for you
1. Confirm the real a1954 beam energy (data says ~155 MeV / 12.9 MeV/u).
2. Refine gates (draw a proton gate; tune the deuteron gate top edge).
3. Deuteron KE offset (~0.5 MeV on 11Be g.s.) — energy-loss calibration.
4. GENFIT backward-angle KE bias — investigate.
5. Angular distributions / cross sections (DWBA) once gates finalized.
