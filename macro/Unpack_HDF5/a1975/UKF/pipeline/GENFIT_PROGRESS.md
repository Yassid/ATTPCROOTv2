# AtGenfitter (a1975) — progress & TODO

Clean genfit track fitter `EventFit::AtGenfitter` + pipeline for a1975 (¹⁶C+p).
Status tracker for the validation effort.

## Done

- **Clean fitter** `AtReconstruction/AtFitter/AtGenfitter.{h,cxx}`: genfit core only
  (KalmanFitterRefTrack + RKTrackRep + AtSpacepointMeasurement), modern `EventFit::AtFitter`
  base → plugs into `AtFitterTask`. No merging/reclustering/legacy machinery.
- **Cluster ordering by drift z** — monotonic along any helix, so curlers/spirals sequence
  correctly (fixed ev1684 spiral).
- **Deterministic vertex** = highest z_digi (lowest z_lab), the PRA/UKF convention. No more
  near-axis forward-track flips.
- **Realistic measurement covariance** — was a hard-coded 10 mm placeholder (seed-dominated,
  χ²/ndf≈0.01); now tunable via `SetMeasSigma` (default 4 mm), data-driven, χ² meaningful.
- **Upstream PID gate** `SetPIDGate(json)` — computes Spyral PID itself, fits only gated species.
  Uses the same `*_band.json` gates as the validated UKF chain.
- **θ acceptance window** `SetThetaWindow(10,170)` — drops unphysical near-beam/backward tracks.
- **Post-fit quality flag** `AtFitTrackMetadata::fGoodFit` — keeps everything inspectable, flags junk.
- **Fitted trajectory stored** (`SetSmoothedPositions`) for display/QA overlays.
- **Inspection tooling**: `inspectGenfitter.C`, `inspectFitEvent.C`, `seedVsFit.C`,
  `contactSheet.C`, `cycleEvents.C`, `cycleEventsGui.C` (button GUI).
- **Material effects** investigated: NO hang here (this genfit iteration-guards all RK loops →
  throws rather than spins; we catch). But mat-ON breaks protons (genfit internal Bethe-Bloch →
  degenerate fits) and only marginally helps deuterons → **default mat-OFF, σ=4 mm**. Per-track
  no-mat fallback added for thrown fits. Proton-in-H CATIMA table generated
  (`pp/make_catima_table_proton.C` → `resources/energy_loss/proton_H2_catima.txt`).
- **¹⁶C(p,p) elastic validation**: proton hypothesis (PDG 2212) + `proton_band.json`, mat-OFF →
  596 good fits / 759, 94% converged. KE-vs-θ_lab overlays the analytic elastic recoil line
  (`pp/overlay_pp_elastic.C`) → absolute calibration (B, angle, energy) confirmed.

- **¹⁶C excitation-energy spectrum vs UKF** (`pp/ex_genfit_vs_ukf.C`): same events, identical
  selection (IC + proton gate + χ² + θ window) and kinematics. run_0106/3000, IC-gated 1774:
  genfit 395 protons, elastic peak at Ex≈0, RMS 1.94 MeV; UKF 351 protons, RMS 2.17 MeV.
  → genfit gives **more statistics and a sharper elastic peak** than the UKF (resolution win,
  as expected). Small +0.3–0.65 MeV centroid offset in both (nominal beam E, no vertex E-loss).

## TODO

- [ ] **Continuity merging** — implement once we are happy with genfit. (Merge track segments
      that are continuous across PRA breaks before/within the fit.)
- [ ] **Material effects in the UKF** — add a model for **energy-loss straggling** and
      **multiple scattering** to the UKF propagation (the UKF currently has neither; genfit does
      via MaterialEffects). Needed for the UKF to match genfit's transport on longer/curved tracks.
- [ ] Tighten/justify the measurement covariance per detector resolution (currently effective 4 mm
      absorbing MS since genfit material effects are off by default).
- [ ] Vertex energy-loss correction for the beam (Ex currently uses nominal 192 MeV → broadens/shifts Ex).
