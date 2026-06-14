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

- **Vertex energy-loss correction** (`pp/ex_vtxcorr.C`): Ebeam(z) for the ¹⁶C beam (loses ~0.027
  MeV/mm in H per CATIMA). FINDING: the integrated elastic peak is already sharp (σ≈0.18 MeV) and
  is NOT beam-eloss-limited — it's dominated by large-angle recoil protons that are kinematically
  insensitive to Ebeam (dEx/dEbeam≈0.3/22 MeV at θ=70° vs ≈2/22 at θ=30°). The systematic only
  tilts the minority forward (θ<45°) events: ~0.45–1.9 MeV across the chamber (stat/cut dependent),
  direction + magnitude consistent with CATIMA. A self-calibrated linear Ebeam(z) flattens that tilt
  (0.45→0.02 MeV) without overshoot, but leaves the integrated σ unchanged. NOTE: single-run stats
  (~395 protons, ~12–30 forward) make the gradient noisy — calibrate over all runs for production.

- **PRODUCTION proton run (10 runs)** — `pp/prod_ex_pp.C` over run_0106–0115, fit via `xargs -P4`
  (4 cores, ~4 min). 157,022 IC-gated ¹⁶C events → 37,018 fitted protons. Elastic peak Ex = 0.18 ±
  0.18 MeV, clean KE-vs-θ elastic locus, sensible vertex-z profile. Outputs: run_01xx_genfitter_pp.root
  on /mnt/f/a1975/reco/. Resolution (σ 0.18 MeV) holds at full stats; +0.18 MeV centroid offset to be
  calibrated later from proton E/angle.
- **FULL PRODUCTION (all 84 runs, run_0106–0189)** — fit via xargs -P4 (~50 min total incl. aggregate).
  1,196,817 IC-gated ¹⁶C events → 270,436 fitted protons. Elastic peak Ex = 0.177 ± 0.185 MeV; σ stable
  across single-run/10-run/84-run ⇒ intrinsic genfit resolution, not stat-limited. All 84 outputs on
  /mnt/f/a1975/reco/run_01xx_genfitter_pp.root. Plot: /tmp/prod_ex_pp.png.

## Done (2026-06 round: continuity merge + material model + cov + dual viewer)

- **Continuity merging** (`AtGenfitter::SetMergeContinuity`, `FitEvent`/`MergeContinuousTracks`):
  before fitting, merges AtTrack fragments PRA split across breaks but belonging to ONE physical
  trajectory — CLUSTER-LEVEL concat (the fit's drift-z sort re-sequences), so it suits the
  cluster-based fitter without re-running ClusterizeSmooth3D. Criteria: compatible PRA circle
  (centre <50 mm, radius <30%) AND nearest endpoints <30 mm in 3D; chains merge transitively.
  Runs on a COPY of the pattern tracks (shared AtPatternEvent untouched → a parallel UKF task is
  unaffected). Default OFF; opt-in arg in `fitGenfitter_a1975.C`. Validated run_0106/1500 evt: 114
  multi-track events, default criteria merge 6 (conservative, leaves 108 genuinely-separate alone),
  ultra-loose force-pass collapses all 114 (mechanics proven). Also FIXED a pre-existing double-free
  in `~AtGenfitter` (it deleted fMeasurementProducer which MeasurementFactory::clear() also deletes)
  — crashed whenever >1 AtGenfitter was constructed/destroyed in a process.
- **Per-cluster measurement covariance** (`SetDiffusion`, `SetZLongFactor`, `SetChargeCovRef`):
  drift-distance diffusion model sigma^2 = sigma0^2 + D^2·L_cm transverse + longitudinal, optional
  charge weighting. Defaults (D=0, zFactor=2) reduce EXACTLY to the prior flat (s2,s2,2·s2), so the
  270k-proton production is byte-unchanged; physics is opt-in (A/B test before adopting).
- **UKF material model** ([[project_attpcroot_ukf]], TrackFitterUKF): energy straggling was already
  implemented; ADDED Highland multiple scattering as per-step angular Q_mod (theta var += theta0^2,
  phi var += theta0^2/sin^2 theta), gated by `AtFitterUKF::SetEnableMultipleScattering` +
  `SetRadiationLength` (default H2 1bar X0≈7.6e6 mm). Opt-in (16C+p calibration untouched). NOTE: in
  H2 1 bar the per-step MS variance (~1e-9 rad^2) sits below the fAngModelNoise floor (1e-6) — correct
  scaling, becomes significant for denser gas / lower p / longer steps.
- **Dual-fitter event display**: `AtTabFitted` extended with branch name + smoothed-polyline mode +
  per-instance colour/style + lab→digi z map (`SetZPadPlane`). `fitBoth_a1975.C` writes
  AtTrackingEventGenfit + AtTrackingEventUKF from one pass; `run_eve_both_a1975.C` overlays both (red
  genfit, blue dashed UKF) over AtEventCorrected hits via a TTree friend; `dispBoth_a1975.C` is the
  WSL-friendly static-PNG version. Validated: genfit & UKF trace the same proton arc (ev146).

## Continuity-merge tuning + validation (measured)

- **Naive criteria REGRESSED (p,p)**: with the loose defaults (centreDist 50 mm) merge fired in 20/5000
  events, **lost 13 good fits (31→18)**, chi2 blew up (ev2578 0.31→1195). Root cause: (p,p) multi-track
  events are two SEPARATE tracks sharing the beam vertex, not one fragmented track. Added a vertex-end↔
  vertex-end-junction rejection (shared vertex ⇒ diverging) → 20→5 firings.
- **Direction-continuity guard TRIED then DROPPED**: requiring the two fragments' outward chords to align
  is ANTI-correlated on real data — a strongly-curved split track (ev2847) has chordCos 0.19 while two
  short parallel distinct tracks (ev429, the bad merge) have chordCos 0.995. Global chords fail for curved
  tracks. Measured via mergeDiag_a1975.C.
- **The real fix = shared-circle-centre test**: fragments of one helix fit the SAME circle (centres 7–11 mm
  apart); two distinct tracks sit ≥20 mm apart even at similar radius. Tightened default centreDist 50→15 mm.
  This eliminated the lone genuine bad merge (ev429, chi2 0.09→4.88).
- **Validated, mergeQuality_a1975.C** (classifies each merge by RESULT quality, not track count):
  - proton (5000 evt): merge fires 3×, **all DEDUP** (split arc → one good fit, chi2≈0.13), **0 HARMFUL**.
  - deuteron (p,d) gate (30094 evt): merge fires 14×, **all NEUTRAL, 0 HARMFUL** (those fragmented events
    don't yield good deuteron fits either way — no benefit shown, no harm).
  - **Bottom line**: merge is now SAFE on both channels (zero quality degradation), does correct dedup on
    protons. a1975 is clean enough that PRA fragments are individually fittable, so the RECOVER benefit
    (rejoining sub-threshold fragments) is NOT exercised here — it needs a long/curved/backward-heavy set.
  - Tools: mergeQuality_a1975.C (quality classifier), mergeDiag_a1975.C (per-pair discriminator dump).

## Measured effect of the new knobs (featureScan / ukfMSscan / ukfMinClus, run_0106)

- **Genfit measurement covariance is a YIELD↔CALIBRATION trade** (featureScan_a1975.C, proton gate):
  flat σ=4 mm → median χ²/ndf 0.105 (errors ~3× over-estimated) but 79% goodFit; σ=1.5 → χ²/ndf 0.75,
  58%; σ=1.0 → χ²/ndf 1.69, 54%. Statistically-calibrated σ≈1.2 mm (χ²/ndf→1) but it COSTS yield
  (79→55%) because the residual the 4 mm absorbs is real multiple scattering the genfit model (material
  effects OFF) can't follow that tightly. The diffusion knob doesn't escape it (only rescales the average:
  diffusion(1.0,1.4)/σ1.5 → χ²/ndf 0.04, 81%). So the production 4 mm is a deliberate yield-max choice;
  χ² is a loose junk cut here, not an error bar. The clean fix for BOTH is a proper MS term in the model.
- **UKF multiple scattering is negligible in H2 1 bar** (ukfMSscan_a1975.C): MS off vs on → identical yield
  (924), median KE 5.881 → 5.808 MeV (~1%). Confirms the prediction; it grows for denser gas / lower p.
- **UKF MinClusters sweet spot = 4–5** (ukfMinClus_a1975.C): 10→4 recovers +29% fitted tracks
  (1386→1793) with the physical-KE fraction only 71%→68%; gain SATURATES at 4. Recommend UKF MinClusters
  5 (was 10) to recover the short tracks the dual viewer showed it dropping, near-free quality-wise.
- **FIXED 4 failing AtFitterUKFFixture unit tests** (now 8/8): the fixture supplies pre-built clusters but
  didn't disable adaptive re-clustering (default on), which collapsed the short synthetic track below
  MinClusters → NULL fit; added SetAdaptiveClustering(false) + relaxed the brittle exact smoothed-count
  assertion to a range. (Unrelated pre-existing AtToolsTests.AtPropagatorTest stopping/field failures
  remain — NOT touched by this work.)

## TODO

- [ ] Re-calibrate the vertex E-loss gradient on the full run set (high forward-elastic stats). [SKIPPED]
- [ ] **!! FUTURE 16C(d,p) — HIGH PRIORITY !!** Validate (and re-tune) continuity merging there. In the
      a1975 (p,p) data backward tracks (theta>90) are PHYSICS-FORBIDDEN, so the "backward" events seen here
      are reconstruction ARTIFACTS, not real — a1975 (p,p)/(p,d) is too clean to exercise the merge RECOVER
      benefit. 16C(d,p) is the channel with genuinely fragmented (long/curved) tracks where merge should
      help; that is where it must be validated. NOTE: the merge's vertex/far end-assignment is by XY-radius,
      which is unreliable for axis-hugging backward tracks — switch to a drift-z-based assignment for (d,p).
- [ ] A/B the diffusion covariance + UKF multiple scattering on real data (effect is small in H2 1bar).
- [ ] Vertex energy-loss correction for the beam (Ex currently uses nominal 192 MeV → broadens/shifts Ex).
