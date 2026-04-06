  # Plan: Compare and Refactor Cluster-Covariance Construction

  ## Summary

  Implement a controlled comparison of three covariance modes:

  1. fixed_sigma: ignore per-cluster covariance and use the fitter’s existing fixed measurement sigma
  2. transformer_direct: keep the current AtTrackTransformer behavior, where cluster centroid and covariance are computed directly from pad/diffusion formulas during clusterization
  3. hit_cluster_online: compute per-hit variance from pad size and diffusion, then build each cluster through AtHitCluster::AddHit(...) so the cluster covariance comes from the existing online update logic

  The first deliverable is not a default switch. It is a reproducible comparison on both synthetic truth-driven tracks and the existing digitized 16C(p,p) workflow. If hit_cluster_online is clearly better,
  it should be left behind a configuration flag first, not made the default immediately.

  ## Implementation Changes

  ### 1. Add an explicit covariance-construction mode in clustering

  Add a small enum-style configuration on the clustering side with exactly these modes:

  - TransformerDirect
  - HitClusterOnline

  Do not include AtHitClusterFull or additional exploratory covariance estimators in this first pass.

  The clustering API should make two decisions explicit:

  - how cluster centroids are formed
  - how cluster covariance is formed

  For this comparison, keep the cluster membership and centroid policy as close as possible to the current behavior. The variable under study should be covariance construction, not a simultaneous change in
  clustering geometry.

  Concretely:

  - preserve the existing cluster membership logic in ClusterizeSmooth3D and ClusterizeByGroup
  - factor the repeated covariance-building code in AtTrackTransformer into a mode-dependent helper
  - for TransformerDirect, keep the current diagonal covariance calculation exactly as the baseline
  - for HitClusterOnline, build an AtHitCluster incrementally from the hits assigned to that cluster and let AddHit(...) produce the covariance

  Important implementation choice:

  - before adding raw hits into AtHitCluster, ensure every raw hit has a physically meaningful PositionVariance
  - if a hit arrives without variance set, derive it from the same detector model used by the transformer baseline:
      - x,y: pad transverse resolution plus transverse diffusion
      - z: longitudinal diffusion plus time-bucket resolution and pad-length contribution if that is part of the current model
  - keep this derivation in one shared helper so both covariance modes consume the same detector parameters

  ### 2. Unify the detector-model inputs used by both methods

  Create one shared helper for per-hit detector uncertainty evaluation, sourced from:

  - transverse diffusion coefficient
  - longitudinal diffusion coefficient
  - drift velocity
  - time bucket duration
  - pad transverse resolution
  - any chosen longitudinal pad-resolution model

  This helper should be used in both places:

  - current transformer-direct covariance construction
  - new hit-based covariance construction

  That ensures the comparison is about aggregation method, not inconsistent detector assumptions.

  Defaults to use:

  - if SetDiffusionParams(...) has been called, use those runtime values
  - otherwise use the existing transformer defaults
  - do not invent a second configuration source for this study

  ### 3. Expose the comparison cleanly in the fitter/scans

  Do not change the core fitter behavior beyond selecting whether it uses per-cluster covariance.

  Comparison matrix to implement:

  - fixed_sigma
      - clusterization unchanged
      - AtFitterUKF::SetUsePerClusterCov(false)
  - transformer_direct
      - cluster covariance built by current transformer method
      - AtFitterUKF::SetUsePerClusterCov(true)
  - hit_cluster_online
      - cluster covariance built via AtHitCluster
      - AtFitterUKF::SetUsePerClusterCov(true)

  Keep fixed fitter settings constant across the comparison unless a test is specifically about covariance sensitivity.

  Recommended constant settings for the comparison harness:

  - same clustering radius/distance within a given scan point
  - same momentum seed policy
  - same MeasurementSigma for the fixed_sigma baseline
  - same straggling setting within each comparison block
  - same ZPadPlane, field, and kinematic cuts

  ### 4. Build quantitative comparison into existing test infrastructure

  Extend the existing scan-style utilities rather than inventing a separate benchmark system.

  Use two layers:

  - synthetic truth-driven comparison
      - extend AtClusteringScanTest-style coverage
      - compare modes on momentum bias, momentum resolution, mean residual, convergence rate, and cluster-count stability
  - digitized end-to-end comparison
      - extend AtClusteringDigiScanTest-style coverage
      - compare modes on convergence rate, fraction of physically plausible fits, fitted KE/theta distributions, and average cluster count

  Add a compact summary table per mode and per clustering configuration.

  For synthetic tests, report at minimum:

  - convergence fraction
  - mean momentum bias
  - RMS momentum error
  - mean spatial residual
  - cluster count

  For digitized tests, report at minimum:

  - tracks tried
  - converged fits
  - “good fit” count using the existing kinematic window logic
  - average KE for accepted fits
  - average cluster count

  ### 5. Add one direct covariance-inspection test

  Add a small deterministic unit test around cluster construction itself, separate from the fitter.

  Use a hand-built micro-track with known per-hit variances and charges and verify:

  - TransformerDirect produces the expected diagonal covariance from the closed-form formula
  - HitClusterOnline produces a covariance consistent with AtHitCluster::AddHit(...)
  - both methods respond monotonically to increasing drift distance
  - both methods respond monotonically to increasing cluster size / summed charge in the expected direction

  This test is important because otherwise fitter-level differences will be hard to interpret.

  ## Test Plan

  ### Synthetic comparison

  Run controlled synthetic tracks with truth available and compare all three modes across multiple random seeds.

  Acceptance criteria:

  - no mode regression that causes widespread fit failure


  Run the extended digi scan on the current 16C(p,p) path with identical clustering settings across modes.

  - results table clearly separates covariance-method effects from clustering-parameter effects
  - hit_cluster_online should not degrade the rate of physically plausible fits
  - if it improves fit quality, the gain should hold across more than one clustering setting

  ### Low-level covariance validation

  Add direct tests for cluster covariance construction.

  Acceptance criteria:

  - no negative diagonal entries
  - covariance grows with diffusion / drift distance as expected
  - covariance shrinks with stronger averaging as expected
  - online AtHitCluster output is reproducible for deterministic hit input

  ## Public Interfaces / Config Additions

  Add one clustering-side selector for covariance construction with two supported values:

  - TransformerDirect
  - HitClusterOnline

  Do not expose AtHitClusterFull estimators in the public comparison interface yet.

  Use this selector in:

  - pattern-recognition clustering path
  - fitter-local adaptive re-clustering path, if that path remains active during comparison

  Add mode labels to comparison output so tables and logs are self-identifying.

  ## Assumptions And Defaults

  - The first study compares only the two realistic production covariance builders plus the fitter’s fixed-sigma baseline.
  - AtHitClusterFull is out of scope for the first round because it broadens the experiment without helping the immediate production decision.
  - The new method should remain behind a flag even if it looks better.
  - Detector-model inputs for hit variance should come from the same diffusion/pad model currently used by AtTrackTransformer, not from a newly invented model.
  - The comparison should reuse the existing synthetic and digitized scan infrastructure rather than introducing a new benchmarking framework.



› Part of this implementation should be writing tests to ensure the refactor does not change the current behavior at all. Do not fix any bugs or inconsistancies in the existing code at this point. That would
  poison the validity of the comparison to the current code.


• Proposed Plan


  # Plan: Add Covariance-Method Comparison Without Changing Current Behavior

  ## Summary

  The work should be split into two tracks that stay deliberately separate:

  1. a behavior-preserving refactor of the current covariance path
  2. a new alternate covariance path for comparison

  The baseline must remain byte-for-byte equivalent in behavior as far as practical. No bug fixes, no parameter cleanup, no inheritance cleanup, no “while we are here” improvements. If the current code has
  inconsistencies, they should remain in place during this work and be documented, not corrected. Otherwise the comparison stops being a comparison to the current code.

  The comparison should cover exactly three modes:

  - fixed_sigma: fitter ignores cluster covariance
  - transformer_direct: current behavior, preserved exactly
  - hit_cluster_online: new path where per-hit variance is set from the same detector model and clusters are accumulated through AtHitCluster::AddHit(...)

  ## Implementation Changes

  ### 1. Freeze the current behavior behind characterization tests

  Before introducing any new covariance mode, add tests that lock down the present behavior of the current clustering path.

  These tests should characterize the current implementation, not the intended implementation. They should assert what the code does today, including awkward details.

  Add characterization coverage for:

  - ClusterizeSmooth3D cluster positions and cluster count on deterministic input
  - ClusterizeSmooth3D cluster covariance values on deterministic input
  - ClusterizeByGroup cluster positions and covariance values on deterministic input
  - fitter behavior with SetUsePerClusterCov(false) versus true on a fixed cluster set
  - adaptive re-clustering path using the current local AtTrackTransformer defaults, without “repairing” parameter propagation

  The purpose is to make sure the refactor leaves the current mode unchanged.

  ### 2. Refactor the current transformer code without changing its outputs

  Refactor AtTrackTransformer only enough to isolate the current covariance calculation into a reusable helper.

  Rules for this refactor:

  - preserve the current formulas exactly
  - preserve the current diagonal-only covariance
  - preserve current use of z /= nHits while x,y are charge-weighted
  - preserve current smoothing and midpoint behavior
  - preserve current defaults and fallback parameter values
  - preserve current handling of unset diffusion parameters
  - do not change ordering, filtering, or cluster membership logic

  The helper should be treated as a frozen baseline implementation. Its job is to make the old path testable and reusable, not better.

  ### 3. Add a second covariance-construction mode alongside the baseline

  After the baseline is frozen by tests, add an explicit covariance-construction selector with two clustering-side modes:

  - TransformerDirect
  - HitClusterOnline

  TransformerDirect must call the preserved baseline helper and remain the default.

  HitClusterOnline should:

  - use the same cluster membership produced by the existing clustering logic
  - derive per-hit position variance from the same detector inputs used by the baseline method
  - create an AtHitCluster and feed member hits through AddHit(...)
  - use the resulting AtHitCluster::GetCovMatrix() as the cluster covariance

  Important constraint:

  - do not alter centroid or membership rules in this phase unless required to make AtHitCluster usable
  - if AtHitCluster naturally produces a different internal position estimate, keep that difference confined to the new mode only
  - the baseline path must not be rerouted through AtHitCluster

  ### 4. Keep the detector model identical between baseline and new mode

  Both covariance methods must consume the same detector-model inputs:

  - transverse diffusion coefficient
  - longitudinal diffusion coefficient
  - drift velocity
  - time bucket duration
  - pad transverse resolution
  - whatever longitudinal pad-resolution assumption the current transformer already uses

  Do not “improve” this model for the new method. The new method should differ only in how hit uncertainties are aggregated into a cluster covariance.

  ### 5. Extend the existing scan infrastructure to compare modes

  Use the existing scan-style tests/utilities to compare:

  - fixed_sigma
  - transformer_direct
  - hit_cluster_online

  The comparison should be run on:

  - synthetic truth-driven tracks
  - digitized 16C(p,p) data

  Keep all non-covariance settings fixed within a comparison block:

  - same clustering parameters
  - same fitter settings
  - same straggling choice
  - same momentum seed policy
  - same event selection

  The goal is to isolate covariance-method effects rather than re-open the entire reconstruction tuning problem.

  ## Test Plan

  ### A. Baseline characterization tests

  These are the most important addition.

  Add deterministic tests that lock down the current transformer_direct behavior:

  - exact cluster count for fixed toy inputs
  - exact or tight-tolerance cluster positions
  - exact or tight-tolerance covariance diagonals
  - zero off-diagonal terms
  - unchanged behavior for both ClusterizeSmooth3D and ClusterizeByGroup

  These tests should fail if the refactor accidentally changes present behavior.

  ### B. Baseline fitter-regression tests

  Add tests that confirm the refactor does not change current fitter behavior when using the baseline mode:

  - same convergence status
  - same number of fitted clusters
  - same fitted momentum/kinematics within tight tolerance
  - same chi2 or fit metadata where applicable

  This should be done for:


  Add targeted tests for hit_cluster_online:

  - cluster covariance is produced through AtHitCluster
  - covariance is symmetric and non-negative on the diagonal
  - covariance changes as expected when hit variance or drift distance changes
  - the new mode differs from the baseline only in covariance construction, not because baseline behavior was altered


  Report at minimum:

  - convergence fraction
  - cluster count

  ## Public Interfaces / Config

  Add one covariance-construction selector on the clustering side with two values:

  - TransformerDirect as the default
  - HitClusterOnline as the alternate mode

  Do not add AtHitClusterFull estimators in this phase.

  The fitter-side comparison still uses:

  - SetUsePerClusterCov(false) for fixed_sigma
  - SetUsePerClusterCov(true) for both cluster-covariance modes

  ## Assumptions And Guardrails

  - No bug fixes or cleanup to existing covariance behavior are allowed in this phase.
  - Existing inconsistencies must remain unchanged and, where relevant, be documented in test names or comments.
  - The baseline mode is the current code path, preserved exactly.
  - The new mode exists only as an alternate comparison path behind configuration.
  - Any future cleanup of the baseline should happen only after comparison results are collected and explicitly separated from this work.