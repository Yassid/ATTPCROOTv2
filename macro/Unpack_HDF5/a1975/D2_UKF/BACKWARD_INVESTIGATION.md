# Backward-track efficiency — investigation & roadmap (16C(d,p)17C, a1975 D2)

Autonomous follow-up to the overnight production. The (d,p) protons go **backward**
(θ_lab > 90°), so backward efficiency *is* the (d,p) efficiency. Two threads flagged:
(1) how good is the clustering for backward tracks, (2) fix the UKF conventions.

State going in: the genfit `SetBackwardSeedFix` already recovers backward protons
(reversal 52%→16%, backward candidate fraction 5.6%→17.8%), but the **good-fit
fraction for backward tracks is still only ~28% vs ~71% forward** — so something
beyond the seed direction still limits backward efficiency. These notes localize it.

---

## Thread #1 — Clustering quality is markedly worse for backward tracks

Tool: `cluster_quality_dp.C` (reads PRA tracks from reco_d2; each AtTrack carries its
raw hits AND its clusters). Run on 0016–0019 (~104k PRA tracks):

| | N tracks | <clusters> | <hits/cluster> | <mean gap mm> | <pathlen/chord> |
|---|---|---|---|---|---|
| **forward** (θ_geo<90) | 68,704 | **42.5** | 3.6 | 15.4 | 3.20 |
| **backward** (θ_geo≥90) | 34,985 | **22.5** | 2.4 | **29.3** | 4.55 |

Backward tracks have **half the clusters**, **~2× the cluster gap**, **fewer hits per
cluster**, and are **more curled** (pathlen/chord 4.55 vs 3.20). Fewer clusters = fewer
fit constraints with a longer lever arm → this is the **second limiter** on backward
fit efficiency, on top of the (now-fixed) seed direction. Plot:
`plots/cluster_quality_dp.png`.

### The clustering deficit is REAL physics, not a binning artifact (corrected)
Initial suspicion was a sharp *algorithmic* step at 90°. **That was checked and ruled
out** (`cluster_step_pra_vs_3d.png`): binning clusters/track by the PRA GeoTheta and by
the independent 3D-geom theta gives the SAME profile (curves overlap) — a *smooth*
minimum at ~90° (step right at 90° is only 1.05×, not a discontinuity). Tracks near 90°
are perpendicular to the beam → shortest extent in the detector → genuinely fewest
clusters. So the backward cluster deficit is **real physics that the fitter must cope
with** (fewer constraints), NOT a clustering bug. Possible mitigation: finer/adaptive
clustering for short backward tracks (the UKF already has `fTargetClusters` re-clustering
— worth A/B testing on backward good-fit fraction).

**Next steps (thread #1):**
- Inspect the clustering + PRA ordering code for a θ=90° (or "going-away-from-pad-plane"
  vs "toward") branch: `AtSampleConsensus`/`AtPRA` ordering, `ClusterizeSmooth3D`, the
  ArcWalk extender (`project_attpcroot_arcwalk`), and `OrderClustersAlongTrack`.
- Re-run reco on 1–2 runs with **persistRaw=true** (intermediate branches) to measure
  the PSA-hit → cluster → PRA chain per angle (current reco_d2 is pattern-only, so the
  upstream hit→cluster loss is not yet visible; the deficit above is at the cluster
  level only).
- Test whether forcing more/finer clusters on backward tracks (smaller cluster radius,
  or the adaptive re-clustering the UKF already has, `fTargetClusters`) raises the
  backward good-fit fraction.

---

## Thread #2 — UKF loses backward tracks → DIAGNOSED & FIXED (existing flags, no code change)

Comparison (10 starter runs, same gate): **genfit 17.8% backward vs UKF 0.5% backward**
(`genfit_vs_ukf_dp.png`).

**The `:233` "convention bug" hypothesis was WRONG (checked & discarded).** Line 233's
`thetaLab = 180−GeoTheta` only feeds a *symmetric* beam-like cut, so the convention flip
doesn't change which tracks it drops. And empirically the UKF output theta is in the
SAME frame as genfit (95.7% agreement, not flipped — `SetZPadPlane` already converts to
lab). So it is neither a `:233` nor an output-frame bug.

**Real cause (empirical, of genfit's validated backward protons on run_0016):**
the UKF **drops 76% (no fit)**, fits 15% **as forward** (flip), 9% bad-chi2, 0.2% correct.
Two independent levers, both EXISTING AtFitterUKF knobs (default off / 10 → so (p,p)/(p,d)
untouched, NO framework code changed):
- **`SetUseClusterDirSeed(true)`** — seeds the direction from the cluster geometry
  (vertex = closest-to-axis cluster, dir = PRA-circle tangent disambiguated by the chord)
  instead of the half-sphere-ambiguous GeoTheta. Fixes the FLIP (15%→5.5%). Its own code
  comment already names this exact failure ("upward tracks see GeoTheta=115° not 65°").
- **`SetMinClusters(4)`** (was 10) — backward (d,p) tracks are cluster-poor (§Thread #1),
  so min=10 DROPPED 76% of them. min=4 (matching genfit) recovers them: no-fit 77%→24%.

**★ VALIDATED, run_0016:** clusterDirSeed + minClusters=4 →
- good-fit of genfit-validated backward protons: **0.2% → 19.2%**
- backward candidate fraction through ex_dp: **0.5% → 15.2%** (genfit 17.8%).

**Landed (no framework change):** `fitUKF_a1975_deuterium.C` now defaults
`clusterDirSeed=true` + `minClusters=4`. Residual (~24% no-fit + ~26% flip on the
shortest tracks) is the real clustering deficit (§Thread #1) — short backward tracks are
genuinely hard; adaptive clustering is the next lever there.

---

## GeoTheta was suspected as a common root — INVESTIGATED & CLEARED (93% reliable)

`AtPRA.cxx:289–315` computes the track's GeoTheta from the **2D XY-projection** of the
RANSAC direction, with the forward/backward sense set by a discrete sign:

```cpp
auto dir2D = XYVector(dir.X(), dir.Y()).Unit();
int sign = (dir2D.X()*dir2D.Y() < 0) ? -1 : 1;
angle = acos(sign * fabs(dir2D.Y())) * RadToDeg;   // sign=+1 -> [0,90], sign=-1 -> [90,180]
track.SetGeoTheta(angle * pi/180);
```

So forward-vs-backward (GeoTheta ≷ 90°) is decided ENTIRELY by `sign(dirX·dirY)` of the
**2D projection** — the 3D drift (z) direction is never used. Consequences:
- The forward/backward assignment is fragile/projection-dependent — exactly the kind of
  thing that makes backward tracks mis-seeded and that both genfit and UKF inherit (each
  then interprets the inherited GeoTheta with its own convention, §Thread #2).
- (It does NOT cause the clustering-vs-angle profile — that was checked and is robust to
  the angle definition, §Thread #1. GeoTheta's damage is to the fwd/bwd SEEDING the
  fitters inherit, not to the clustering metrics.)

**This is the highest-leverage thing to fix:** give PRA a robust 3D GeoTheta (use the
drift-z direction / the z-ordered cluster sequence to set the forward/backward sense),
and both the clustering-metric discontinuity and the cross-fitter backward seeding
become well-defined.

### ⚠ RETRACTED: PRA GeoTheta is actually GOOD (93%) — do NOT fix it
A follow-up check against the **validated truth** (genfit *fitted* theta on 19,270
good-fit run_0016 tracks) measured forward/backward agreement:
- **PRA GeoTheta: 93.0% correct** | my 3D-geom reference: 81.5%.

So the PRA GeoTheta is the BETTER estimator — it is **not** the bottleneck. The 32.7%
"disagreement" below was my crude 2-point 3D-geom reference being WORSE, not GeoTheta
being wrong. The current GeoTheta = `acos(slope/√(1+slope²))` with `slope` from the
**outlier-robust RANSAC fit of z-vs-arclength** (`AtPRA.cxx:289-315`) is a sound 3D
pitch; the genfit seed-fix already relies on it successfully (17.8% backward, physics-
consistent). **Recommendation reversed: leave PRA GeoTheta alone.** The real fixable
backward issue is the per-fitter seeding/conventions (§Thread #2, esp. the UKF), plus
the real clustering deficit (§Thread #1). The naive geotheta_check below is kept only
as the cautionary example that motivated this verification.

### (superseded) geotheta_check.C, 102,704 tracks, runs 0016–0019
Compared PRA GeoTheta to an INDEPENDENT 3D geometric theta (clusters in lab frame
z_lab=ZPadPlane−z_digi, vertex = cluster nearest the beam axis):

- **32.7% of tracks DISAGREE on forward-vs-backward** (>90°) between the two methods.
- The confusion is concentrated as a bright cross at **θ ≈ 80–100°** in the 2D plot
  (`plots/geotheta_check.png`) — exactly where the 2D sign is degenerate (|dirY|→0)…
  …and **exactly where the (d,p) protons live**. So the backward protons sit in the
  single worst angular region for the current GeoTheta. This is quantitative proof that
  the forward/backward determination is NOT robust, and it is upstream of both fitters.
- Caveat: the 3D-geom reference is not ground truth either (vertex-end ambiguity for
  forward near-axis tracks), so 32.7% is the *disagreement* rate, not a pure GeoTheta
  error rate — but a 1/3 disagreement concentrated at the physics-relevant angle is
  decisive that the classification needs a robust 3D treatment.

**Conclusion (revised after verification above): PRA GeoTheta is 93% reliable — do NOT
change it** (risks the (p,p)/(p,d) regression for no gain). The real backward fixes are
downstream: the per-fitter seed/conventions (§Thread #2 — the UKF `:233` handedness is
the concrete one), and the real clustering deficit (§Thread #1). Priority order:
**(1) UKF convention fix, (2) adaptive clustering for short backward tracks.**

---

## Quick-start for resuming

```bash
cd .../a1975/D2_UKF && source .../build/config.sh
root -l -b -q 'cluster_quality_dp.C("run_0016,run_0017,run_0018,run_0019")'  # thread #1
# thread #2: edit AtFitterUKF.cxx:233 convention -> rebuild AtReconstruction ->
root -l -b -q 'fitUKF_a1975_deuterium.C("run_0016",-1,"proton",-1)'          # re-fit UKF
root -l -b -q 'genfit_vs_ukf_dp.C("proton_kin_gf10.root","proton_kin_ukf10.root")'  # backward % check
```
