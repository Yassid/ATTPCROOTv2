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

### ★ Key lead: a sharp STEP at exactly θ_geo = 90°
All metrics step discontinuously at 90° (clusters 42→22, max gap ~50→100 mm,
hits/cluster drops). A *clean step at exactly 90°* is an **algorithmic boundary**, not
smooth physics (physics would be gradual). Strong suspicion: the hit-finding / cluster
ordering / PRA arc-walk has a forward-vs-backward branch keyed on 90° that under-serves
the backward branch.

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

## Thread #2 — UKF conventions lose backward tracks

Comparison (10 starter runs, same gate): **genfit 17.8% backward vs UKF 0.5% backward**
(`genfit_vs_ukf_dp.png`). The UKF (`AtFitterUKF`) mis-handles backward tracks the way
genfit did pre-fix. Root cause localized in `AtReconstruction/AtFitter/AtFitterUKF.cxx`:

- **`:233` — `double thetaLab = 180.0 - GeoTheta*180/π;`** This is the **simulation**
  z-frame convention (sim digitizes z_digi = ZPadPlane − z_MC, so theta is mirrored).
  EXPERIMENTAL data has the *opposite* z-handedness (the same reason the fit macro flips
  `bFieldSign=-1`). But this `180 − GeoTheta` is **hardcoded** and NOT flipped with the
  handedness, so for exp data the lab angle (and the forward/backward sense) is mirrored.
  This contradicts the genfit-validated mapping where **GeoTheta>90 *is* backward lab**.
- **`:117` vertex = `GetHitClusterArray()->front()`** and **`:128` seed momentum at
  `θ = GeoTheta`** directly. With the mirrored convention, the vertex-end and the seed
  direction are taken in the wrong frame for backward tracks → genfit-style reflection.

Contrast — what genfit does (the fix that works): positions in lab frame
`z_lab = ZPadPlane − z_digi`; `SetBackwardSeedFix` seeds GeoTheta>90 tracks from the
far (highest-z_lab) end. The UKF needs the analogous, handedness-consistent treatment.

**Proposed fix (NEEDS physics validation by Yassid — do not land blind):**
1. Make the GeoTheta→θ_lab convention at `:233` depend on the data handedness (same
   switch as `bFieldSign`), i.e. `θ_lab = GeoTheta` for experimental a1975, not
   `180 − GeoTheta`. Add a `SetExperimentalHandedness(bool)` / reuse the B-sign.
2. Re-derive the vertex-end + seed direction from that corrected θ_lab (mirror the
   genfit `SetBackwardSeedFix`: backward = corrected θ_lab > 90 → seed from far end).
3. **Validate with the objective metric**: backward candidate fraction on run_0016
   should jump from ~0.5% toward genfit's ~17–18%. `genfit_vs_ukf_dp.C` already
   measures it; `fitUKF_a1975_deuterium.C` + `ex_dp_a1975.C(...,"_UKF")` produce the
   UKF spectrum to compare.

This was left UNCHANGED on purpose — conventions/handedness are physics calls that
need your sign-off, and the objective test above makes validation a 1-run check.

---

## ★★ Likely COMMON ROOT of both threads: how GeoTheta is computed

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
- It discretely partitions tracks into the two sign branches, which is almost certainly
  the **sharp step at θ_geo=90°** seen in every clustering metric (§Thread #1) — the step
  is an artifact of binning by this projection-based angle, on top of the real backward
  cluster deficit.

**This is the highest-leverage thing to fix:** give PRA a robust 3D GeoTheta (use the
drift-z direction / the z-ordered cluster sequence to set the forward/backward sense),
and both the clustering-metric discontinuity and the cross-fitter backward seeding
become well-defined. Worth checking first whether `dir2D` even carries reliable sign
info for near-90° (|dirY|→0) tracks. Recommend starting here before the per-fitter fixes.

---

## Quick-start for resuming

```bash
cd .../a1975/D2_UKF && source .../build/config.sh
root -l -b -q 'cluster_quality_dp.C("run_0016,run_0017,run_0018,run_0019")'  # thread #1
# thread #2: edit AtFitterUKF.cxx:233 convention -> rebuild AtReconstruction ->
root -l -b -q 'fitUKF_a1975_deuterium.C("run_0016",-1,"proton",-1)'          # re-fit UKF
root -l -b -q 'genfit_vs_ukf_dp.C("proton_kin_gf10.root","proton_kin_ukf10.root")'  # backward % check
```
