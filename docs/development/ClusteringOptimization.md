# Clustering Optimization for UKF Track Fitting

## Problem

The UKF fitter converges on 83% of events, but the remaining 17% fail.
The clustering step (`ClusterizeSmooth3D`) is suspected as a contributing factor:

1. **Cluster positions may not accurately represent the track** — the smoothing
   pass inserts midpoints and re-clusters with half radius, potentially moving
   cluster centers away from the actual track
2. **Cluster ordering is fragile** — the nearest-neighbor walk works well but
   depends on a good starting point (vertex end identification)
3. **Cluster covariance may not reflect actual resolution** — hardcoded pad
   resolution and diffusion parameters affect the UKF's trust balance
4. **Brho momentum seed uses the RANSAC circle fit radius** — which is computed
   in the XY projection of the digi frame, not the full 3D helix

## Current Clustering Pipeline

```
Raw hits (AtHit, ~280 per track)
  ↓
ClusterizeSmooth3D(radius=15, distance=30.5)
  ├── Pass 1: Walk hits by distance, collect within radius, charge-weight XY
  ├── Pass 2: Insert midpoints, re-cluster with radius/2
  └── Output: ~65 AtHitCluster per track
  ↓
AtFitterUKF (nearest-neighbor walk + min spacing filter)
  └── ~50 ordered clusters for UKF
```

## Diagnostics Plan

### Step 1: Visualize clustering quality in the GUI

Add to `AtUKFDisplay`:
- Show **raw hits** (gray, small) alongside **clusters** (blue, large)
- Show **cluster covariance ellipses** or error bars
- Show **inter-cluster distances** as a diagnostic plot
- Color clusters by their **charge** to see the Bragg peak
- Mark the **seed cluster** (vertex) in a different color

### Step 2: Compare clustering methods

Test alternative clustering approaches:
- **Fixed-spacing**: place clusters at uniform arc-length intervals
- **Time-bucket clustering**: group hits by time bucket (natural Z slicing)
- **Charge-weighted centroid only** (no smoothing pass)
- **Different radius/distance** parameters (we already scanned these)

### Step 3: Improve Brho determination

The RANSAC circle fit gives `GeoRadius` in the XY digi projection:
- This underestimates the true helix radius for tilted tracks
- The `GeoTheta` correction `Brho = B*R/sin(theta)` helps but theta is also
  from the digi frame projection
- Alternative: fit a helix to the 3D cluster positions directly

### Step 4: Evaluate with GUI

For each change, use the GUI to:
1. Fit with MC truth seed → see clustering effect in isolation
2. Fit with Brho seed → see combined seed + clustering effect
3. Compare residuals and chi2 between methods

## File References

| File | Role |
|------|------|
| `AtTools/AtTrackTransformer.cxx` | ClusterizeSmooth3D implementation |
| `AtReconstruction/AtFitter/AtFitterUKF.cxx` | Nearest-neighbor ordering + filtering |
| `AtReconstruction/AtFitter/AtUKFDisplay.cxx` | GUI for interactive testing |
| `AtReconstruction/AtFitter/AtClusteringScanTest.cxx` | Compiled clustering parameter scan |
| `AtReconstruction/AtPatternRecognition/AtPRA.h` | SetClusterRadius/Distance |
| `AtReconstruction/AtPatternRecognition/AtTrackFinderHC.cxx` | Hardcoded radius=15, distance=30.5 |
