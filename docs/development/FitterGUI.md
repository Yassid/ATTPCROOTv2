# Interactive UKF Fitter GUI

## Goal

Build a standalone interactive tool for visualizing and tuning the UKF track
fitter. Inspired by [GenFit's EventDisplay](https://github.com/GenFit/GenFit),
which provides a self-contained singleton with live refitting, error
visualization, and fitter parameter controls — all independent of the
experiment framework.

## GenFit EventDisplay Reference

Key design features to adopt:

- **Singleton pattern**: `EventDisplay::getInstance()` — one global instance
- **Standalone**: only needs `TEveManager` + geometry, no FairRoot dependency
- **Event storage**: vector of track collections, navigate with `next/prev/goto`
- **Live refit**: user selects fitter type and parameters, clicks redraw, track
  is refitted and redisplayed instantly
- **Display options as flags**: `"THDPME"` string toggles tracks, hits, detectors,
  planes, markers, errors
- **Self-contained GUI**: builds its own `TGCheckButton`, `TGNumberEntry`,
  `TGButtonGroup` panels via `makeGui()`
- **Error visualization**: `TEveGeoShape` ellipses/ellipsoids for covariance,
  `TEveTriangleSet` for error bands between measurement planes
- **Track rendering**: `TEveStraightLineSet` segments between planes, color-coded
  by charge sign

## Design

### `AtUKFDisplay` — Standalone Singleton

```cpp
class AtUKFDisplay : public TNamed {
public:
   static AtUKFDisplay *getInstance();

   // Event management
   void addEvent(const AtPatternEvent &patEvt, const AtTrackingEvent *fittedEvt = nullptr);
   void loadFile(const char *digiFile, const char *fittedFile = nullptr);
   void next(int step = 1);
   void prev(int step = 1);
   void gotoEvent(int id);

   // Fitting
   void fitCurrentTrack();  // Run UKF with current parameters
   void setParticle(int Z, int A, double mass_MeV, double charge_C);
   void setBField(double Bx, double By, double Bz);

   // Display options
   void setOptions(const char *opts); // "CHFSREM"
   // C=clusters, H=hits, F=filtered, S=smoothed, R=residuals,
   // E=errors, M=MC truth

private:
   void makeGui();
   void drawEvent(int id);
   void drawClusters(const AtTrack &track);
   void drawFittedTrack(const AtFittedTrack &fitted);
   void drawResiduals();
   void drawMomentumProfile();
   void drawCovariance();
};
```

### GUI Layout

```
┌─────────────────────────────────────────────────────────────────┐
│                    TEveManager 3D View                          │
│  ● blue = measured clusters                                    │
│  ● red  = smoothed track                                       │
│  ● green = filtered track                                      │
│  ○ covariance ellipses at each measurement plane                │
├──────────────────────┬──────────────────────────────────────────┤
│  UKF Controls        │  Diagnostics Canvas                     │
│  ┌────────────────┐  │  ┌──────────────┬───────────────┐       │
│  │ Event: [__] ←→ │  │  │ p vs index   │ residual vs i │       │
│  │                │  │  ├──────────────┼───────────────┤       │
│  │ Particle: [p ] │  │  │ theta/phi    │ chi2/ndf      │       │
│  │ B field:  [2.85│  │  └──────────────┴───────────────┘       │
│  │ Alpha:    [1e-3│  │                                          │
│  │ Meas σ:   [2.0]│  │  Fit Result:                             │
│  │ Min clust: [10]│  │    p = 45.2 ± 1.3 MeV/c                 │
│  │ ☐ Straggling   │  │    θ = 108.3 ± 0.8 deg                  │
│  │ ☐ Per-cl cov   │  │    KE = 1.09 MeV                        │
│  │ E-loss: [1.0]  │  │    chi2/ndf = 1.23                      │
│  │ ZPadPlane:[1000│  │    Convergence: YES                      │
│  │                │  │                                          │
│  │ [Fit Track]    │  │                                          │
│  │ ☐ Auto-fit     │  │                                          │
│  └────────────────┘  │                                          │
└──────────────────────┴──────────────────────────────────────────┘
```

### Display Options

| Flag | What it draws | TEve class |
|------|--------------|------------|
| C | Measured clusters | TEvePointSet (blue) |
| H | Raw hits (before clustering) | TEvePointSet (gray, small) |
| F | Filtered states | TEvePointSet (green) |
| S | Smoothed states | TEveStraightLineSet (red) |
| R | Residual vectors (smoothed→measured) | TEveLine (yellow) |
| E | Covariance ellipses | TEveGeoShape (semi-transparent) |
| M | MC truth track (if available) | TEveLine (white, dashed) |

## Implementation Plan

### Phase 1: Standalone Visualization Macro

ROOT macro that reads `output_digi.root` + `output_ukf_only.root`:
- Draws clusters and fitted track in TCanvas (XY, XZ, YZ projections)
- Momentum profile and residual plots
- Quick prototype, no interactive fitting

### Phase 2: `AtUKFDisplay` Compiled Class

- New directory: `AtEventDisplay/AtUKFDisplay/`
- Files: `AtUKFDisplay.h`, `AtUKFDisplay.cxx`
- Standalone singleton, no AtViewerManager dependency
- Self-contained GUI via `makeGui()` with ROOT TG widgets
- Live refitting: change parameters → click Fit → instant redraw
- 3D visualization with TEvePointSet, TEveStraightLineSet, TEveGeoShape
- 2D diagnostics in embedded TCanvas (momentum, residuals, chi2)

### Phase 3: Advanced Features

- Load MC truth from attpcsim.root for comparison overlay
- Error band visualization between measurement planes
- Multiple fitter comparison (side-by-side UKF vs GenFit)
- Export fit results to text/ROOT file
- Batch mode: fit all events, save summary histograms

## File References

| File | Role |
|------|------|
| GenFit `EventDisplay.h/cc` | Reference design (singleton, makeGui, drawEvent) |
| `AtReconstruction/AtFitter/AtFitterUKF.h` | Fitter to integrate |
| `AtData/AtPatternEvent.h` | Input: tracks with clusters |
| `AtData/AtTrackingEvent.h` | Output: fitted tracks |
| `AtData/AtFittedTrack.h` | Fitted result with kinematics |
| `macro/Simulation/ATTPC/16C_pp/run_ukf_only.C` | UKF via FairTask reference |

Sources:
- [GenFit GitHub](https://github.com/GenFit/GenFit)
- [GENFIT paper (arXiv:1410.3698)](https://arxiv.org/abs/1410.3698)
