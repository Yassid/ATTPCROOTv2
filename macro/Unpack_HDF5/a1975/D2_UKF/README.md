# a1975 16C(d,p)17C — D2-target UKF/GENFIT workspace

Reconstruction of the **deuterium (D2) target** runs of a1975 for the transfer
channel **16C(d,p)17C** (inverse kinematics: `16C + d -> p + 17C`). The detected
light ejectile is the **PROTON**, which goes largely **BACKWARD** in the lab — so
this pipeline is explicitly extended to backward tracks.

Distinct from the sibling `../UKF/` workspace (proton/H2 target, (p,p′) & (p,d)).
Source the env first: `source build/config.sh` from the repo root, then run macros
from this folder.

## Run set & parameters

- **D2 (d,p) physics runs = 0016–0103** (47 runs). Runs 0106–0189 are the H2/proton
  set (handled by `../UKF/`).
- Raw HDF5: `/mnt/f/a1975/h5/run_XXXX.h5` (drvfs, ~10–20 GB/run).
- Parameters: `parameters/ATTPC.a1975_deuterium.par` — identical to `ATTPC.a1954.par`
  except **DriftVelocity 1.15 cm/µs** (the D2 gas value; H2 uses 1.30).
- Reco/fit outputs: `/mnt/f/a1975/reco_d2/`.

## Pipeline

| Macro | Does | Output |
|-------|------|--------|
| `unpackReco_a1975_deuterium.C(run, nEv, persistRaw, outDir)` | unpack + PSA + space-charge + PRA | `<run>_reco.root` (AtPatternEvent) |
| `fitGenfitter_a1975_deuterium.C(run, nEv, ioDir, suffix, outDir, B, ...)` | **PROTON-hyp** genfit (θ-window 5–178° for backward protons) | `<run>_genfitter_p.root` (AtTrackingEvent + AtPIDEvent) |

### Batch drivers (resumable)

- `reco_batch.sh [Nparallel]` — unpack the 10-run starter set → `reco_d2/`.
  **Use 2-parallel**: 4-parallel thrashes drvfs on the 16 GB HDF5 files (confirmed —
  zero events in ~4 min at 4-way; 2 is the cap, page cache then accelerates).
- `fit_batch.sh [Nparallel] [bField]` — proton-hyp genfit over the reco files.
  CPU-bound (reads the smaller `_reco.root`) → **4-parallel is fine**.

## Backward-track notes

- The clean `EventFit::AtGenfitter` already accepts backward tracks
  (`SetThetaWindow`, default 10–170°; opened here to 5–178°) and seeds the vertex by
  z-ordering (closest-to-axis convention). 31% of PRA tracks in the test run are
  θ_geo ≥ 90°.
- **B-sign is unverified for D2**: experimental handedness `z_lab = ZPadPlane - z_digi`.
  Proton-target genfit (p,d) used B=+2.85; this macro defaults to −2.85. Verify on the
  fitted vertex (near beam axis) and physical proton KE; flip if needed.
- No proton PID gate yet — first pass fits all tracks; build the gate from the PID
  plane (`AtPIDEvent`: sqrt(dEdx) vs Brho) afterwards.

## Next

1. Validate the fit on one finished reco (vertex, KE, θ distribution; pick B sign).
2. Build the proton PID gate.
3. Two-body kinematics → **17C excitation-energy spectrum** (m1=16C, m2=d, m3=p, m4=17C).
