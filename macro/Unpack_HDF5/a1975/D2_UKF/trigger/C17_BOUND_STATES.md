# 16C(d,p)17C bound states: how the angular distribution was resolved

Experiment a1975, D2 target. Adopted 2026-11-09. Supersedes the result adopted on 2026-09-24.

## Summary

The bound-state angular distribution could not be described by DWBA (chi2 84/3). The cause was the
trigger. Events are recorded only if the mesh signal fires, and whether it fires depends on where the
heavy residual goes. A fiducial cut that keeps only events where the residual always fires the mesh
removes the distortion, and the result is described by an l=2 shape.

**The three bound states are not separated in this measurement.** The g.s. (3/2+), 0.217 MeV (1/2+) and
0.331 MeV (5/2+) lie inside one peak, with a reconstructed width (sigma about 0.25-0.3 MeV) comparable to
their spacing. Every number below is the sum of the three. "Resolved" in the title refers to the shape
problem, not to the states.

| theta_cm (deg) | Y | Background | eps (genfit) | Spyral validity | dsigma/dOmega (mb/sr) |
|---|---|---|---|---|---|
| 25-30 | 21 | 0.89 | 0.998 | 1.00 | 11.00(2.56) |
| 30-35 | 67 | 0.00 | 0.976 | 1.00 | 6.92(85) |
| 35-40 | 72 | 2.67 | 1.006 | 1.00 | 3.70(46) |
| 40-45 | 23 | 8.02 | 0.873 | 0.32 | 1.98(74) |

Errors are statistical only. Values are in `dp_bound_angdist_ADOPTED.csv`.

Fit of one normalisation S to the FRESCO l=2 sum (state 1 + state 3, equal weight) over the four bins:

| Curve | S | chi2/ndf |
|---|---|---|
| l=2, adopted | 0.294(24) | 0.90/3 |
| l=2, without the environment factor | 0.247(20) | 0.90/3 |
| l=0 | 17.2(1.6) | 42.0/3 |

## What was wrong

The result adopted on 2026-09-24 (vz 200-400, theta_cm 20-40) rose to 32.5 deg while every DWBA curve
falls. In (vz, theta_cm) the data showed a diagonal loss band: 20-30 deg survived only at the front of the
chamber, 10-20 deg only at the far end, and 30-40 deg everywhere. The simulation had no such band.

These explanations were tested and excluded: Spyral fitting, reconstruction configuration, drift
attenuation, event environment (flat in angle), beam spot and offset, the vz shape of the data, and the
FRESCO coupling bug in the decks. The right field sign and the data Ex smear in the acceptance did not
fix it either.

## Cause: the mesh trigger

The mesh and IC waveforms are in the raw files (`frib/evt/evt<N>_1903`) and were never read by ATTPCROOT.
Measured on 47 runs and 8222 gated backward tracks:

- No track was recorded without a mesh crossing.
- The fraction of events where the mesh fires early (within 100 bins of the IC) rises with the residual
  radius r_res at the pad plane: 0.07 at 5-10 mm, 0.46 at 14-17 mm, 0.78 at 20-23 mm, 0.92 at 26-30 mm.
  This is the residual leaving the ~20 mm beam hole.

Interpretation: when the Z=6 residual leaves the hole near the pad plane the mesh always fires. Otherwise
the proton has to fire it alone, which needs a tall pulse, so only tracks near 90 deg or stopping protons
do. Bound-state protons that exit the gas with low dQ/dt are lost at vz > 400 mm unless near 90 deg.

The trigger logic itself is inferred from these measurements and has not been confirmed.

## Method

Input: genfit proton fits, r280 geometry, CATIMA, -2.85 T, `/mnt/f/a1975/sim_rebuild/trig/gf_data.txt`;
simulation `gf_sim.txt` and `keys.txt`.

- Selection: chi2 <= 13.5, theta_lab 89-161 deg, KE >= 1.2 MeV, vz 200-900 mm.
- Beam energy at the vertex: E(vz) = 183.871 - 0.011071 vz MeV. Simulation vertices are shifted by 19.1 mm.
- **Fiducial cut:** r_res >= 28 mm, with r_res = (1000 - vz) tan(theta_res), computed from the proton
  kinematics. The same cut is applied to generated, reconstructed and data events. No trigger model is used.
- Yield: Ex -1.5 to 0.729 MeV. Background: sideband Ex -4 to -1.5 MeV, scaled by width.
- Acceptance: reconstructed over generated in r_res >= 28, with the Ex window applied through a Gaussian
  smear of 0.238 MeV.
- Luminosity: 316.98 mb^-1 over 971.73 mm, times the window length and the generated fraction of each angle
  bin passing the cut.

Two corrections are included:

- **Spyral validity.** The data are PID-gated and so need a successful Spyral fit, but the acceptance was
  computed without that condition. At 40-45 deg the protons are flat to the pad plane and Spyral fails on
  about two thirds. The simulated success fraction (0.32) agrees with the data (0.33(8)). It is about 1 at
  25-40 deg, so it only changes the last bin, from 0.53 to 1.66(62) mb/sr before the environment factor.
- **Environment, 0.84.** Overlaying simulated events on random data costs 22 %, but the loss grows with
  event size. Re-weighted to the event sizes of real bound (d,p) events, survival is 0.84. This is a lower
  bound on survival, because the environment is slightly overestimated.

## Stability

The four-bin fit gives S 0.245-0.253 (no environment factor) and chi2 <= 1.8/3 at R 28 and R 30 for each of:
chi2 < 5, Ex top 1.0 MeV, smear 0.15 and 0.30 MeV, vz 200-700 mm. For vz 300-900 mm, S is 0.23(3).
Inside the fiducial region data and simulation Spyral efficiencies agree at 25-40 deg.

## Limits

- Four bins only. Below 25 deg cannot be reached by this cut. The 2026-09-24 points at 20-25 deg are not
  part of this result.
- The R cut was chosen after looking at the fit.
- The 25-30 deg bin has 21 counts and has not clearly levelled off: 5.94(1.01) mb/sr at R >= 26 against
  9.24(2.15) at R >= 28, before the environment factor. R 24 and R 26 give S about 0.207 with chi2 6-8/3.
- **The absolute scale is open.** The mesh efficiency is measured only among triggered events, so its
  plateau is assumed to be about 1.
- S is one number for an equal-weight sum of the state 1 and state 3 curves. It is not a spectroscopic
  factor for any single state. An l=0 + l=2 mixture was not fitted on these points, so the 0.217 MeV
  contribution is not constrained.
- PLB 811 (Pereira-Lopez et al.) is GANIL/TIARA at 17.2 MeV/u with a different reaction model. The ratio
  to its C2S (0.03, 0.64, 0.62) on our zero-range FRESCO curves is 0.69, 0.75, 0.70, 0.57. Part of that
  difference may be the model.
- The unbound states cannot be made trigger-safe this way: r_res >= 28 mm empties Ex 2-8 MeV.

## Open

- 16C(d,d) elastic scale check on the same data: deuteron fits running in `/mnt/f/a1975/gf_dd_r280`.
- Confirmation of the trigger logic (coincidence, gate, threshold), and the 22 % of triggers with no mesh
  pulse above 60 ADC.

## Reproduce

    cd macro/Unpack_HDF5/a1975/D2_UKF
    ~/Spyral/venv/bin/python -I trigger/adopt_fid_dp.py $PWD/dp_bound_angdist_ADOPTED

`adopt_fid_dp.py` runs `/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py` (same as `trigger/dsigma_fid.py`) and
reads `/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root`. The evidence for the Spyral term, the environment
factor and the stability tests comes from `fid_spy.py`, `census_sv2.py`, `env_rew.py`, `ev_size.C`,
`fid_spy_sys.py` and `unb_fid.py` in `trigger/`, which are not committed at the time of writing.
The 2026-09-24 result is kept as `dp_bound_angdist_ADOPTED_2026-09-24.csv`.
