# a1975 — the procedure moved

This document, the analysis macros and the data caches now live together in one repository:

    ~/a1975_analysis

with `ANALYSIS_PROCEDURE.md` at its root. Read it there; this copy is not maintained.

## What stayed here

Only the parts that genuinely need ATTPCROOT — unpacking, PSA, pattern recognition and the
genfit/UKF fitters under `UKF/` and `D2_UKF/`. Everything downstream of the kinematics caches is
plain ROOT and was moved out, so the analysis can be run without building this framework.

The caching macros themselves are mirrored into each channel's `production/` directory over there,
next to the analysis they feed.
