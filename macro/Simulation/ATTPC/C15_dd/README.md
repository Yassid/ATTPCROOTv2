# 15C(d,d') simulation — acceptance and response for `levels_15C`

Why: every open item in `macro/Unpack_HDF5/C15p/levels_15C/LEVELS_15C.md` needs the same missing
ingredient — the **reconstruction** efficiency and the **response function**. The *kinematic*
acceptance is already computed exactly (`levels_15C/scripts/binacc.h`) and is ~1 inside a θ_CM
bin, so geometry is not what is missing.

```bash
source build/config.sh
cd macro/Simulation/ATTPC/C15_dd
# nEvents, thetaCM min/max, engine, B field (kG), output, Ex of the residual, seed, Ebeam
root -b -q 'C15_dd_sim.C(5000, 2., 70., "TGeant4", -28.5, "/home/yassid/Data/C15dd_sim/gs_sim.root", 0.0, 5001, 191.25)'
root -b -q truth.C        # does the analysis formula return the Ex that was generated?
root -b -q ebzscan.C      # what Ebeam(vz) slope minimises the truth-level spread?
```

## Validated so far

**The two-body missing-mass formula returns the generated Ex**, with a `+0.21` MeV offset at
Ex = 0 and `+0.29` at 6.68. That offset is the beam losing energy *before* it reacts
(`maxELoss = 10 MeV`, sampled uniformly, so ~5 MeV on average). It is **not** removable by
`Ebeam(vz)` — `<vz>` = 514 mm sits at the correction's 500 mm pivot, so the correction only takes
out the vz-*dependent* part. The data has the same effect, and the analysis absorbs it by tuning
Ebeam so the g.s. sits at zero: **191.25 is an effective value at the mean vertex, not the
injection energy**. Calibrate the simulation the same way; acceptance and response are both
insensitive to a 0.2 MeV scale offset.

**★ The simulation independently confirms `ebzSlope = 11 MeV/m`.** Scanning the slope and taking
the one that minimises the truth-level Ex spread:

```
slope:      0      2      4      6      8     10     12     14     16     18     20    MeV/m
Ex 0.00  0.183  0.167  0.153  0.143  0.136  0.133  0.135  0.142  0.152  0.166  0.182   rms (MeV)
Ex 6.68  0.192  0.168  0.147  0.129  0.117  0.112  0.115  0.126  0.142  0.163  0.186
```

Minimum at **10 MeV/m at both Ex values independently**, against the 11 MeV/m set empirically on
the data. The correction removes 27-42 % of the truth-level spread.

**The beam-energy spread alone contributes σ ≈ 0.11-0.13 MeV to the Ex resolution**, before any
detector effect. The measured g.s. is FWHM 0.91 (σ ≈ 0.39) for ATTPCROOT, so this is ~10 % of the
variance — real, but not the dominant term.

## Traps

- `FairParRootFileIo` **segfaults** rather than reporting a missing directory if the par file path
  does not exist. The par and geometry files are derived from `outFile` here for that reason.
- `ATTPC_D300torr_r280.root` (used by the a1975 `_r280` macros) is **not** in this repo's
  `geometry/`; use `ATTPC_D300torr_v2.root`, which pairs with the `_geomanager` file the C15d
  reconstruction uses.
- The residual 15C is transported as a **stable ion** even above Sn = 1.218 MeV. Correct for the
  deuteron kinematics, which is all the acceptance is measured on, but this simulation must not be
  used for the neutron or any coincidence.
