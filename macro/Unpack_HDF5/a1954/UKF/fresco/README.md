# FRESCO for a1954 `14C(p,p')`

DWBA reference curves for the a1954 elastic and inelastic angular distributions. The macros in
`../pp/` read `outputs/*.dat` from here by a path relative to themselves, so nothing needs
configuring.

## Regenerating

```bash
export PATH=$HOME/.local/bin:$PATH      # needs the `fresco` binary
python3 kd_params.py                    # writes inputs/*.nin
./run_fresco.sh p14C_el_161             # writes outputs/p14C_el_161{.out,_dsdo.dat}
```

`run_fresco.sh <stem>` runs one input in a scratch directory (so `fort.*` never lands here) and
splits the text log into per-state `_dsdo[_exN].dat` files of `mb/sr` vs degrees. `_dsdo.dat` is
the elastic channel (state 1) and `_dsdo_ex2.dat` the excited state.

Only the `.dat` files are tracked; the full `.out` logs are gitignored because `run_fresco.sh`
regenerates them from the `.nin`.

## What is calculated

Normal kinematics, proton on `14C`, at the lab energy giving the same `E_cm` as the beam:

| beam `14C` | `E_cm` | `E_lab(p)` |
|---|---|---|
| 161 MeV (the analysis value) | 10.804 MeV | **11.581 MeV** |
| 155 MeV (mean vertex energy) | 10.401 MeV | 11.149 MeV |

The predicted shape is insensitive to that difference, which is why both are kept — comparing
them is the cheapest test of beam-energy sensitivity.

Elastic prediction: steep forward fall, **deep minimum at ~62 deg**, **secondary maximum of
~19.7 mb/sr at ~93 deg**, second minimum near 140 deg. This reproduces the DWBA curve in
Ayyad et al., Eur. Phys. J. A **59**:294 (2023), Fig. 8 (lower), for the same data.

Inelastic inputs use a `type=11` deformation at multipole slot `L` with a reference deformation
length `beta_L R = 0.281 fm`, one excited state each: 6.094 (`1-`, L=1), 6.728 (`3-`, L=3),
7.012 (`2+`, L=2). Only the shape is used; the deformation length sets an arbitrary scale.

## Optical model

Koning-Delaroche 2003 global **proton** parameters, computed analytically in `kd_params.py`
for A=14, Z=6.

**Do not copy the a2091 `15C` potential.** Its real depth (`V_V = 46.89` at 13 MeV) comes from
the *neutron* asymmetry term `v1 = 59.30 - 21(N-Z)/A` and carries no Coulomb correction. For
protons the sign is `+21(N-Z)/A` plus `V_C`, which for p+`15C` gives ~55 MeV, not 46.9. Using the
neutron form here would bias the diffraction pattern that the whole comparison rests on.

At `E_lab(p) = 11.581 MeV`:

| | | | |
|---|---|---|---|
| `V_V` 54.807 | `r_V` 1.1357 | `a_V` 0.6757 | `W_V` 0.972 |
| `W_D` 8.880 | `r_D` 1.3042 | `a_D` 0.5260 | `V_SO` 5.510 |
| `r_SO` 0.9170 | `a_SO` 0.5900 | `r_C` 1.4778 | `E_F` -8.215 |

`W_SO` is -0.047 MeV and is omitted.
