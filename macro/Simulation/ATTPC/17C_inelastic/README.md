# 17C(p,p') and 17C(d,d') — the M_n/M_p proposal simulation

The primary measurement of `C17p_FRIB_Proposal.pdf` (in this directory): the first determination of
M_n/M_p in 17C, from inelastic scattering of a 17C beam on hydrogen and on deuterium in the AT-TPC
inside SOLARIS. The (d,p) arm of the same proposal — the transfer channel the deuterium day gives
"for free" — is `../17C_dp/`.

**Results: [RESULTS.md](RESULTS.md).**

## The measurement, and why it is hard

The proposal wants dσ/dΩ for the **1/2⁺ at 217 keV** and the **5/2⁺ at 332 keV**. S_n(17C) = 733 keV,
so those two plus the 3/2⁺ ground state are the entire bound spectrum.

- The two states are **115 keV apart**. (The proposal text says 130 keV; 332 − 217 = 115.)
- The ground state is 217 keV below the 1/2⁺, so **all three sit inside one 300 keV resolution
  width**, and the elastic channel is far stronger than either inelastic.
- So the two peaks to be extracted sit on the flank of a much larger one, less than one resolution
  width away. The 14C(p,p') precedent the proposal cites (Ref. [24]) does not cover this: there the
  states of interest are at 6–7 MeV, about twenty resolution widths from elastic.

The proposal's own claim is correspondingly modest — not that the doublet is resolved, but that
"even with the 300 keV energy resolution achievable in the AT-TPC, we expect to be able to fit a
doublet in order to extract counts from the individual states". That is a claim about a
**constrained fit**, and no resolution number confirms or refutes it. `decompose_C17.C` tests it
directly.

## Running it

```bash
root -b -q 'inel_kinematics_C17.C'          # the gate: minutes, no campaign needed
./run_C17inel_campaign.sh -j 8 16000        # 12 samples, ~2 h on this box
root -b -q 'inel_summary_C17.C'             # resolution, floor, vertex correction, separation
root -b -q 'decompose_C17.C("...","pp_b285")'   # the headline: can the yields be extracted
root -b -q 'kine_lines_C17.C("...","pp_b285")'   # the loci, raw vs vertex-corrected
root -b -q 'export_gui_C17.C' && gui/build_gui.sh   # the interactive explorer, for figure captures
```

`accumulate_C17inel.sh <channel> <state> <field_T> <seed> [nEvents]` is one sample, and is the
thing to call by hand to re-run a single cell. Every stage is resumable on evidence that it
finished — a written file is not a finished job.

Products live in `/media/yassid/Seagate Hub/ATTPC/C17_inel/` (override with `INEL_ROOT`).

## The matrix

| | |
|---|---|
| channels | `pp` 17C(p,p') on 300 torr H₂ · `dd` 17C(d,d') on 300 torr D₂ — one beam day each |
| levels | `gs` 0 (3/2⁺) · `ex217` 0.217 (1/2⁺) · `ex332` 0.332 (5/2⁺) |
| fields | 2.85 T (what Ref. [24] ran at) · 4.00 T (SOLARIS's design field) |

**There is no pad-pitch axis.** This proposal runs the conventional AT-TPC pad plane, so the 2 mm
comparison the 14C campaign made is not on offer.

**There is a field axis, but the field turned out not to be the lever.** The axis was added on a
pre-campaign extrapolation that propagated 14C(p,p')'s single angle-averaged σ(KE) through this
channel's leverage and predicted a 4× gain at 4 T. That was wrong twice over — σ(KE) varies 20×
across the angular range rather than being one number, and the constant-E_beam term dominates until
it is corrected — so 4 T changes the as-analysed σ(E_x) by nothing and costs acceptance. It is worth
1.6× once the vertex correction is applied, which matters most on the deuteron day. See
[RESULTS.md](RESULTS.md); `inel_kinematics_C17.C` now states the *requirement* instead of predicting
the outcome.

What the gate does get right, and is worth running first, is the leverage: dEx/dKE is flat at
−0.533 (p) / −0.563 (d) across every angle, so the 115 keV level gap always appears as ~217 keV of
recoil energy and σ(E_x) is set by σ(KE) alone.

## What is held fixed

The 14C(p,p') reference campaign's configuration (`../14C_pp/highfield/`), which was debugged
against a1954 data: 1 m drift, 3 cm beam hole, `AtPSAMultiFit` + `AtDirDeDxCleaner` +
`AtTrackFinderHDBSCAN(mcs 20, ms 8, ε 10 mm, mover)`, genfit with material effects and native
CATIMA dE/dx, acceptance and resolution at χ²/ndf < 5 on `GetKinematicsXtr`. Only the beam, the
target and the field differ, so a difference against that campaign is the reaction and not a
setting. The elastic acceptance came out at 0.827 against the 14C campaign's 0.826, which is the
cheapest available check that nothing drifted.

Par files are inherited, not generated: `ATTPC.a1954_C14_hf_b285/b400.par` (H₂) and
`ATTPC.C17dp_D300torr_b285.par` / `ATTPC.a1954_C14dp_b400.par` (D₂). All four are already
gas- and field-matched, with CoefT at 4 T the Magboltz-scaled value anchored on the a1954-tuned
9e-4 at 2.85 T. Hence no `make_*_par.sh` here.

## Angular distributions

`fresco/` holds the DWBA calculations supplied with the proposal — `c17pp_217keV.out` and
`c17pp_332keV.out`, FRESCO, 17C + p at a lab energy of 136 MeV. That is not a different beam energy
from the proposal's 8.37 MeV/u: 142.29 MeV enters the chamber and ~14.5 MeV is lost crossing the
metre, so 136 is the mid-chamber mean.

Integrated, σ(217) = 8.58 mb and σ(332) = 3.60 mb. At 300 torr H₂ over 1 m with 940 pps for one
day that is **1379 and 579 counts at 4π** — the proposal's "~1000 scattered protons" is right for
the 1/2⁺ and about twice optimistic for the 5/2⁺, which is therefore the state that sets the
statistical error.

The **generator is flat in cos θ_cm** and the FRESCO shapes are applied as weights in the analysis,
never at generation. With a flat generator the per-bin acceptance is a clean ratio and the Ex
resolution is measured per angular slice, so both are independent of the true shape and a revised
DWBA calculation reweights the result instead of requiring a regeneration.

**θ_cm is the projectile angle**, so the light recoil comes out at θ_lab = (180 − θ_cm)/2. That is
what `acceptance_C14.C:44` computes and what FRESCO tabulates, so the three agree; the macros assert
it numerically rather than trusting it, because the (d,p) arm shipped a mirrored table for exactly
this reason.

## Missing inputs

1. **The elastic FRESCO curve at 136 MeV.** It is the single number that decides whether the
   decomposition works, and no calculation was supplied. `decompose_C17.C` therefore scans
   R = N_elastic/N_217 rather than assuming it. Obtaining this is the highest-value thing anyone
   can add to this study.
2. **The (d,d') distributions.** The proposal asserts "a similar amount for deuterons" with nothing
   behind it. The deuteron day also carries (d,p), (d,t), (d,³He) and breakup as competition.

## Files

| file | what |
|---|---|
| `inel_kinematics_C17.C` | the analytic map and the go/no-go gate; `plots/kinematics_C17inel.png` |
| `C17_inel_sim.C` | the simulation, both channels, selected by argument |
| `accumulate_C17inel.sh` | one sample, all stages, resumable |
| `run_C17inel_campaign.sh` | the 12-sample matrix, in two waves, with a PID-checked driver lock |
| `inel_summary_C17.C` | σ(Ex), the method floor, the vertex correction, level separation |
| `decompose_C17.C` | the headline: fixed-position three-component fit at one day of beam |
| `kine_lines_C17.C` | the reconstructed kinematic lines and per-window E_x spectra; `plots/kinelines_*.png`, `plots/kinemap_*.png` |
| `export_gui_C17.C` | dumps the accepted events to `gui/C17_gui_data.json` for the viewer |
| `gui/build_gui.sh` | inlines that JSON into `gui/template.html` -> one self-contained HTML file |
| `fresco/` | the DWBA angular distributions supplied with the proposal |

Reco, genfit, acceptance and Ex resolution reuse the 14C macros unchanged — they already take the
target, ejectile, residual and beam masses as arguments. `check_beam_C17.C` is shared with the
(d,p) arm and gained an ejectile-PDG argument here, because with the proton default its kinematics
check finds only stray secondaries on the (d,d') channel and reports a spurious failure.
