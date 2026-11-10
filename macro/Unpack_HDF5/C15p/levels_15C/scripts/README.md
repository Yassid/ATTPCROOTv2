# Scripts moved here from C15d/statefit on 2026-11-09 — these are 15C LEVELS, seen in both
# 15C(p,p') and 15C(d,d'), so they live with the 15C+p analysis.

# (d,d') state fit — scripts

Reproduces `../DD_STATE_FIT.md`. `exlib.h` re-implements the explorer's transforms exactly, so
offline numbers match what is on screen; the two JSONs are the saved explorer parameters the fit
uses verbatim.

```bash
R=/home/yassid/fair_install/FairSoft/install/bin/thisroot.sh; source $R
root -b -q ppstates2.C      # 15C(p,p') spectrum, g.s. position vs Ebeam
root -b -q ppfind2.C        # its inelastic structures + significance  -> the reference states
root -b -q ddbuild.C        # both (d,d') spectra with the saved parameters -> dd_ex.root
root -b -q statefit2.C      # the fit: background-only vs +2 peaks, fixed then free positions
root -b -q syst.C           # fit-range systematic on the state energies
root -b -q statefig.C       # the figure
root -b -q psrun.C          # 3-body phase space templates -> ps_alt.root  (needed by globalfit.C)
root -b -q globalfit.C      # THE GLOBAL FIT: g.s. + 5/2+ + A + B + background, 3 background forms
root -b -q plotglobal.C     # the global figure -> global_ex_fit.png
```

`globalfit.C` is the one to read: it fits the whole spectrum for both codes and prints the state
energies with a background systematic taken from the spread over three background parametrisations.
Its header documents the three single-minimisation attempts that FAILED first and why, so they are
not retried.

⚠ `plotglobal.C` trap: **`TF1::Draw` is lazy** — the function is evaluated when the canvas is
painted (at `SaveAs`), not when `Draw` is called. Any draw-only function that reads a global will
be painted with whatever that global holds at the END of the macro, so the first panel silently
gets the second panel's parameters. Every draw function there carries the background mode, Kb and
sigma_res as explicit TF1 parameters for this reason.

⚠ `ddbuild.C` and the rest write/read `dd_ex.root` in the CURRENT directory; run them from a
scratch dir (they were developed in `~/Data/C15d_statefit`).

⚠ The pp reference must come from `pp/plots_pgate/` (proton-GATED, 142k tracks), not
`pp/plots_hdb/` (3.3k) — the latter has too few tracks for a peak search and its "peaks" are noise,
including one below the ground state.
