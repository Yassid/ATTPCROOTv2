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
```

⚠ `ddbuild.C` and the rest write/read `dd_ex.root` in the CURRENT directory; run them from a
scratch dir (they were developed in `~/Data/C15d_statefit`).

⚠ The pp reference must come from `pp/plots_pgate/` (proton-GATED, 142k tracks), not
`pp/plots_hdb/` (3.3k) — the latter has too few tracks for a peak search and its "peaks" are noise,
including one below the ground state.
