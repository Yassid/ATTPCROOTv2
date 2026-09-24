# C15p — 15C(p,p') and 15C(p,d)14C, a2091 proton-target runs

a2091 runs **0138–0182** are the **proton-target** set: 15C beam on H2 at 300 torr. "Proton or
deuteron" is never a property of the runs — it is the **ejectile**. Both channels live on the same
PID plane: `(p,p')` elastic/inelastic, and `(p,d)` neutron pickup leaving 14C.

## Pipeline

```
reco_batch_C15p.sh      raw HDF5 -> AtPSAMultiFit + AtDirDeDxCleaner + HDBSCAN + AtPIDTask
ic_batch_C15p.sh        raw HDF5 -> unpackFRIB -> icsum  (writes evtid: SEE THE TRAP BELOW)
overnight_hdb_C15p.sh   the whole chain up to the PID plane, for drawing gates
pid/open_gate_draw.sh   the gate GUI       (OUT=, POINTS=, EBEAM=, REFD=, LOCUS_NBX/NBY)
regate_C15p.sh          refit a redrawn gate and print the g.s. / cluster fits
pp/build_explorer_C15p.sh   the explorer, COMPLETE (see the second trap)
```

Constants: **Ebeam 198 MeV** (confirmed; the elastic ridge independently gives 199.8–200.4),
B = −2.85 T, geometry `ATTPC_H300torr_RT`, par `ATTPC.a2091_C15.par`, gas density 3.308e-5.

## ★ Trap 1 — the IC join, and why it came back

`pid/gate_events_C15p.C` writes a reco holding only the gated events, **renumbered from zero**. A
fit on that input therefore reports a *gated-file index* as its event number, and every downstream
join on `(run, event, trackID)` mismatches. The cache still fills and the spectrum still looks
plausible, so nothing announces it.

**This has happened twice**: C15d (commit `0d220fc3`, 92 % of 353,860 tracks hidden) and then C15p
(62 % of 18,264), because each workspace keeps its own copy of these macros and a new workspace is
ported from an older, unfixed one.

The fix is five pieces and **all of them are needed**:

| file | what it must do |
|---|---|
| `pid/icsum_C15p.C` | write `evtid` from `AtRawEvent::GetEventID` |
| `pid/gate_events_C15p.C` | `SetEventID` before filling; join the IC on `evtid` |
| `AtReconstruction/AtFitterTask.cxx` | propagate the id through the fit (framework, already in) |
| `dump_kine_C15p.C` | read `GetEventID`, falling back to the tree index for ungated input |
| `pid/make_points_C15p.C` | join on `evtid`, size arrays by the largest event number |

Two guards now make a recurrence visible instead of silent — **keep them when porting**:

* `pp/kin_pp_C15p.C` prints `IC matched : N of M (x %)` and raises a red **IC JOIN BROKEN** alarm
  below 50 %, naming the three things to check.
* `pp/make_explorer_pp_C15p.C` omits the `ic` **and** `npulse` columns together when fewer than
  half the tracks carry a value, and says why. Dropping only `ic` is what let the stale
  multiplicity cut through last time.

`npulse == 0` does **not** mean "no beam" — it means the join never reached that track. The page's
default `npLo = npHi = 1` then discards those silently, which is how a viewer ends up showing a
quarter of its sample.

## ★ Trap 2 — the explorer loses its panels on every rebuild

`pp/make_explorer_pp_C15p.C` produces a **bare** page. The IC panel, the per-run Ex map, the
Ebeam(vz) knob and the control defaults are four *separate* post-processing steps, each silently
lost on a rebuild while the page still looks fine. **Always rebuild through
`pp/build_explorer_C15p.sh`**, which runs all five stages and reports each.

## ★ Trap 3 — do not tighten a gate in Brho

Brho is momentum, so a Brho cut is a cut on ejectile energy and therefore on Ex. Measured here:
tightening to `Brho > 0.898` halved the sample and displaced the g.s. by 0.6 MeV while nearly
doubling its width. A genuine species cut narrows the peak and **leaves the centroid alone** —
watch the g.s. mean, not just the width. For contamination, use fit-quality levers instead:
`chi2ndf`, `maxVtxR`, a minimum `nclusters`.

To identify which band a gate encloses, compare its locus occupancy to the **plane baseline**
((p,p') 28.3 %, (p,d) 4.7 %), never the raw percentage — the raw numbers look ambiguous because
protons are abundant enough that their locus sweeps most of the plane.

## Known gaps

* **run_0180** reconstructs 1 event from a 26 GB raw file, and failed identically in the previous
  analysis. Its `/frib` side is fine (51,790 events), so the failure is purely on `/get`. Cause
  unknown; ~36,000 events.
* **run_0151** the FRIB DAQ stopped early (evtid 0..38,593 against 39,993 reco). `evtid` recovers
  96.5 % of it; the trailing ~1,400 events honestly carry `ic = -1`.
* The **Ebeam(vz) knob** in the explorer defaults to 11.0 MeV/m, which is the a1975 16C(d,p)-in-D2
  value. It has **not** been measured for 15C+p in H2 at 300 torr — `a1975/D2_UKF/exvz_dp.C` is the
  macro that measures it.
