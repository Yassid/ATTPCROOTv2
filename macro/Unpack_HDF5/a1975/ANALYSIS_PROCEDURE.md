# a1975 — common analysis procedure

The one procedure every a1975 channel follows: `16C(p,p)`, `16C(p,d)`, `16C(p,t)`, `16C(d,t)`.
It exists because the channels kept diverging in ways that were invisible until they produced a
wrong number — a missing ion-chamber cut, a window inherited from a different target, a run
silently truncated to one event. Each rule below is here because something went wrong without it.

Anything marked **MEASURED** was determined from the data and the measurement is named. Anything
marked **ASSUMED** has not been verified and should be treated as an open item.

---

## 0. Working point

Source it, never retype it:

```bash
source macro/Unpack_HDF5/a1975/working_point.sh
```

| setting | value | note |
|---|---|---|
| beam energy, (p,d) | 185.0 MeV | |
| beam energy, (d,t) | 184.25 MeV | with the kinematics correction |
| geometry | `ATTPC_H300torr_RT` | 3.308e-5 g/cm3 — **not** `ATTPC_H1bar` (8.27e-5, 2.5x too much material) |
| material effects | **ON everywhere** | CATIMA MSC + straggling + native dE/dx |
| B field | 2.85 T | sign: +1 in simulation, −1 in data |

The geometry default is the dangerous one: with material effects OFF the density cannot reach the
fit, so a wrong geometry is harmless and stays unnoticed. Turning material effects on makes it live.

## 1. Run blocks

**MEASURED** (from the production run lists):

| target | runs | channels |
|---|---|---|
| D2 | 0016–0103 (47 runs) | `(d,t)` |
| H2 | 0106–0189 (84 runs) | `(p,p)`, `(p,d)`, `(p,t)` |

Never carry a number derived on one block over to the other without re-measuring it. The IC window
is the standard example — see §3.

## 2. The five stages

Every channel is the same chain. Each stage has a check that must pass before the next one runs.

| stage | produces | how to tell it worked |
|---|---|---|
| 1. unpack | `<run>_reco.root`, `<run>_FRIB.root` | both exist and open cleanly |
| 2. genfit production (`prod_*.sh`) | `<run>_genfitter_*.root` + `.done` | banner reports CATIMA MSC/straggling ON, correct geometry, back-extrapolation ON |
| 3. cache (`cache_*_run.C`) | flat tree, one row per track | entry-count check silent (§4); per-run counts non-zero |
| 4. PID gate | `*.json` polygon | scored against the sim plane, not reused across planes |
| 5. analysis | spectra, cross sections | §6 |

**The `.done` marker is the completion guard, not file existence.** A killed job leaves a
short, openable ROOT file; `[ -s file ]` calls it present. Use `pd/root_ok.C`, which opens the
file, checks it closed cleanly and holds a TTree, and rejects `kRecovered`.

## 3. The ion chamber — the one selection every H2 channel shares

**The IC is not a stored variable.** It lives only in `<run>_FRIB.root`, inside `AtRawEvent`, as
8 generic traces x 2048 ADC samples. Reaching the one number we want (max ADC in a time window)
costs a full deserialisation: ~1 GB per run, ~88 GB over the H2 block, ~45 min sharded 7 ways.
Pay it **once**, into a cache, and gate from the cache thereafter.

The amplitude is `max(trace[0].ADC[b])` for `b` in `[1000, 1350)`. **MEASURED**: the pulse sits at
tb ~1150 (`macros/ic_peaks_pt.C`, full-length trace profile), so the window contains it.

### What the beam actually is

**MEASURED** over all 84 H2 runs, 1.95M events (`macros/ic_spectrum_pt.C`). The IC measures dE/dx,
so peak positions must scale as Z^2 with one constant. Anchoring the tallest as carbon predicts
the rest with no free parameters:

| Z | element | predicted | observed | diff | fraction |
|---|---|---|---|---|---|
| 3 | Li | 289 | 288 | −0.6% | 1.4% |
| 4 | Be | 514 | — | absent | — |
| 5 | B | 804 | 778 | −3.3% | 0.9% |
| 6 | **C** | 1158 | **1158** | anchor | **61.3%** |
| 7 | N | 1575 | 1578 | +0.1% | 1.1% |
| 8 | O | 2058 | 2058 | −0.0% | **33.1%** |

A third of the beam is oxygen. The carbon valleys sit at **950 and 1380**.

### Windows in use

| window | used by |
|---|---|
| **1000–1350** | `pp/angdist_gs_pp.C`, `pp/xsec_pd.C` — the cross-section chain |
| 950–1350 | `pd/ex_pd_a1975.C`, `pid/pid_ic_a1975.C` — the PID / Ex macros |
| 900–1400 | the `(d,t)` chain, on **D2** runs |

Both H2 windows sit inside the carbon peak and neither clips it (valleys at 950/1380). They are
nonetheless different numbers, and which one a macro uses should be deliberate.

**RULE: the same window on both sides of a ratio.** The `(p,d)` normalisation is the elastic yield
against DWBA, so the window enters numerator and denominator. **MEASURED**: both use 1000–1350,
so it cancels — verified in `angdist_gs_pp.C` and `xsec_pd.C`.

**The IC selects Z, not A.** It removes oxygen and the light contaminants. It cannot separate 16C
from other carbon isotopes, and it does not remove reactions on the 16C beam itself — protons and
deuterons from `(p,p)` and `(p,d)` leak into a triton PID gate regardless of the IC.

## 4. Event alignment — and the check that must be there

**MEASURED**: entry `i` of the FRIB tree is the same event as entry `i` of the reco and genfit
trees. There is no event number to join on — the FRIB `AtRawEvent` carries a sequential ID, but
`AtPatternEvent` and `AtTrackingEvent` both report **−1**. Index order is the only handle.

Every macro that joins the two streams does `N = min(fit, FRIB)`, which means **a short FRIB file
does not fail — it silently truncates the run**. `cache_pp_run.C`, `cache_pd_run.C` and
`macros/cache_pt_run.C` now report:

- **FRIB shorter** — red warning with the count and percentage dropped. This is a bug.
- **FRIB longer** — a note. Six H2 runs end with a junk entry (event ID −1); it sits at the END,
  so indices still align for every real event. Harmless.
- **equal** — silent, so healthy runs do not fill a driver log.

### Known bad run

**`run_0148`: its FRIB tree holds 1 entry against 23514 in reco.** Every event gets `ic = −1`,
so any IC gate discards the whole run. **MEASURED**: it is the only run in 106–189 with zero
entries in `pd_kin_catima.root` (147 → 49, 148 → 0, 149 → 113). It is already absent from both
the `(p,d)` yield and the elastic calibration by the same mechanism, so the ratio is unbiased —
but the effective H2 run list is **83, not 84**, and that should be stated wherever it matters.

## 5. Cache conventions

Two tree names, fixed, so downstream tools work on any channel unchanged:

| tree | one row per | branches |
|---|---|---|
| `pts` | pattern track (no fit needed) | `sqrtdedx, brho, polar, vertexz, vertexr, ncl, ic, run` |
| `pk` | fitted track | `ke, theta, vz, chi2ndf, ic, run, kefit, thetafit` |

- `ke`/`theta` are the **back-extrapolated** slot (`GetKinematicsXtr`), `kefit`/`thetafit` the raw
  fit (`GetKinematics`). Identical for any production run without back-extrapolation.
- `polar` in `pts` is the AtSpyralPID convention: **theta_lab = 180 − polar**. `theta` in `pk` is
  already theta_lab. Getting this backwards mirrors the data about 90 deg.
- `chi2ndf = 1e9` is the collapse sentinel. Drop it explicitly; plotted as data it looks like a
  60% background.

A cache that costs a FRIB pass should be built by a **script in the repo**, never an inline
`root -e`. An inline command is not reproducible and silently drops columns — that is exactly how
the `(p,t)` chain lost its `ic`.

## 6. Before believing a result

- **Medians, not means.** A heavy fit-failure tail once faked an entire "genfit is broken" result.
  Check any large ratio against a control and a quality cut.
- **Acceptance binning at least as fine as the yield.** A coarse bin once let a zero-acceptance
  bin be corrected by x2.60.
- **A fitted parameter sitting exactly on its bound is not a measurement.** The `(d,t)` g.s.
  centroid was railed at every limit, with a 0.001 MeV error on a 250-count peak.
- **`SetRangeUser` beyond the axis draws the OVERFLOW bin** and makes a good fit look catastrophic.
- **Plot components with the total and the continuum**, or a fine fit looks broken.
- **A peak that walks with angle is not a state.** Slice in theta_lab and check the centroid is
  stationary. This is what rejected `(p,t)`: an ~11.5 MeV walk from 8 to 40 deg, unchanged by the
  IC gate.
- **`(Brho, theta_lab)` identifies the reaction**, the PID plane does not. A two-body channel is a
  single curve there; `macros/pt_brho_theta.C` overlays the loci of every competing channel.

## 7. Where things live

| path | contents |
|---|---|
| `/mnt/f/a1975/` | **all** a1975 data — `reco/`, `reco_*_catima/`, `caches/`, logs |
| `/mnt/g/` | a1954 (14C, 12Be) + gas data — **not** a1975 |
| `/mnt/c/Users/Yassid/Desktop/` | plot and explorer copies for viewing |

## 8. Open items

- **ASSUMED**: the elastic cache side of §3's cancellation was verified from shared code, not
  measured — `pp_kin.root` was on `/tmp` and is gone. Rebuilding it would close this.
- **MEASURED**: the `(p,d)` simulation is over-gained ~2.5x (117 vs 46 hits/track against run_0106).
  `(d,t)` uses a measured gain of 35000; `(p,d)` inherited 150000.
- Spyral and ATTPCROOT differ by 2–3x in normalisation; Spyral has no channel-specific acceptance.
- `min(fit, FRIB)` truncation is **reported but not changed**, so existing caches keep their
  meaning. It matters only to an un-gated analysis.
