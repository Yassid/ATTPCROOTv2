# AT-TPC December 2014 α runs — analysis and simulation manual

Everything needed to reproduce the α+α analysis from raw `.graw` files to an excitation
function, and to run the simulation that validates the reconstruction.

Branch: **`FairRootv18.00-fairroot19-port`**, worktree
`/home/yassid/fair_install/ATTPCROOTv2_fr19port`.

> **Work on this branch only.** `AtEvent` is `ClassDef(3)` here and `ClassDef(6)` on
> `OpenKF-Claude`. Opening a file written by the other branch yields **zero hits with no
> error message at all** — not a warning, not an exception, just empty events. If every
> event looks empty, check which branch wrote the file before debugging anything else.
> This cost a full day once.

---

## 1. Setup

### 1.1 Environment

```bash
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port
source setup_fr19port.sh
```

This sets `VMCWORKDIR`, `SIMPATH`, the ROOT/FairRoot paths, **and the Geant4 data
variables** (`G4ENSDFSTATEDATA`, `G4LEVELGAMMADATA`, …). Those last ones are set *nowhere*
in `build/config.sh`; without them every simulation dies with
`G4ENSDFSTATEDATA environment variable must be set` before the first event.

### 1.2 Build

```bash
cd build && cmake .. && make -j8
```

The event display is built PCL-free: `AtEventDisplay` is no longer gated behind
`PCL_FOUND`, and the RANSAC/Hough drawing paths were removed from `AtEventDrawTask`. You
do **not** need PCL to visualise events.

### 1.3 Python

The analysis scripts need `numpy`, `scipy`, `matplotlib` and **`sklearn` ≥ 1.3** (for
`sklearn.cluster.HDBSCAN`, which is part of scikit-learn — *not* the standalone `hdbscan`
package, which is not installed and is not needed).

The system `python3` (3.8) has numpy and matplotlib but **no scipy and no sklearn**. Use:

```bash
PY=/home/yassid/attpc_spyral/.venv/bin/python     # sklearn 1.5.0, numpy 1.26.4, scipy 1.13.1
```

Any Python ≥3.9 environment with those four packages works; nothing here is Spyral-specific.

### 1.4 Directory layout

| Path | Contents |
|---|---|
| `/media/yassid/Cris_OneD/Dec2014_alphas/run_XXXX/` | raw `CoBo*.graw` (read-only source) |
| `~/dec2014_links/run_XXXX/` | symlink trees (scratch, regenerable) |
| `$VMCWORKDIR/runfiles/NSCL/Dec2014_alphas/` | run lists — **not tracked in git**, regenerate |
| `~/dec2014_alphas_reco/lowP/` | unpacked hits, runs 128–139 (3.1 GB) |
| `~/dec2014_calib/` | analysis scratch: CSVs, pickles, plots |
| `macro/Unpack_GETDecoder2/calib/` | all analysis scripts (tracked) |

The scripts read and write intermediate products in `~/dec2014_calib`. Create it and run
from there. Script *imports* are self-locating, but several data paths are still hardcoded
to that directory — grep for `dec2014_calib` if you relocate.

---

## 2. Detector parameters and where they came from

`parameters/ATTPC.alpha_150torr.par`. The values that were **measured in this analysis**
rather than inherited:

| Parameter | Value | Provenance |
|---|---|---|
| `DriftVelocity` | 2.25 cm/µs | 2.251 ± 0.011 (run 100, B=0) and 2.253 (run 128) from the 90° opening angle |
| `TiltAng` | 6.47° | run 100 with the magnet **off**, where the beam direction measures the tilt directly |
| `TBEntrance` | 320 | beam's closest approach to the detector axis (misses by 2.2 mm). **Was 280** — an error of 40 tb = 144 mm in *z* |
| `EField` | 12000 V/m | from the beam geometry. The cathode (15 kV / 1 m) implies 15000; the beam axis-crossing strongly prefers 12000. The ~20 % gap is within GARFIELD-level accuracy and is accepted |
| `BField` | 0.5691 T | run log |
| `ThetaPad` | 110.9° | pad-plane rotation; a fit of the Lorentz shear gives 112 |
| `ThetaRot` | −161.9° | azimuth of **B** in the pad plane, from the run-100 beam |
| `SamplingRate` | 6.0 | → `GetTBTime()` = **160 ns** per time bucket |

Two gotchas in this file:

- **`TB0` (98) is marked DEPRECATED and must not be used as the drift-time origin.** See §3.4.
- **`Gain: 0.28` is not a gain.** It is harmless when unpacking data (unused) but makes
  simulated ADC ~20000× too small. The simulation therefore uses a separate file.

### Parameter sets

| File | Use |
|---|---|
| `ATTPC.alpha_150torr.par` | experimental data, 149–150 torr (runs 128–139) |
| `ATTPC.alpha_150torr_sim.par` | **simulation only** — `Gain: 5000`, `CoefL/CoefT: 0.0009` |
| `ATTPC.alpha_300torr_noB.par` | magnet-off calibration runs 96–113, 299.5 torr |

The diffusion coefficients differ because the experimental file says `0.025` labelled
`[cm^-0.5]` while `AtClusterizeTask` consumes the number as `[cm²/µs]`. The label is wrong,
not the code; the simulation file carries a value in the units actually used.

---

## 3. The Lorentz / Langevin correction

This is the heart of the reconstruction. Read this section before changing anything.

### 3.1 Why it is needed

The detector is tilted by ~6.5° so that the beam does not run down the pad-plane normal.
That means **E** and **B** are *not parallel*. Ionisation electrons therefore do not drift
straight down **E**: the **E**×**B** term pushes them sideways, and a track is reconstructed
sheared away from its true direction. The shear grows linearly with drift time, so hits far
from the pad plane are displaced most.

### 3.2 The equation

From the AT-TPC commissioning paper (Bradt *et al.*, NIM A **875** (2017) 65, eq. 4):

```
v_D = v/(1+ωτ²) · [ Ê + ωτ (Ê × B̂) + ωτ² (Ê·B̂) B̂ ]        with  ωτ = (B/E)·v
```

Two departures from the published form, both necessary here:

1. **The paper assumes B lies in the y–z plane** (tilt azimuth = 90°). For the Dec 2014
   data the beam — and hence **B** — sits at about **−162°** in the pad plane. The azimuth
   is therefore kept general and read from `ThetaRot`. Passing `azimDeg = 90` reproduces
   the published expression exactly.

2. **The result must be rotated into the pad frame by `ThetaPad`.** The Langevin solution
   is derived in the *field* frame, but it is applied to *pad* coordinates. Omitting this
   rotation was the original bug: the shear came out with roughly the right magnitude
   pointing about 110° away from where it belonged, which looks like a plausible-but-wrong
   correction rather than an obvious failure.

### 3.3 Where it lives

`AtParameter/AtLangevin.h` — header-only, `AtTools::LangevinDrift(...)`. **One** definition,
used by both directions of the calculation:

| Direction | Caller | Action |
|---|---|---|
| forward (simulation) | `AtClusterizeTask` | displaces electrons by **+v_xy · t** |
| reverse (reconstruction) | `AtPSA::CalculateXCorr/YCorr` | subtracts the same |

These were once two independent copies. That is dangerous, not merely untidy: the
simulation's entire claim to validating the correction is that truth and reconstruction
disagree by nothing. If the copies drift apart, the residual silently measures the
difference *between the copies* instead of testing the physics. Keep the single helper.

### 3.4 The drift-time origin

The transverse and longitudinal corrections must agree about where *z* = 0 is:

```cpp
Double_t AtPSA::DriftTimeUs(Int_t tb) const {
   Double_t tbZero = fEntTB - fZk * 100.0 / (fTBTime * fDriftVelocity);
   return (tb - tbZero) * fTBTime * 1E-3;
}
```

`CalculateZGeo` puts *z* = 0 at that `tbZero`; `fTB0` (98) is **not** that value, which is
why it is deprecated. With the calibrated parameters `tbZero = 42.22`.

Note the origin only **translates** the event — it cannot rotate it. It moves the vertex,
never a track direction. An earlier hypothesis that `TB0` caused the over-rotation was
disproved numerically: the beam tilt is identical to four decimal places for `TB0` from
−100 to +200.

Sign matters: drift time is `T = tb − tbZero`, **not** `tbZero − tb`. The inverted form
goes negative for beam hits (tb up to 289) and produced a bogus v_D of 4.05.

### 3.5 Numbers for the calibrated parameter set

With v = 2.25 cm/µs, B = 0.5691 T, E = 12000 V/m, tilt 6.47°, azimuth −161.9°, `ThetaPad`
110.9°:

```
ωτ         = 1.0671
v_D        = (0.18272, −0.02463, 2.23479) cm/µs      [pad frame]
|v_xy|     = 0.1844 cm/µs
shear      = atan(|v_xy| / v_z) = 4.716°
full drift = 44.75 µs  →  82.5 mm transverse displacement
```

The reconstruction logs exactly these numbers once per run:

```
==== Lorentz drift vector (pad frame) : (0.182718, -0.024625, 2.23479) cm/us,  omega*tau = 1.06706
```

### 3.6 How it was validated

**Against Monte-Carlo truth**, commit `8f9d645e`. Simulated α+α events go
truth → Langevin drift → digitisation → *the same PSA the experimental data uses*, so the
residual is measured against a truth no experiment can supply.

| drift [mm] | N tracks | uncorrected | corrected |
|---|---|---|---|
| 255 | 40 | 3.60° | 1.79° |
| 595 | 94 | 3.70° | 1.31° |
| 935 | 143 | **4.73°** | **0.22°** |

What makes this convincing is the *quantitative* prediction: a transverse v_xy alongside a
drift v_z tilts an apparent track by exactly `atan(v_xy/v_z) = 4.716°`, so the uncorrected
error must converge to that at full drift. It measures 4.73° — 0.3 %. The correction then
removes it to 0.22°.

Improvement shrinks at short drift (2.0× at 255 mm vs 21.5× at 935 mm) simply because there
is less shear to remove. That is proportionality to drift time, not a weakness.

> **Use track DIRECTION, not hit-to-point distance, to test a position correction.** A hit
> sits at a pad centre, so for an oblique track the true point can be most of a pad away
> however perfect the correction is. That granularity floor made the hit-distance test look
> like a failure at middle drift (residuals of 60–90 mm against a shear of only 26–36 mm).
> Direction averages the granularity down over many hits and is the observable the physics
> rests on anyway.

> **Fix the pass criterion before looking at the numbers.** Here: the uncorrected error had
> to *grow* with drift, and the corrected error had to be small. Stating that in advance is
> what made the result trustworthy — two earlier rounds got the *z* convention backwards and
> then misread the outcome as a code bug (see `71df2f1f`, which retracts a drift-length
> "bug" that never existed). If the advance prediction fails, the setup is wrong; do not
> reinterpret the numbers.

### 3.7 Re-checking the correction on real data

The correction actually written into the production files can be verified in closed form —
dump the hits (§4.4) and compare `PositionCorr − Position` against the prediction:

```python
tbZero = 320 - 1000*100/(160*2.25)          # = 42.22
t      = (tb - tbZero) * 160e-3             # µs
dx_pred, dy_pred = -0.182718*t*10, -0.024625*t*10   # mm
```

On run 128 this closes to **max |applied − predicted| = 0.0011 mm over 19 465 hits**.

### 3.8 Known limitation

`ThetaPad` and `ThetaRot` are **not individually identified** — only their net effect is
constrained by the data. Changing `ThetaPad` alone silently breaks the correction. Resolving
the degeneracy needs the pad-map geometry and has not been done.

---

## 4. Experimental chain

### 4.1 Build the symlink trees and run lists

`AtGRAWUnpacker` maps files to decoders with a `%i` pattern, so the inputs need contiguous
`file<N>_` indices. It also matters that chunk 0 comes first: plain `ls` sorts `".1.graw"`
*before* `".graw"`, which desynchronises the per-CoBo event counters and **silently kills
PSA** — you get a file with no hits and no error.

```bash
cd $VMCWORKDIR
for r in $(seq -w 128 139); do bash macro/Unpack_GETDecoder2/calib/make_links.sh run_0$r; done
```

Prints `<runfile path> <number of CoBos>`; expect **9** CoBos for these runs. The list is
written to `$VMCWORKDIR/runfiles/NSCL/Dec2014_alphas/alpha_run_XXXX.txt`, which is where the
unpack drivers look. Override `ATTPC_RAW` / `ATTPC_LINKS` if your paths differ.

> Run lists are **not tracked in git** — always regenerate them. They contain absolute
> paths to the symlink trees.

### 4.2 Unpack

```bash
bash macro/Unpack_GETDecoder2/calib/unpack_lowP.sh ATTPC.alpha_150torr.par 128 129 130
```

The driver counts events first with `count_events.C`, then unpacks. **The event count must
be exact**: asking for more events than the run holds crashes at end-of-data (`rc=129`) and
*truncates the output*.

Output: `~/dec2014_alphas_reco/lowP/alpha_run_XXXX_hits.root`, ~20 MB per 1000 events.
`AtRawEvent` is not persisted (~30 kB/event instead of ~1.3 MB); you get hits only.

Expect **9216 pads/event**. If you see 2304, you are hitting the AsAd-merge bug — that is
fixed on this branch and validated digit-for-digit against the 2014 reference analysis.

Runs 128–139 are already produced (3.1 GB, ~15 min/run, all `rc=0`):

| run | events | run | events |
|---|---|---|---|
| 128 | 11361 | 134 | 13103 |
| 129 | 11880 | 135 | 12776 |
| 130 | 15778 | 136 | 10594 |
| 131 | 15734 | 137 | 14337 |
| 132 | 12915 | 138 | 14674 |
| 133 | 10523 | 139 | 10115 |

### 4.3 Look at events

```bash
bash macro/Unpack_GETDecoder2/calib/open_eve_fr19.sh ~/dec2014_alphas_reco/lowP/alpha_run_0128_hits.root
```

Needs a display. `run_eve_Dec2014_alphas.C` is the underlying macro.

### 4.4 Dump hits to CSV

The python analysis reads CSV, not ROOT:

```bash
cd $VMCWORKDIR/macro/Unpack_GETDecoder2        # so ./rootlogon.C is picked up
root -l -b -q 'dump_hits.C("'$HOME'/dec2014_alphas_reco/lowP/alpha_run_0128_hits.root","'$HOME'/dec2014_calib/prod_128.csv","hits",-1)'
```

Columns: `event,pad,tb,x,y,z,xc,yc,zc,q,qtot` where `x,y,z` are raw pad coordinates,
`xc,yc,zc` the Lorentz-corrected ones, and:

- **`q` = `GetCharge()`** — pulse *peak height*
- **`qtot` = `GetQHit()`** — pulse *integral*

Use `qtot` for anything involving deposited energy or total charge. Mode `"qtot"` writes one
summed `GetQHit()` per event, which is what the trigger study consumes.

### 4.5 Cluster the tracks

Plain HDBSCAN on (x, y, z) merges the two outgoing arms **at the vertex** — exactly where you
least want it to, because the vertex is a real density maximum. `cluster_tracks.py` therefore
clusters in a *direction-aware* space: position augmented by the sign-invariant orientation
tensor **d**⊗**d** of each hit's local neighbourhood. Two arms crossing at a vertex have
different orientations even where they touch, so they separate.

```bash
$PY macro/Unpack_GETDecoder2/calib/batch_scan.py     # clusters run 128, writes candidates.pkl
```

Clustering runs **once per event** at a reference drift velocity, and the assignment is then
reused while v_D is scanned — otherwise HDBSCAN would run 11k times per scan point.

> **Watch for isolated noise hits.** One stray hit (`tb=224, x=−51.5, q=22`) gave a 22 mm
> arm a fake 142 mm extent and moved the fitted v_D from 2.8 to 4.6. Cluster-size and
> extent cuts are not cosmetic here.

### 4.6 Calibrate

```bash
$PY macro/Unpack_GETDecoder2/calib/calibrate_run.py all_run0100.csv --B 0.0 --E 30000 --tilt 7.0
$PY macro/Unpack_GETDecoder2/calib/fit_params.py
```

The method: for equal-mass non-relativistic elastic scattering the two outgoing α's always
separate by **90°**. The pad coordinates are fixed by geometry, but the *z* coordinate scales
with v_D — so v_D is the value that puts the opening-angle peak at 90°.

Restrictions that matter:

- The 90° constraint holds **only for elastic scattering off He**. Scattering off the C and
  O in the CO₂ populates smaller opening angles, so fit the *peak position*, never the mean.
- **Use a magnet-off run (96–113) to measure the tilt.** With B on, the solenoid steers the
  beam, so the beam direction describes the beam, not the detector. For the same reason the
  axis-crossing check is only meaningful with B off.
- **Regress x and y *on* drift time**, and average slope *components* before forming an
  angle. A 3D total-least-squares line fit inflates the angle (6.89° vs 6.40°) because pads
  are 8×12 mm while time is finely sampled; and taking a median of per-event angles inflates
  it again (6.96° vs 6.42°) because `hypot()` is positive-definite and noise cannot cancel.

### 4.7 Verify the production

```bash
$PY macro/Unpack_GETDecoder2/calib/verify_production.py
$PY macro/Unpack_GETDecoder2/calib/plot_production.py
```

Checks `GetPositionCorr()` per run: beam polar angle, opening-angle peak, vertex position.

### 4.8 Excitation function

The beam enters at 7.8 MeV and slows continuously, so **the vertex position *is* an energy
measurement** — a single run scans E_cm from 3.9 MeV down to 0. That is the whole point of
thick-target inverse kinematics.

Energy loss comes from CATIMA (ported from OpenKF-Claude) for **He:CO₂ 90–10** at 150 torr,
tabulated in `calib/eloss_alpha_heco2.txt`.

```bash
$PY macro/Unpack_GETDecoder2/calib/excitation.py
$PY macro/Unpack_GETDecoder2/calib/kinematics.py     # acceptance-independent closure test
```

`kinematics.py` is the consistency check worth understanding: acceptance changes how densely
events populate a kinematic locus, but **it does not move the locus**. So E₁ = E cos²θ₁ and
θ₁+θ₂ = 90° can be tested even though the efficiency is unknown.

> The elastic selection is not optional. Omitting it leaves a sample dominated by split
> tracks and θ₁+θ₂ comes out at 27° instead of 90°.

---

## 5. Simulation chain

`macro/Simulation/ATTPC/Dec2014_alphas/`.

### 5.1 Geometry and gas

- Medium `heco2_150` in `geometry/media.geo` — He:CO₂ 90/10, 7.27e-5 g/cm³
- `geometry/ATTPC_HeCO2_150torr.C` generates `ATTPC_HeCO2_150torr.root`

> The geometry generator **hangs in `CheckOverlaps` after writing its files** — just kill it.

### 5.2 Generate

```bash
cd $VMCWORKDIR/macro/Simulation/ATTPC/Dec2014_alphas
root -l -b -q 'He4He4_sim_el.C(4000)'
```

α+α elastic, beam 4He at 1.95 MeV/u (p_z = 0.241296 GeV/c per nucleon), field 5.691 kG.
~23 s for 4000 events.

> `AtTPC2Body` via `AtVertexPropagator` **alternates beam-only and reaction events**, so a
> 4000-event run gives 2000 reactions. Outgoing energies well under half the beam energy are
> also expected, not a bug: vertices are spread along a beam that has already slowed, so a
> late vertex gives products of only tens of keV.

#### Running in parallel

For the statistics the trigger study needs, use the driver rather than a loop:

```bash
bash bulk_sim.sh 8 5000 ~/dec2014_sim_bulk     # 8 streams x 5000 = 40000 events
```

It generates, digitises and dumps the per-event charge in one pass, then concatenates the
streams into `qtot_sim_bulk.txt`. Two reasons it cannot be a simple `for` loop in one
directory:

- **Both macros write `./data/attpcsim_in.root` by name.** Streams sharing a working
  directory overwrite each other's files mid-run. Each stream gets its own directory.
- **The default seed comes from `time(NULL)`**, which has one-second granularity — streams
  launched together would draw *byte-identical* events and the extra CPU would buy no extra
  statistics at all. `He4He4_sim_el.C` now takes an explicit seed as its third argument
  (`0` keeps the old time-seeded behaviour), and the driver gives each stream a distinct
  one. The seed is echoed as `==== Generator seed : N`, so a run can be reproduced.

> **Do not put `set -u` in a script that sources `setup_fr19port.sh`.** That script, and
> FairSoft's `thisroot.sh` underneath it, reference unset variables as a matter of course,
> so `set -u` kills the shell *at the source line* — before anything is echoed. The symptom
> is a script that produces no output whatsoever and creates none of its directories.

### 5.3 Digitise

```bash
root -l -b -q 'rundigi_sim.C("./data/attpcsim_in.root","./data/digi.root",kFALSE)'
```

`clusterize → pulse → AtPSASimple2` — **the same PSA the experimental data goes through**,
with the same threshold (20). That identity is what makes the comparison meaningful.

The third argument is `keepElectrons`. Leave it `kFALSE`: persisting the drifted electrons
costs ~6 MB **per event**. With it off, 4000 events digitise in ~7 min into 82 MB.

Uses `ATTPC.alpha_150torr_sim.par` (see §2).

### 5.4 Validate the Lorentz correction against truth

This is the run that produced the table in §3.6 — the only test of the correction that has
access to a truth the experiment cannot supply.

```bash
root -l -b -q 'resdir.C("./data/digi.root")'
```

`resdir.C` compares **track directions**, binned by drift distance: reconstructed direction
vs the MC truth direction, with and without the correction. Read §3.6 before interpreting
the output — in particular, do not substitute a hit-to-point distance metric, and fix the
pass criterion (uncorrected error must *grow* with drift) before looking.

`res.C`/`res2.C`/`res3.C` are earlier position-residual versions kept for reference;
`res4.C` tests whether the middle-drift residual is an artefact of labelling a hit by
`MCSimPoint[0]` when several MC points hit the same pad. It is not: `<points/hit> = 1.0`,
so the weighted centroid equals `MCSimPoint[0]` exactly, and that hypothesis is refuted.

### 5.5 Trigger efficiency

The simulation has **no trigger**, so comparing total charge per event between simulation and
data isolates exactly that: where the data is deficient relative to simulation, the trigger
was firing inefficiently.

```bash
root -l -b -q 'dump_hits.C(".../digi.root","'$HOME'/dec2014_calib/qtot_sim.txt","qtot",-1)'
root -l -b -q 'dump_hits.C(".../alpha_run_0128_hits.root","'$HOME'/dec2014_calib/qtot_exp_128.txt","qtot",-1)'
$PY macro/Unpack_GETDecoder2/calib/trigger_eff.py
```

Result (4000 simulated events vs all 11361 of run 128):

| ΣQ [ADC] | exp / sim |
|---|---|
| 7.1e3 | 0.16 |
| 1.2e4 | 0.38 |
| 2.4e4 | 0.54 |
| 3.9e4 | 0.80 |
| 6.2e4 | 0.99 |

Two things fell out that were not put in: the **absolute charge scale agrees to ~7 %** in the
median even though the simulated gain is an untuned guess; and above 4e5 ADC the data shows a
tail the simulation cannot produce — **beam pile-up**, since the generator makes one reaction
at a time.

Normalisation is taken on the bulk (1.2–3.0e5), *not* on total area, precisely so the pile-up
tail cannot distort the low-charge comparison.

> **This is not yet a usable efficiency curve.** Between 3e4 and 1.3e5 the ratio sits at
> 0.75–0.85 rather than climbing cleanly to 1, and above 1e5 it scatters between 0.5 and 2.0
> on only ~2000 effective simulated reactions. Whether the sub-unity plateau is real or an
> artefact of the normalisation window needs more statistics. **Do not use it to correct
> anything yet.** ~40 000 events (about an hour) should settle it.

---

## 6. Traps

Ordered by how much time each one cost.

1. **Cross-branch `AtEvent` ClassDef 3 vs 6 → silently zero hits.** No error of any kind.
2. **A ROOT macro's entry function must match its filename.** Copying a macro and renaming
   the file makes ROOT fall back to `.L` and do *nothing*. Bit us three separate times.
3. **Libraries must be loaded from `rootlogon.C`, not `#include`d.** Loading from inside a
   macro function is too late — cling parses the whole function first and every class comes
   out as "unknown type name". `AtTpcPoint` (MC truth) lives in **`libAtTpc`**; omitting it
   makes the macro fail to parse and produce *no output at all*, which reads like an empty
   loop rather than an error. `libAtSimulationData` does not exist on this branch.
4. **Geant4 data paths are set nowhere in `build/config.sh`** — always source
   `setup_fr19port.sh`.
5. **`ls` orders `.1.graw` before `.graw`**, desynchronising CoBo event counters and killing
   PSA silently. `make_links.sh` handles it; do not hand-roll the list.
6. **Overrunning the event count crashes and truncates** the output file.
7. **`cd ... && nohup ... &` backgrounds the whole compound**, so the `cd` happens in a
   subshell and the job launches from the wrong directory.
8. **`pgrep`/`pkill` patterns match the matcher's own command line.** Bracket them
   (`unpackReco_C15[d]`) and verify with `ps` before believing — or killing — anything.
9. **`FairRunAna::SetInputFile` is gone in FairRoot 18.6** — use
   `SetSource(new FairFileSource(...))`.
10. **`CalcLorentzVector` runs per event, not once.** An unguarded `cout` there fired
    599 923 times in 1000 events. It is now behind a `static bool`.

---

## 7. Open issues

- **Energy reconstruction** is not done. The range route is broken by the `arc_length`
  estimator; the curvature route `p = 0.3qBR` is preferred and the validated simulation can
  now check it against truth momenta.
- **Trigger efficiency** needs ~10× the statistics before it can correct anything (§5.5).
- **Pile-up** is absent from the simulation, so the high-charge region cannot be compared.
- **`ThetaPad`/`ThetaRot` degeneracy** unresolved (§3.8).
- **Bradt thesis eq. 3.14** (detector → beam frame, left-handed) not implemented.
- **Angular distributions per reaction energy** vs published data — blocked on acceptance.
- **GENFIT port** deferred; its precision cannot yet be exploited.

---

## 8. Script reference

`macro/Unpack_GETDecoder2/calib/`

| Script | Purpose |
|---|---|
| `make_links.sh` | symlink trees + run lists (§4.1) |
| `unpack_lowP.sh` | bulk unpack, 149–150 torr runs |
| `unpack_noB.sh` | magnet-off calibration runs 96–113 |
| `unpack_one.sh` / `unpack_all_fr19.sh` | single / bulk unpack |
| `open_eve_fr19.sh` | event display (this branch) |
| `cluster_tracks.py` | direction-aware HDBSCAN (§4.5) |
| `batch_scan.py` | cluster a whole run, collect elastic candidates |
| `fit_params.py` | fit v_D (and ωτ, tilt) from the opening angle |
| `calibrate_run.py` | full calibration chain for one run |
| `opening_angle.py` | v_D from the 90° constraint |
| `langevin.py` | Langevin drift for arbitrary tilt azimuth (python mirror of `AtLangevin.h`) |
| `verify_production.py` | per-run checks on `GetPositionCorr()` |
| `plot_production.py`, `plot_opening.py` | production plots |
| `excitation.py` | excitation function via CATIMA energy loss |
| `kinematics.py` | acceptance-independent closure test |
| `trigger_eff.py` | trigger turn-on from charge spectra |
| `view_event.py`, `analyse_event.py`, `vertex_tangent.py`, `validate130.py` | single-event tools |
| `eloss_alpha_heco2.txt` | CATIMA α energy loss, He:CO₂ 90-10 @ 150 torr |

`macro/Unpack_GETDecoder2/`

| Macro | Purpose |
|---|---|
| `run_unpack_Dec2014_alphas.C` | decode + PSA |
| `count_events.C` | exact event count (required, §4.2) |
| `dump_hits.C` | hits → CSV, charge → txt (§4.4) |
| `run_eve_Dec2014_alphas.C` | event display |

`macro/Simulation/ATTPC/Dec2014_alphas/`

| Macro | Purpose |
|---|---|
| `He4He4_sim_el.C` | α+α elastic generator + Geant4 |
| `rundigi_sim.C` | clusterize → pulse → PSA |
| `resdir.C` | **direction** validation of the correction vs MC truth (§3.6, §5.4) |
| `res.C`, `res2.C`, `res3.C`, `res4.C` | earlier position-residual versions, kept for reference |

Plots live in `~/dec2014_calib/plots/`.
