# Reproducing the a1954 12Be(p,t)10Be analysis from scratch

Both fitters, on a fresh machine. Everything below has been run; the numbers in
§7 are what you must land on. Companion to `doc/a1954_Be12_guide.tex` (the full
12Be guide) — this file is only the (p,t) path.

**What the answer is.** There are no 10Be levels in this data set. The yield sits
on a constant-triton-KE ridge (~16 MeV), so the Ex spectrum shows a fake ~6 MeV
bump and the median Ex *slides* with angle at −0.57 MeV/deg (both fitters; §7).
A real state is **flat** in Ex-vs-θ. If you reproduce a peak, check the slope
before you believe it. Reproducing the pipeline is the goal; a discovery is not.

---

## 1. Build (≈ half a day)

```bash
export FI=$HOME/fair_install && mkdir -p $FI && cd $FI
# 1. FairSoft + FairRoot  (standard ATTPCROOT instructions, install to $FI/FairSoft/install
#    and $FI/FairRootInstall)
# 2. GenFit -- MUST be the fork with the CATIMA dE/dx backend
git clone -b catima-scattering https://github.com/Yassid/GenFit.git    # c54d9cb
cmake -S GenFit -B GenFit/build -DCMAKE_INSTALL_PREFIX=$FI/GenFit && cmake --build GenFit/build -j8 --target install
# 3. ATTPCROOT fork, branch OpenKF-Claude
git clone -b OpenKF-Claude git@github.com:Yassid/ATTPCROOTv2.git ATTPCROOTv2-OpenKF
cmake -S ATTPCROOTv2-OpenKF -B ATTPCROOTv2-OpenKF/build -DGENFIT=$FI/GenFit && cmake --build ATTPCROOTv2-OpenKF/build -j8
```
CATIMA itself is pulled by CMake (FetchContent); nothing to install by hand.
Verify GenFit really has it — an empty result means the fits will silently run
with **zero** stopping power below βγ = 0.05 (KE < 3.5 MeV for a triton):

```bash
nm -DC $FI/GenFit/lib/libgenfit2.so | grep -ci catima     # must be > 0
```

Environment, every session (`config.sh` does **not** set `ROOT_INCLUDE_PATH`, and
cling crashes without it):

```bash
export REPO=$FI/ATTPCROOTv2-OpenKF
export HERE=$REPO/macro/Unpack_HDF5/a1954_Be12/UKF
source $FI/FairSoft/install/bin/thisroot.sh
source $REPO/build/config.sh
export ROOT_INCLUDE_PATH=$REPO/build/include:$FI/FairRootInstall/include
```

Geometry (once per clone, not tracked in git):

```bash
cd $REPO && root -b -q geometry/ATTPC_H600torr.C
```

## 2. Inputs and settings that are not negotiable

| | |
|---|---|
| raw data | `run_XXXX.h5`, a1954 remerged; runs **0142–0224** (0160, 0177–0181, 0183 absent) |
| gas | H2 **600 torr**, 293 K, ρ = **6.616e-5 g/cm³**, geometry `ATTPC_H600torr` |
| par file | `parameters/ATTPC.a1954_Be12.par` (B 2.85 T, drift 1.00 cm/µs, 3.0 MHz, 512 tb) |
| IC beam window | **500–800** on the FRIB IC amplitude, tb 1050–1250 |
| `minPoints` | **30**, and it must be identical in the PID dump and in the gating step |
| `Ebeam` | 155 MeV, anchored on the (p,p') elastic peak — see §8 |

**The gas really was 600 torr** (confirmed 2026-09-25). A 2026-08-25 rebuild
wrongly assumed 300 torr (3.308e-5, `ATTPC_H300torr_RT`) — half the material, and
live the moment CATIMA/matEffects is on: it pushed the (p,p') elastic peak to
+0.2–0.3 MeV. The drivers now pass the gas to both fitters **explicitly**
(`GASDENS`, `GEONAME`); do not rely on the macro defaults, which still say 300 torr.
(The 14C setting of a1954 *was* 300 torr — do not copy its gas across.)

## 3. Stage 1 — reconstruction, FRIB, slim cache

```bash
cd $HERE
./reco_hdb_batch.sh "run_0143 run_0144 ..." 2     # unpack + AtPSAMultiFit + cleaner + HDBSCAN
./frib_batch.sh     "run_0143 run_0144 ..." 1     # ion chamber -> <run>_FRIB.root
./slim_cache_batch.sh "run_0143 ..." 4            # AtPatternEvent-only cache, ~250x faster reads
```
Edit the input/output paths at the top of each script for your machine. The slim
cache (`<run>_slim.root` + `<run>_FRIB.root` in one directory) is the only thing
the rest of the analysis reads. Cost: hours per run for reco, minutes for the rest.

## 4. Stage 2 — the triton gate

The (p,p') gated inputs are **proton-only** — every triton was thrown away there,
so (p,t) must be re-gated from the full reco.

```bash
root -b -q "pid/dump_pid_Be12.C(\"run_0143,run_0144,...\",\"$SLIM\",\"pid/pid_Be12_data_mp30.csv\",\"$SLIM\",500,800,1050,1250,2.85,false,30)"
root -b -q 'pid/make_points_Be12.C("pid_Be12_data_mp30.csv","plots/points_Be12_mp30.root")'
root -l   'pid/gate_draw_pt_Be12.C'      # interactive: needs a display, never -b
```
The drawer overlays the 12Be(p,t) g.s. Brho(θ) locus plus the existing proton and
deuteron gates. **Use [Locus check]**: a gate in (√(dE/dx), Brho) is only
trustworthy if what it selects lands on the kinematic curve. Saves to
`pid/triton_12Be.json`. The committed gate is deliberately broad — the IC window,
not the polygon, does the beam selection (contaminant beam near IC 1900 outnumbers
12Be ~20:1 inside the gate).

## 5. Stage 3 — gate + both fits, one command

```bash
./fitpipe_pt_Be12.sh "run_0142 run_0143 ... run_0224" 4
```
75 runs, ~21 min at 4-parallel (~6 min at 6-parallel when the gated inputs already
exist). Output goes to `$PTDIR`, default `~/a1954_Be12_fit_pt_600torr/`. Per run it
does three things:

1. `pipeline/gate_events_Be12.C` — IC 500–800 + triton polygon + `minPoints=30`
   → `$PTDIR/in/<run>_reco.root`
2. `pipeline/fitUKF_Be12.C(run,-1,"triton",-1,2.85,6.616e-5,...)` → `<run>_ukf.root`
3. `pipeline/fitGenfit_Be12.C(run,-1,in,"",out,-2.85,2,5,"",4.0,3.0,170.0,kTRUE,kFALSE,"triton","ATTPC_H600torr")`
   → `<run>_genfit.root`

Note in (3): `bField = −2.85` (experimental handedness), θ window **3–170°**
(tritons are forward — the proton default of 10–170 slices the band), and
`matEffects = kTRUE` with the macro defaults `catimaELoss=kTRUE`,
`catimaELossFull=kFALSE`, `matFallback=kFALSE`. Those defaults are the point of
the whole exercise; do not override them. `run_0161` has an empty FRIB file and
is skipped automatically. Check any log: UKF prints `gasDensity : 6.616e-05 g/cm^3`,
and GENFIT prints a red WARNING if its geometry is not a 600-torr one.

**Gate-reuse stamp.** The script writes `<run>.gatecfg` (IC window, tb, θ, mp, gate
md5) and reuses a gated input only when that string matches. If you change a knob
that is not in `gatecfg()`, the change is a silent no-op.

## 6. Stage 4 — excitation energy and the explorer

```bash
R="run_0142,run_0143,...,run_0224"
D=$HOME/a1954_Be12_fit_pt_600torr/
root -b -q "pp/ex_Be12.C(\"$R\",\"$D\",155,5.0,\"_pt600_ukf\",3.016049,10.013534,\"12Be(p,t)10Be\",\"ukf\")"
root -b -q "pp/ex_Be12.C(\"$R\",\"$D\",155,5.0,\"_pt600_genfit\",3.016049,10.013534,\"12Be(p,t)10Be\",\"genfit\")"
./pp/open_explorer.sh pt     # both fitters in one browser page (GAS=300 for the superseded set)
```
Masses: triton 3.016049 u, 10Be 10.013534 u. χ²/ndf cut 5.0. Caches land in
`pp/plots/proton_kin_pt600_{ukf,genfit}.root` (the `pt800` caches are the superseded
300-torr set).

The page has an **Ebeam(vz)** control (off by default): Ebeam = 155 − loss·(vz − 500 mm),
loss default **9.1 MeV/m** = CATIMA for 12Be at 155 MeV in 600-torr H2. The beam enters
at low vz (the (p,p') elastic peak rises with vz). 9.1 MeV/m removes only ~60 % of that
drift — the (p,p') data would need ~14 MeV/m to go flat, and the excess is unexplained.

## 7. Did it work? — the checks, in order

| check | expected |
|---|---|
| tracks after χ²/ndf < 5 | **2169 UKF / 2631 GENFIT**, no collapsed fits |
| gate cleanliness | 0 tracks shared with the proton gate, 1 with the deuteron |
| IC window matters | 500→800 vs 500–900 gains ~18 % of tritons at equal purity |
| **median Ex vs θ_lab** | slope **−0.57 (UKF) / −0.57 (GENFIT) MeV/deg** — the real result |
| control: (p,p') | run `fitpipe_Be12.sh` too; elastic peak **+0.003 / −0.010 MeV**, σ 0.210 / 0.181 |

How the slope is measured: median Ex in **5° θ_lab bins from 5 to 45°**, bins with
≥ 20 tracks, straight-line fit. Binning moves it: the same 300-torr production gave
−0.42/−0.48 on the old uneven bins and −0.60/−0.55 on these, so compare like with like.
The slope is **gas-insensitive** (300 torr gave −0.60/−0.55 on the same bins; (p,d)
gives −0.32/−0.31 at both pressures) — fast ejectiles barely lose energy in the gas.
The elastic peak is a gaussian over ±1.2σ around the maximum of a 30 keV-binned Ex
histogram, Ebeam 155.

The (p,p') control is what proves the reconstruction is sound: it gives a real,
flat peak while both transfer channels slide. If your (p,t) is flat and (p,p') is
not, something upstream is wrong, not the physics.

**Use medians.** Per-angle-bin gaussian fits rail on about half the bins here
(+38.98 ± 19.74 and similar) and produce a slope built on garbage.

## 8. Traps that have each cost a day

1. **Density.** 600 torr / 6.616e-5 / `ATTPC_H600torr`. Wrong values are inert
   with matEffects off and catastrophic with it on — which is exactly this analysis.
   At 300 torr the (p,p') elastic peak sits at +0.265 / +0.210 instead of ~0.
2. **`minPoints` mismatch.** A gate drawn on an mp15 plane and applied at mp30 is a
   different cut, silently. `AtSpyralPID::fMinPoints` defaults to 30.
3. **Ebeam = 155 MeV is in-sample.** It was found in July by scanning until the
   (p,p') elastic peak sat at zero, at 6.5e-5 g/cm³ — i.e. essentially the real
   600-torr gas. So the peak landing at ~0 now is *consistency*, not an independent
   measurement of the beam energy. dEx/dEbeam = +0.0088 MeV/MeV for (p,p'), so the
   elastic peak constrains Ebeam only weakly. Quote Ex on this scale with that caveat.
4. **θ window.** 3–170° for tritons. The proton 10–170 cuts the most forward tritons,
   where the cross section is largest.
5. **matEffects co-requisites** — correct geometry, `matFallback=kFALSE`, no manual
   dE/dx, and a dE/dx source (CATIMA). Each fails silently on its own.
6. **The two fitters disagree by 0.012 MeV** on the (p,p') elastic peak at 600 torr
   (it was 0.05–0.085 at the wrong 300 torr). Quote it as a systematic; it is still
   ~4× the statistical error on either.

## 9. Where the reference output lives (this machine)

```
$HERE                                   macros, drivers, gates
~/a1954_Be12_reco_hdb_slim/             slim + FRIB inputs (76 runs)
~/a1954_Be12_fit_pt_600torr/            the (p,t) production (75 runs, ukf + genfit)
~/a1954_Be12_fit{,_pd}_600torr/         the (p,p') and (p,d) productions
$HERE/pp/plots/proton_kin_pt600_*.root  kinematics caches
~/a1954_Be12_pt_explorer_600torr.html   browser explorer
~/a1954_Be12_fit_pt/, *_pt800_*         SUPERSEDED 300-torr production, kept for comparison
~/a1954_Be12_pt_share/                  collaborator bundle — STILL THE 300-TORR SET
~/a1954_analysis_runs/2026-08-25_Be12_300torr_rebuild/   logs, params, plots of the campaign
```
