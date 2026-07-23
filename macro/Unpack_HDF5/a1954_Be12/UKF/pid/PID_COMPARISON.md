# a1954 ¹²Be PID — Spyral 1.0.0 vs ATTPCROOT, HDBSCAN vs Triplclust

Overnight study (2026-07-21), runs 0147–0150 for the 3-way comparison. Goal: the cleanest
¹²Be(p,p′) PID plane (√dE/dx vs Bρ).

## TL;DR — recipe for the cleanest PID

**HDBSCAN clustering  +  arclen > 200 mm  +  IC(¹²Be) gate 625–750.**

- `arclen > 200 mm` is the single biggest lever — it turns a smeared continuous band into
  clean, separated compact blobs (the smear was short tracks with unreliable dE/dx).
- HDBSCAN edges out Triplclust; Spyral 1.0.0 and ATTPCROOT-HDBSCAN are comparable.

## Quantitative band tightness (fractional FWHM of the proton band; lower = cleaner)

| sample | Spyral 1.0.0 | ATTPCROOT HDBSCAN | ATTPCROOT Triplclust |
|---|---|---|---|
| all beam            | 0.485 | **0.422** | 0.413 |
| 12Be-gated (IC 625–750) | 0.272 | **0.215** | 0.279 |
| 12Be-gated + arclen>200 | 0.111 | **0.109** | 0.147 |

(The all-beam Triplclust number edges HDBSCAN, but on the physics-relevant 12Be-gated
sample HDBSCAN is consistently tighter. The metric is a coarse peak-FWHM; trust the trend.)

## Findings

1. **Clustering — HDBSCAN > Triplclust** for PID band definition. `praType="hdbscan"`
   (mover join) in `unpackReco_Be12.C`. Triplclust leaves the band a bit more diffuse.

2. **Framework — Spyral 1.0.0 ≈ ATTPCROOT-HDBSCAN.** Same band quality. The ATTPCROOT
   in-framework `AtSpyralPID` reproduces Spyral's PID shape.

3. **⚠ dE/dx SCALE MISMATCH.** Spyral's dE/dx ≈ **4× ATTPCROOT's** (√dE/dx ≈ 2×: Spyral
   median ~26, ATTPCROOT ~13). A normalization difference (charge/arc-segment convention),
   not a bug — but it means **a PID gate from one framework does NOT transfer to the other**.
   Build the gate on whichever framework you'll analyze in. The **IC amplitude scale DOES
   match** (~1890 dominant in both), so the beam gate is portable.

4. **arclen > 200 mm resolves two clean populations** in the all-beam plane (√dE/dx≈5–8,
   Bρ≈0.6 and √dE/dx≈13–22, Bρ≈0.2). `npts` cuts add nothing beyond arclen.

5. **Statistics-limited.** ¹²Be is only ~7% of the beam cocktail → ~2000 gated tracks / 4
   runs, ~850 after arclen. Overnight reco of +16 runs (0143–0163) is running to sharpen it.

## How to reproduce

```bash
# ATTPCROOT HDBSCAN reco (praType=hdbscan) -> _reco.root
./reco_hdb_batch.sh "run_0147 run_0148 ..." 2
# dump PID to CSV (IC recorded from _FRIB.root)
root -b -q 'pid/dump_pid_Be12.C("run_0147,...","/home/yassid/a1954_Be12_reco_hdb/","/tmp/pid_hdb.csv","/home/yassid/a1954_Be12_reco/",-1,-1)'
# Spyral 1.0.0
cd ~/attpc_spyral_1.0.0 && .venv/bin/python run_spyral_a1954_12Be.py 147 150 8
# 3-way comparison (IC gate + arclen)
python3 ~/attpc_spyral_1.0.0/compare_pid_final.py 625 750 200 _12Be
```

## Open items
- Identify which of the two arclen-cut blobs is protons vs deuterons (needs kinematics /
  UKF hypothesis test).
- Optionally match the ATTPCROOT `AtSpyralPID` dE/dx normalization to Spyral's (×4) so gates
  become portable.
- Rebuild with the full +16-run statistics (in progress).

## dE/dx scale mismatch — root cause (diagnosed)
The AtSpyralPID estimator formula is IDENTICAL to Spyral's: dEdx = (charge summed over the
first-arc inner big-pad segment, r<152mm) / (1000-pt spline arclength). So the ~4× dE/dx
difference is NOT in the estimator — it is in the **per-point charge** fed to it: ATTPCROOT's
`AtPSAMultiFit` hit charge vs Spyral's pointcloud charge (peak amplitude vs integrated
charge convention). To make gates portable, match the PSA charge definition (≈×4). Until
then, build/apply PID gates within a single framework.

## Trace-integral charge experiment (result)
Added opt-in `AtSpyralPID::SetUseTraceIntegral(true)` (uses AtHit::GetTraceIntegral instead of
GetCharge). Effect on √dE/dx median: amplitude 11.6 -> integral 20.0 (Spyral 25.9) — recovers
most of the scale gap (residual ~1.3x is charge calibration). BUT band tightness is UNCHANGED
(fracFWHM 0.075 both) because fractional width is scale-invariant. Conclusion: charge definition
affects only the absolute scale / gate portability, NOT band cleanliness. Cleanliness is set by
clustering (HDBSCAN) + arclen>200 cut. Flag is default-off; rebuilt in build_genfit.

## FINAL high-statistics result (35 runs, 0142-0185)
Reco'd all 20 runs with HDBSCAN (+FRIB). Combined PID: 118,814 tracks; 11,611 12Be-gated
(IC 625-750); 5,592 with arclen>200. The 12Be-gated proton band is now clean and
well-defined (`plots/pid_Be12_HIGHSTATS.png`).

**Built a 12Be proton PID gate** `pid/proton_band_Be12.json` (AtParticleID format, ATTPCROOT
AtSpyralPID scale) — wraps the proton band, 58.7% of IC-gated tracks in-gate
(`plots/proton_gate_Be12.png`). NOTE: this gate is on the ATTPCROOT amplitude-charge scale;
it does NOT transfer to Spyral (use SetUseTraceIntegral + rescale, or rebuild on Spyral).

## Bottom line / recommendation for the cleanest 12Be PID
1. Reco with **HDBSCAN** (`praType="hdbscan"`), multipeak PSA.
2. Beam gate: **IC 625-750** (12Be), event-ID matched from `<run>_FRIB.root`.
3. Track cut: **arclen > 200 mm** (drops short tracks with unreliable dE/dx).
4. Particle gate: **proton_band_Be12.json**.
Spyral 1.0.0 gives equivalent band quality; ATTPCROOT-HDBSCAN matches it and keeps everything
in one framework.

## Updated final numbers (35 runs, 0142-0185)
Reco'd 35 HDBSCAN runs + FRIB. Combined PID: **204,520 tracks; 15,755 12Be-gated (IC 625-750);
7,502 with arclen>200**. Definitive plots: `plots/pid_Be12_FINAL.png`,
`plots/proton_gate_Be12_FINAL.png`. Refined gate `pid/proton_band_Be12.json` (28 vertices,
58.8% of IC-gated tracks in-gate). NOTE: `dump_pid_Be12.C` segfaults on near-empty runs
(run_0161, 6MB) — exclude such runs; a GetEntries()<10 guard would harden it.
