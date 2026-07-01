# AtPSA hit-extraction map — for the multipeak upgrade (2026-06-28)

Goal: resolve MULTIPLE peaks per pad (two tracks crossing the same pad at different drift
times) so overlapping/backward tracks don't lose hits = the point-cloud continuity issue.

## Framework is ALREADY multipeak-ready (the bottleneck is only the leaf peak-finder)
- `AtPSAtask` -> `AtPSA::Analyze(rawEvent, event)` (AtPSA.cxx:127): loops pads, calls
  `AnalyzePad(pad)`, adds ALL returned hits to the AtEvent, records per-pad multiplicity.
- `virtual HitVector AnalyzePad(AtPad*) = 0;`  HitVector = vector<unique_ptr<AtHit>>.
- Each hit: (x,y) = pad coord; **z = CalculateZGeo(peakTB)** (AtPSA.cxx:90); charge.
  => N peaks on a pad -> N hits, same (x,y), different z. Pipeline handles it end to end.

## Leaf analyzers — what each produces (AtReconstruction/AtPulseAnalyzer/)
| Class            | how it extracts | hits/pad | state |
|------------------|-----------------|----------|-------|
| **AtPSAMax**     | std::max_element over [20,end-12] -> single global peak | 1 | PRODUCTION (a1975 macros). LOSES 2nd track on shared pads = the continuity issue |
| AtPSASpectrum    | TSpectrum::SearchHighRes -> 1 hit/peak | many | **BROKEN: floatADC.fill(0) (line ~44) zeroes the trace BEFORE SearchHighRes -> searches zeros**; per-peak charge also reads the zeroed/centroid array |
| AtPSAFull        | intervals above threshold, long ones sliced into 50-TB chunks | many | not true peak separation; short-branch leaves charge=0 |
| AtPSADeconv      | FFT deconvolve response (/R(k)+Butterworth) -> sharp current "Qreco"; chargeToHits loops getZandQ() (z,q) vector | multipeak-READY | base `getZandQ` returns the single charge-weighted CENTROID (1 hit) |
| AtPSAIterDeconv  | iterative deconv (sharpens further) then chargeToHits | inherits getZandQ | same: 1 hit unless getZandQ overridden |

a1975/D2 pipelines instantiate **AtPSAMax** (AtPSASimple2 is commented out).

## Insertion points for multipeak (best -> quick)
1. **Override `AtPSADeconv::getZandQ(trace)`** (AtPSADeconv.cxx:263) to find MULTIPLE peaks on
   the DECONVOLVED current and return one ZHitData per peak. Deconvolution un-smears the
   electronics response that merges nearby peaks -> cleanest two-track separation. Downstream
   (chargeToHits:235 -> Analyze) already loops the vector. **Right place; most physics-correct.**
2. **Fix `AtPSASpectrum`**: delete the `floatADC.fill(0)` so TSpectrum searches the real trace,
   and assign each peak its local charge (not the total). Quick multipeak-on-raw to validate.

## Validate with
The viewer (../../../Simulation/ATTPC/16C_dp_gnn/viewer/) + compare pre/post point-cloud
continuity; then re-check clustering (Spyral/GNN) on the improved clouds.
