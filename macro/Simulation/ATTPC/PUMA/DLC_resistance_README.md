# PUMA DLC: measured sheet resistance → PRF σ

This branch already spreads charge over a resistive pad plane via
`AtPulse::SetChargeDispersion(sigma_mm)` (the `fPRFSigma` grid PRF, used by the
`prfSigma` argument of `run_digi_ukf_genfit_test8.C` and the sweep scripts).

What was missing is the link from the **measured DLC sheet resistance
(1.2–1.5 MΩ/□)** to that σ. On a resistive anode the avalanche charge disperses
(Dixit/telegraph) as:

```
sigma = sqrt( 2 t / (R C) )
  R = DLC sheet resistance      [Ohm/square]   = 1.2-1.5e6 (measured)
  C = areal capacitance DLC<->pads [F/mm^2]
  t = charge spreading time     [us]           ~ shaping/peaking time (0.5)
```

## Added
- **`AtPulse::SetChargeDispersionFromDLC(R, C, t_us)`** — sets `fPRFSigma` from the
  physical DLC properties (feeds the same PRF machinery as `SetChargeDispersion`).
- **`dlc_prf_sigma.C`** — prints the σ↔R table so you can pass the right `prfSigma`
  to the existing macros/sweeps.

## Result (t = 0.5 µs)
| R [MΩ/□] | σ (d=50 µm ERAM) | σ (T2K-empirical) |
|---|---|---|
| 1.20 | 1.18 mm | 1.83 mm |
| 1.35 | 1.11 mm | 1.72 mm |
| 1.50 | 1.05 mm | 1.63 mm |

σ_DLC ≈ 1.0–1.8 mm (central ~1.1 mm) → charge over ~1–2 pads (2 mm pitch). At the
4 T field gas transverse diffusion is suppressed, so this DLC term dominates.

## Usage
```cpp
atPulse->SetChargeDispersionFromDLC(1.35e6, 6.0e-13, 0.5); // -> sigma ≈ 1.11 mm
// or, equivalently, pass prfSigma ≈ 1.1 to run_digi_ukf_genfit_test8.C
```

## The one thing to confirm
`C` (areal capacitance) is set by the DLC insulator gap/dielectric and is the
only remaining unknown (swings σ by ~2×). Default `6.0e-13 F/mm²` ≈ 50 µm insulator,
ε_r≈3.4. Set it from the detector build, or from a measured RC, via the `C` argument.
