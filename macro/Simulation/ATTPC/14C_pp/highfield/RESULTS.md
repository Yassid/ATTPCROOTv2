# Results — 14C(p,p') at 2.85 / 4 / 7 T with AT-TPC and 2 mm pads

Campaign of 2026-08-28: 3 fields × 2 pitches × 5 levels, 8000 events per sample, in
`/mnt/f/a1954_C14_hf`. Reproduce the tables with

```bash
root -b -q 'summary_hf_C14.C("/mnt/f/a1954_C14_hf")'
root -b -q 'hf_figures_C14.C("/mnt/f/a1954_C14_hf")'
```

θ_lab window 20–90°, χ²/ndf < 5, truth-matched fits. Median and IQR throughout; σ ≡ IQR/1.349.

## The headline: the measurement is not detector-limited

Under the adopted a1954 analysis — one constant beam energy for every vertex — **σ(Ex) is flat
across the entire matrix**, 0.19–0.27 MeV, and the 6.728/7.012 pair is one bump in all six
configurations (separation 0.71–0.76). No pad pitch and no field changes that.

The reason is not the detector. Rebuild Ex from **perfect truth kinematics** under the same
constant-beam-energy assumption and the spread is still 0.25 MeV. The 14C beam loses 11.3 MeV
crossing the chamber — the profile extracted from this simulation's own truth is

    E(z) = 160.30 − 0.01070 z − 5.78e-7 z²   MeV,  z in mm
    E(0) = 160.30, E(500) = 154.81, E(1000) = 149.03

— and the vertex is uniform in z, so one number is wrong by up to ±5 MeV in a way perfectly
correlated with a quantity the detector already measures well. Beam-energy straggling about that
profile is σ = 0.32 MeV, but dEx/dE_beam ≈ 0.05, so the straggling contributes only ~17 keV: the
**method floor is 0.013–0.019 MeV**, essentially nothing.

Evaluating the beam energy at the *reconstructed* vertex is therefore a pure analysis change with
no hardware attached, and it is worth more than any hardware in this matrix.

## σ(Ex) [MeV], the two reconstructions

| config | as analysed (gs) | E_beam(z_reco) (gs) | as analysed (6.728) | E_beam(z_reco) (6.728) |
|---|---|---|---|---|
| b285_attpc (a1954 as run) | 0.265 | **0.176** | 0.188 | **0.047** |
| b285_2mm | 0.221 | **0.077** | 0.201 | **0.035** |
| b400_attpc | 0.216 | **0.074** | 0.189 | **0.036** |
| b400_2mm | 0.232 | **0.049** | 0.202 | **0.032** |
| b700_attpc | 0.214 | **0.053** | 0.197 | **0.040** |
| b700_2mm | 0.230 | **0.043** | 0.199 | **0.031** |

Once corrected, the numbers are detector performance: the elastic channel improves 4.1× from the
worst cell to the best, and 2.3× from switching the pad plane alone at unchanged field.

## Separation of the multiplet, ΔE/(σ_a+σ_b)

| config | 6.094–6.728 | 6.728–7.012 | 7.012–8.317 |
|---|---|---|---|
| | as analysed → corrected | | |
| b285_attpc | 1.68 → 6.44 | 0.76 → **3.02** | 3.50 → 16.5 |
| b285_2mm | 1.58 → 8.87 | 0.71 → **4.23** | 3.33 → 21.8 |
| b400_attpc | 1.72 → 8.26 | 0.74 → **3.97** | 3.45 → 20.1 |
| b400_2mm | 1.59 → 9.77 | 0.71 → **4.61** | 3.32 → 23.3 |
| b700_attpc | 1.60 → 7.69 | 0.74 → **3.45** | 3.42 → 17.2 |
| b700_2mm | 1.58 → 9.93 | 0.72 → **4.71** | 3.23 → 24.4 |

The 284 keV pair goes from unresolvable to comfortably resolvable **on the existing hardware**,
and the pad plane then improves it by a further 40 %.

## Tracking: this is where the hardware shows

σ(KE_reco − KE_true) [MeV], elastic:

| b285_attpc | b285_2mm | b400_attpc | b400_2mm | b700_attpc | b700_2mm |
|---|---|---|---|---|---|
| 0.343 | 0.095 | 0.067 | 0.038 | 0.055 | 0.038 |

A factor 9 between the worst and best cells. Angle resolution moves much less (0.099 → 0.071°),
and at 7 T with AT-TPC pads it gets *worse* for the excited levels (0.14–0.19°): the helix radius
shrinks by 2.5× against 2.85 T and 8 × 12 mm pads under-sample it. **7 T only works with a fine
pad plane.**

## The cost of high field: the curl-up hole

Acceptance is not flat in θ_cm. Low recoil energies (forward θ_cm) curl into small helices that
pattern recognition and the fit lose, and the hole widens with the field:

| field | hole in θ_cm | depth |
|---|---|---|
| 2.85 T | 25–35° | 0.49–0.66 |
| 4 T | 20–30° | 0.16–0.25 |
| 7 T | 20–50° | 0.06–0.26 |

Overall acceptance (elastic): 0.826 / 0.863 / 0.791 / 0.823 / 0.657 / 0.689 in matrix order. 7 T
costs ~20 % of everything, concentrated exactly where an inelastic angular distribution needs
forward θ_cm points. 2 mm pads *gain* a few per cent at every field.

## CORRECTION, 2026-08-28 evening: the angular weighting changes the conclusion

Everything above is computed on samples generated **flat in theta_cm**. That is the right choice
for an acceptance simulation and the wrong one for reading off what a configuration would do to
a1954, because dEx/dE_beam is not a constant:

    theta_lab   20    28    36    44    52    60    68    76    84
    theta_cm   140   124   108    92    76    60    44    28    12
    dEx/dEb  .118  .104  .087  .069  .050  .033  .019  .008  .001

A flat sample puts most of its weight at the top of that range. The a1954 elastic peak sits at
theta_cm 21-38 deg (quartiles) -- at the bottom. Applying the correction to the data
(`pp/ex_ebeamz_C14.C`) accordingly changes the elastic peak by nothing: width 0.146 -> 0.148 MeV,
Ex-vs-z drift +0.010 MeV over the metre before and after. On the cache with no empirical theta
correction the drift does fall from +0.109 to +0.037 MeV, so the correction is real and correctly
signed -- there is simply almost nothing for it to remove where the elastic yield is.

It DOES pay on the inelastic group, which sits at mean theta_cm 80 deg: IQR 0.475 -> 0.411 at
theta_cm 70-95, 0.574 -> 0.494 at 95-140, and on the raw cache the 6.094 peak narrows 27 % with
its centroid moving 113 keV toward the true energy.

**Acceptance re-weighted by the real angular distributions** (fold of each configuration's
acceptance(theta_cm) with the a1954 yield):

| config | flat in theta_cm | a1954 elastic | a1954 inelastic |
|---|---|---|---|
| b285_attpc | 0.826 | 0.783 | 0.883 |
| b285_2mm | 0.863 | **0.881** | **0.910** |
| b400_attpc | 0.791 | 0.594 | 0.851 |
| b400_2mm | 0.823 | 0.616 | 0.896 |
| b700_attpc | 0.657 | **0.227** | 0.737 |
| b700_2mm | 0.689 | **0.260** | 0.807 |

The curl-up hole is centred exactly on the a1954 elastic peak. **7 T would keep a quarter of the
elastic normalisation and 4 T under two thirds**, which the flat-weighted column hides completely.
The inelastic states are barely affected because they sit at 80 deg.

And in the angular region the inelastic yield occupies (theta_lab 35-65 deg), the simulation gives
sigma(Ex) = 0.151-0.160 MeV as analysed and 0.031-0.039 MeV corrected, **in every configuration** --
the hardware makes no difference there either. The data, corrected, gives 0.150 MeV for the same
peak. So a1954 carries roughly 0.15 MeV of excitation-energy width that the simulation does not
model at all, and that term is larger than every effect in this matrix.

## Recommendation

(revised after applying the correction to the data -- read the correction section above first)

1. **Find the ~0.15 MeV the simulation does not contain.** In the angular region of the inelastic
   yield every configuration predicts 0.031-0.039 MeV after the vertex correction and the data
   delivers 0.150. That gap dominates everything else measured here.
2. **Apply the vertex beam energy anyway.** Free, correctly signed, worth 13-27 % on the width of
   the states the analysis measures and 113 keV of centroid on the 6.094 line in the uncorrected
   cache, and it leaves the elastic peak -- and therefore the energy calibration -- untouched.
3. **2 mm pads at 2.85 T.** Better resolution everywhere and better acceptance in BOTH channels
   (0.881 vs 0.783 on the elastic). The "4 T sweet spot" of the first pass was an artefact of the
   flat-in-theta_cm weighting: at 4 T a1954 loses a quarter of its elastic normalisation and at
   7 T three quarters of it.

## Checks that were done rather than assumed

- **Pattern recognition is not what changed at 2 mm.** `cluster_eval_C14.C` on a 4 T sample:
  97.5 % clean recoil-proton clusters, 0 % split, mean purity 0.988. HDBSCAN's ε = 10 mm spans
  five pads at 2 mm instead of one, and it does not matter.
- **genfit's assumed per-hit error is not hiding the fine-pitch gain.** At 7 T + 2 mm, measSigma
  4 → 2 → 1 mm gives σ(Ex) 0.045 → 0.044 → 0.043 MeV and identical efficiency
  (`meassigma_check.sh`).
- **Diffusion was measured, not assumed.** Magboltz for H2 at 300 torr, 293 K, 50 V/cm, B ∥ E:
  v_d and D_L are field-independent (0.436 cm/µs; 406/377/387/397 µm/√cm at 0/2.85/4/7 T) while
  D_T falls 490 → 167 → 124 → 73. Anchored on the a1954-tuned CoefT = 9e-4 at 2.85 T, that gives
  σ_T at full drift of 3.72 / 2.75 / 1.64 mm at 2.85 / 4 / 7 T — i.e. only at 7 T does diffusion
  stop being the thing that limits a 2 mm pad.

## Caveats

- The corrected numbers use an energy-loss profile fitted to this simulation's truth. On data the
  profile has to come from an energy-loss calculation (CATIMA is in the tree) or from the data
  itself; the shape is a smooth quadratic over 1 m, so this is not a hard step, but it is a step
  that has not been taken yet.
- The a1954 vertex-z resolution enters through E(z_reco). Here it costs nothing — the
  E_beam(z_true) and E_beam(z_reco) columns agree to a few keV — because dEx/dE_beam is small.
- Relative level intensities are not modelled; the summed spectra in `hf_multiplet_*.png` use
  equal weights.
