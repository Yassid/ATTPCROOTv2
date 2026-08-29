# 14C(d,p)15C at 2.85 / 4 / 7 T with AT-TPC and 2 mm pads

Built 2026-08-29 to answer a question the (p,p') matrix could not: that campaign found almost no
gain from finer pads or a higher field, but its recoil protons never pass θ_lab 90°. A (d,p)
transfer at the same beam energy puts its ejectile **backward** — θ_lab 95–125° with 3–8 MeV over
the angular range where the cross section peaks. Does the matrix pay off there?

Products: `/mnt/f/a1954_C14dp_hf` (uniform θ_cm 2–178) and `/mnt/f/a1954_C14dp_back` (θ_cm 2–60,
extra statistics where transfer peaks). Reproduce every table with `summary_C14dp.C`.

## Setup

14C at 161 MeV on **D₂ at 300 torr** (`ATTPC_D300torr_v2`). Everything not gas-related is identical
to the (p,p') campaign so the two channels are directly comparable. Levels: the 15C ½⁺ ground state
and the 5/2⁺ at 0.740 MeV — the only two bound states, S_n = 1.218 MeV.

Density is set to **6.5643e-5 g/cm³, matched exactly to the transport geometry** rather than to
a1975's 0.0661 mg/cm³, which disagrees with its own geometry by 0.7 %. Diffusion from a dedicated
Magboltz D₂ scan (v_d and D_L field-independent; D_T 526 → 197 → 148 → 89 µm/√cm at 0/2.85/4/7 T),
anchored the same way as the (p,p') campaign.

## Two configuration errors found and fixed — both faked "backward doesn't work"

1. **`backwardSeedFix` was off.** `AtGenfitter` seeds from the lowest-z end with momentum toward
   +z, reflecting a backward track into the forward hemisphere; the truth match then drops it. The
   first run had **zero** reconstructed protons below θ_cm 45°. With the flag on, backward tracks
   reconstruct to θ_lab 168° and the angle comes back with median 0.0°, IQR 0.36°.
2. **Gain 35000 starved the 2 mm plane.** A 2×2 mm pad collects 1/24 the charge of an 8×12 mm one.
   Usable protons were 1870 (AT-TPC) vs 848 (2 mm); at gain 150000 they are 898 vs 897. A pitch
   comparison at that gain measures the threshold. **Design consequence: a 2 mm plane needs ~4×
   the gain, or a proportionally lower threshold.**

Also fixed: a `CM_RANGE_NOT_APPLIED` guard that compared "2 - 60" against "2.0 - 60.0" and failed
every correct run — the same string-vs-number trap the 46Ar script documents.

## Kinematics (`dp_kinematics_C14.C`)

Reconstructed plane per configuration with the two-body loci over it, drawn parametric in θ_cm
(a transfer locus can be double-valued in θ_lab) and at three beam energies (entrance, mid, far
end) so the beam-loss band is to scale. Two things read directly off it: the scatter off the locus
is all at forward θ_lab and tightens with field, and **the two level loci are separated by less
than the line width** — the doublet is a displacement smaller than the reconstruction scatter
almost everywhere.

## Acceptance (`dp_acceptance_C14.C`)

The overall acceptance is 0.63–0.71 in every cell, which makes the field look nearly free. Per
θ_cm it is not:

| θ_cm | 5–10 | 10–20 | 20–30 | 30–45 | 45–60 | 60–80 | 80–110 |
|---|---|---|---|---|---|---|---|
| 2.85 T, AT-TPC | 0.68 | 0.82 | 0.77 | **0.50** | 0.79 | 0.78 | 0.75 |
| 7 T, AT-TPC | 0.00 | 0.63 | **0.18** | **0.08** | 0.54 | 0.82 | 0.80 |

**Two separate features.** A dip at θ_cm 30–45 is present at *every* field, and the 7 T curl-up
hole sits on top of it at θ_cm 20–50, with total loss below 10°.

I first attributed the field-independent dip to θ_lab ≈ 90°, the AT-TPC's blind direction where a
track perpendicular to the drift axis spans almost no z. **The lab-angle acceptance does not
support that** — the hole moves with the field:

| θ_lab | 80–86 | 86–90 | 90–94 | 94–100 | 100–110 |
|---|---|---|---|---|---|
| 2.85 T | 0.82 | 0.67 | 0.66 | **0.22** | 0.72 |
| 4 T | 0.72 | **0.38** | **0.30** | 0.79 | 0.75 |
| 7 T | 0.59 | **0.08** | **0.00** | 0.11 | 0.10 |

A perpendicular-track effect would sit at 90° and not care about the field. This one is at 94–100°
at 2.85 T, 86–94° at 4 T, and spans 86–110° at 7 T — i.e. it tracks a proton *energy* (6–7, 8–10,
4–10 MeV) rather than an orientation. Cause not yet established; the higher-statistics backward
samples are being generated to characterise it. **Do not quote the perpendicular-track
explanation.**

## Resolution (`summary_C14dp.C`)

σ(Ex) for the ground state with the beam energy at the reconstructed vertex:

| θ_cm slice | 2.85 T / 8×12 | 7 T / 2 mm |
|---|---|---|
| 8–30° (backward 125°, 3 MeV) | 1.15 | 0.28 |
| 30–60° (85°, 10 MeV) | 0.58 | 0.069 |
| 60–100° (60°, 22 MeV) | 0.66 | 0.050 |
| 100–180° (34°, 40 MeV) | 1.33 | 0.094 |

The 0.740 doublet goes from unresolvable (separation 0.28–0.65) to 4–8. **The field does this, not
the pads**: σ(KE) at the transfer peak is 0.43 vs 0.46 MeV for the two pad planes at 2.85 T — the
pitch changes nothing — and 0.10–0.14 at 7 T. A backward 3 MeV proton has an 8.8 cm helix radius at
2.85 T, nearly straight over its short length; at 7 T it is 3.6 cm and curls enough to measure.
Fine pads cannot help at 2.85 T because transverse diffusion (3.7 mm at full drift) is already
larger than the pad; only at 7 T (1.6 mm) does the 2 mm plane add anything.

## OPEN: a −20 % energy bias on backward tracks

**SOLVED 2026-08-29 — it is the spiral, and the cause is quantitative.**

At 2.85 and 4 T, protons at θ_lab 92–125° come back with KE about 20 % low. The chain:

1. A backward proton of 3–6 MeV has a **range of 1.5–3 m** in D₂ at 300 torr. I first read that as
   "it cannot stop in a 50 cm chamber" — **wrong**: a spiralling proton travels that path while
   staying inside. It stops.
2. It therefore spirals many turns, shedding a large fraction of its energy on the way. The
   projected hits lie on a **tightening spiral**, not a circle. Measured, window by window from the
   vertex, local radius / true radius at the vertex: **0.964, 0.931, 0.890, 0.843, 0.803, 0.734**.
3. The pattern recognition fits **one circle to all of it** and returns **0.895** — the average.
4. genfit seeds its momentum from that radius and does not recover, so
   **KE_seed/KE_true = 0.895² = 0.80**, since KE ∝ p² at these energies. That is the −20 % exactly.

Everything else was excluded on the way: the vertex correction (`Xtr − raw` = +0.01 to +0.04 MeV,
right sign and negligible), a mis-seeded subpopulation (the distribution is one shifted lump, not
two peaks), and merged clusters — the fat tracks are **100 % pure proton**, purity 1.000, so the
extra hits are the proton's own.

The signature is the hit count, not the angle: at θ ≈ 107° and KE ≈ 5 MeV, a 173-hit track gives
R_geo/R_true = 1.007 and a 1617-hit track gives 0.893. Binned by hit count over all backward
tracks: 170 hits → 0.991, 704 → 0.944, 1471 → 0.888.

At 7 T the bias is ≈ −3 % — the tighter helix stops sooner and the sampled spiral spans less energy
loss — but there the backward sample is also small and acceptance-selected.

## FIXED: the measurement order in AtGenfitter

genfit's `getFittedState()` returns TrackPoint 0 and the back-extrapolation starts there. The
clusters were added in ascending z_lab **whatever the direction of travel**, while the seed is
placed at the vertex — which for a backward track is the *highest* z_lab. Point 0 was therefore the
**stopping end**: the momentum read out was the one the proton had after spending most of its
energy, and the extrapolation ran from the wrong end, integrating the energy loss the wrong way.

The fix adds the measurements in travel order (`addSeq` = `order` reversed when `backwardSeed`).
Forward tracks are byte-for-byte unchanged. Same reconstruction, no truth match applied:

| | θ_lab | median KE_reco/KE_true | within 5 % | fits / truth protons |
|---|---|---|---|---|
| before | 92–110° | 0.834 | 15.7 % | 198 / 267 |
| **after** | 92–110° | **0.999** | **88.4 %** | **233 / 267** |
| before | 110–140° | 0.893 | 39.3 % | 140 / 157 |
| **after** | 110–140° | **0.999** | **90.5 %** | **148 / 157** |

With the truth match: bias/spread **−17.4 % / 19.1 % → −0.1 % / 0.7 %**, while forward 60–85°
stays at −0.3 % / 3.4 %. Backward tracks now measure *better* than forward ones, which is what a
500–2000 hit spiral should give once its curvature is read at the right end.

**Everything below this line was measured with the defect present and is void in the backward
region.** The campaigns are being re-fitted.

### Two candidate fixes tried before the cause was found

Both knobs existed in `AtGenfitter` and neither was reachable from `fitGenfit_C14.C`. Both are now
exposed (trailing arguments, defaults off, no existing caller changes). Same reconstruction each
time, so only the fit differs:

| arm | 60-85 deg forward | 92-110 deg backward | 110-140 deg backward |
|---|---|---|---|
| default seed, matFX on | -0.3 % | **-17.4 %** | -10.2 % |
| default seed, matFX **off** | -0.5 % | **-10.2 %** | -5.9 % |
| + range constraint | -0.3 % | -17.4 % | -9.4 % |
| + **first-arc seed** | -0.1 % | **-17.0 %** | -10.1 % |

**The first-arc seed changes nothing - and that is the important result.** The seed *is* biased
(0.895 of the true radius, 0.80 in energy) but replacing it with a first-arc circle leaves the
answer where it was. So genfit is not sitting at its seed: **the fit itself converges to the biased
value**, and the bias lives in how the hits are fitted, not in how the fit is started.

The pointer to where is the matFX row: switching material effects **off** halves the bias
(-17.4 to -10.2 %) and halves the spread. A constant-momentum helix fitted to a track whose radius
runs 1.00 to 0.73 should land near the middle, about -24 %; it lands at -10 %. Turning the energy
loss on should recover the vertex momentum and instead makes it worse. That is the signature of the
energy loss being integrated along a direction the backward track does not travel -
`backwardSeedFix` corrects the seed's direction but evidently not the propagation used for material
effects. **That is where the next work goes; it is inside `AtGenfitter`, not in a macro argument.**

### The range constraint is mis-gated and, here, dangerous

Its gates replicated offline on the clusters genfit uses:

| direction | r < 240 mm | z margin | len > 20 | Bragg > 1.3 | all pass |
|---|---|---|---|---|---|
| forward 20-85 deg | 29 % | 100 % | 100 % | 29 % | **9 %** |
| backward 92-140 deg | 73 % | **33 %** | 100 % | **29 %** | **4 %** |

It reaches 4 % of the tracks it was designed for (backward spirals fail on the z margin - many
genuinely leave through the entrance window - and on the Bragg test, whose median ratio is 1.14
against a 1.3 threshold), while 9 % of *forward* tracks pass and are destroyed: the 82 tracks whose
fit changed went from KE_fit/KE_true = 0.993 to **0.017**, because a truncated path is read as a
tiny energy. They then fail the truth match and vanish from the selection, so the summary looked
unchanged. **A constraint that destroys tracks can present as a null result.** Do not enable it for
this channel without re-tuning the gates.

An alternative stopping test I proposed - require the curvature to shrink along the track - was
measured and **fails backwards**: escaping tracks give R_last/R_first = 0.547 against 0.831 for
stopping ones, because a track leaving the chamber has a short final arc the circle fit handles
badly. It measures fit quality, not stopping. Recorded so it is not retried.

### The matEffects A/B, run on identical reconstruction

| arm | θ_lab | n | bias | spread |
|---|---|---|---|---|
| matFX ON (campaign) | 60–85° forward | 726 | −0.3 % | 3.4 % |
| matFX OFF | 60–85° forward | 728 | −0.5 % | 3.4 % |
| matFX ON | 92–110° backward | 137 | **−17.4 %** | 19.1 % |
| matFX OFF | 92–110° backward | 159 | **−10.2 %** | 10.6 % |
| matFX ON | 110–140° backward | 128 | −10.2 % | 15.2 % |
| matFX OFF | 110–140° backward | 132 | −5.9 % | 4.0 % |

Forward tracks are unbiased to 0.3–0.5 % with 3.4 % spread in **both** arms, so the fitter is
sound where (p,p') exercises it. Backward, material effects roughly **double** the bias and
double-to-quadruple the spread — a strong pointer that genfit is integrating energy loss along a
direction the backward track does not travel.

But turning material effects off is a diagnostic, not a fix (it removes the energy-loss correction
altogether), and a residual **−6 to −10 %** bias survives it. So there are two contributions: one
identified in the material model, one still open in the momentum determination itself.

**This qualifies the headline.** The backward σ(Ex) quoted for 2.85/4 T is inflated by a fitter
bias rather than a detector limit, so the *size* of the 7 T gain is not settled — part of it is
7 T avoiding a problem 2.85 T has. The ordering (7 T better) is not in question; the factor is.
