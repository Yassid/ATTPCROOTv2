# 17C(d,p)18C — SOLARIS + AT-TPC proposal simulation

Built 2026-09-02 for the 17C M_n/M_p proposal (`C17p_FRIB_Proposal`). That proposal's primary
measurement is 17C(p,p') and 17C(d,d'); it notes that the deuterium day also delivers (d,p), (p,d),
(d,3He) "for free". This is the (d,p) arm.

**Angular distributions are FLAT** (isotropic in the CM — `AtTPC2Body` samples flat in
cos θ_cm, `AtTPC2Body.cxx:169`). No DWBA input yet. Acceptance is a per-bin ratio and resolution
is measured per slice, so both are independent of the true shape: when the DWBA calculations
arrive they *weight* these results, they do not require regenerating anything.

## Configuration

| | |
|---|---|
| beam | 17C at 8.37 MeV/u = 142.29 MeV, ReA6, 940 pps |
| gas | D2, 300 torr, 293 K (`ATTPC_D300torr_v2`, ρ = 6.5643e-5 g/cm³) |
| field | 2.85 T (nominal SOLARIS + AT-TPC) |
| pads | the real AT-TPC pad plane, drift 1000 mm |
| par | `ATTPC.C17dp_D300torr_b285.par` |
| products | `/mnt/f/C17dp/`, 4 levels × 24000 events |

Everything except the beam and the reaction is held at the 14C(d,p) reference campaign's
configuration (`../14C_dp/RESULTS.md`), so a difference between the two channels is the reaction
and not a setting.

**Q(17C(d,p)18C) = +1.9594 MeV.** 18C has S_n = 4.184 MeV and, per ENSDF, exactly four bound
states — so this is the whole bound spectrum, with no unbound padding:

| | E_x [MeV] | Jπ | |
|---|---|---|---|
| `gs` | 0 | 0⁺ | |
| `ex1588` | 1.588 | 2⁺ | 15.5 ps, B(E2) = 0.000364 e²b² — the state the proposal discusses |
| `ex2515` | 2.515 | (2⁺) | < 3.2 ps |
| `ex3972` | 3.972 | (2,3)⁺ | |

## Verified from MC truth, not assumed

`check_beam_C17.C`: beam momentum exact (2.129007 GeV/c), KE 142.29 MeV = 8.370 MeV/u,
Q = +1.9594 MeV, two-body kinematics close to machine precision, vertices uniform in z.

Independently, `dp_summary_C17.C` solves for the beam energy that makes the *generated* E_x close,
per event, and fits it against z:

    E_beam(z) = 141.51 - 0.014504 z[mm]   ->  14.50 MeV lost over the metre  (n = 8694)
    CATIMA, same gas, independently        ->  14.80 MeV

2 % agreement between the transported truth and the stopping-power model. The constant used for
the constant-E_beam analysis is **135.0 MeV**, the mean over a uniform vertex distribution.

## The kinematics land on top of 14C(d,p) — the reference campaign transfers

Q flips sign against 14C(d,p) (+1.96 vs −1.01 MeV) and the beam is 3 MeV/u slower. **The two
effects nearly cancel**, so the proton kinematics are almost the same, and the transfer peak is
backward-lab and slow in both:

| θ_cm | 17C(d,p) θ_lab / KE | 14C(d,p) θ_lab / KE |
|---|---|---|
| 5° | 166.0° / 2.01 MeV | 163.4° / 1.68 MeV |
| 20° | 130.5° / 3.13 MeV | 124.6° / 3.11 MeV |
| 30° | 113.6° / 4.60 MeV | 107.9° / 4.96 MeV |
| 60° | 80.2° / 11.90 MeV | 76.5° / 14.19 MeV |

`plots/kinematics_C17dp.png` draws this: θ_lab vs θ_cm, the proton KE locus with the **MC truth
overlaid on the analytic curves** (11328 protons, they coincide — which is what validates the
generator, the transport and these formulae against each other), the level separation in proton KE,
and the |dE_x/dKE| leverage. The transfer-peak band is shaded on every panel, because that is where
every quotable number comes from and it is *not* where a flat generator puts most of its events.

**θ_cm CONVENTION.** The standard (d,p) angle, from the deuteron direction in the entrance
channel — the *supplement* of the proton's cm polar angle about the beam. This is what
`acceptance_C14.C:44` uses (`θ_cm = π − acos(...)`), so small θ_cm is where a stripping
distribution has its yield. `dp_kinematics_C17.C` initially used the complement and every table
came out reversed; it is now checked against the 14C campaign, reproducing its documented
"θ_cm 60 → 76°, 14.3 MeV" and "θ_cm 20 → 125°, 3.1 MeV".

## Acceptance

Overall 0.746 / 0.744 / 0.724 / 0.704 for the four levels (14C(d,p) gave 0.71–0.73). Flat at
~0.80 over θ_cm 50–130, falling to zero past 165°.

**OPEN ITEM — a deep, level-dependent acceptance hole at θ_cm 25–45.** It is mild for the g.s.
(0.686 at 45–50) and severe for the 3.972:

    theta_cm    20-25  25-30  30-35  35-40  40-45  45-50
    gs           0.897  0.883  0.849  0.908  0.835  0.686
    ex3972       0.779  0.447  0.437  0.474  0.156  0.749

θ_cm 40–45 is θ_lab ≈ 95–100°, i.e. a slow proton emitted nearly perpendicular to the beam, which
advances almost no z per turn. The 14C(d,p) campaign found the same feature at every field and
left the cause open; here there is only one field, so it is not further resolved. It matters: it
sits inside the transfer-peak range for the upper levels, and a 0.16 acceptance is a ×6
correction. **Do not quote a cross section in those bins without understanding this.**

## Resolution — the headline

σ(E_x) as IQR/1.349 (robust: the fit-failure tail is heavy enough to fake large effects).
"floor" = the same inversion on the **truth** kinematics at the same constant beam energy, i.e.
what the *method* costs with a perfect detector.

**At the transfer peak (θ_cm 2–40°, θ_lab ≈ 100–175°):**

| level | n | σ const-E_beam | σ vertex-corrected | floor | tail(>1 MeV) |
|---|---|---|---|---|---|
| 0⁺ g.s. | 1180 | 0.168 | **0.052** | 0.167 | 1.7 % |
| 2⁺ 1.588 | 1091 | 0.194 | **0.058** | 0.194 | 3.7 % |
| (2⁺) 2.515 | 961 | 0.182 | **0.055** | 0.177 | 5.9 % |
| (2,3)⁺ 3.972 | 786 | 0.235 | **0.066** | 0.220 | 8.8 % |

**Measured equals the floor to within a few percent at the transfer peak.** The detector is
invisible there — the entire limitation is the constant-E_beam assumption — so no detector change
(field, pad pitch) can improve it, and the *software* correction takes σ(E_x) to **52–66 keV**.
This reproduces the 14C(d,p) result (0.178 → 0.064 MeV) on an independent channel.

σ(E_x) vs θ_lab, g.s., measured / floor: 1.839/0.464 (0–15°), 1.000/0.592 (30–45°),
0.365/0.291 (75–90°), 0.165/0.172 (105–120°). Forward of ~60° the detector dominates and the
resolution is 0.7–1.8 MeV; backward of ~90° measured and floor coincide. The 150–165° bin has
22–58 tracks per level and its widths should not be quoted.

## Do the 18C levels separate? Yes.

Δ(E_x) / (σ₁ + σ₂), vertex-corrected:

| pair | Δ [MeV] | at the transfer peak | all accepted |
|---|---|---|---|
| g.s. → 1.588 | 1.588 | 14.4 | 1.93 |
| **1.588 → 2.515** | **0.927** | **8.2** | 1.29 |
| 2.515 → 3.972 | 1.457 | 12.0 | 2.30 |

The tightest pair in 18C — 927 keV — is separated by 8.2 σ where the yield is. The proposal quotes
~300 keV as the achievable AT-TPC resolution; at the transfer peak this simulation gives ~60 keV
with the vertex correction and ~200 keV without it, so the proposal's claim is conservative.

The "all accepted" column is much worse only because a flat generator puts most events at large
θ_cm / forward lab angles, where resolution is 1–2 MeV. A real stripping distribution is
forward-peaked in θ_cm and therefore sits in the good region; this column is the pessimistic bound,
not the expectation.

## Files

| | |
|---|---|
| `C17_dp_sim.C` | the simulation |
| `make_c17dp_par.sh` | writes `ATTPC.C17dp_D300torr_b285.par`, verifying its inherited values |
| `accumulate_C17dp.sh` | one sample: generate → reco → genfit → acceptance → Ex resolution |
| `run_C17dp_campaign.sh` | the four levels (`DP_ONE_WAVE=1` to merge waves on a resume) |
| `check_beam_C17.C` | truth verification of beam, vertices and two-body closure |
| `dp_kinematics_C17.C` | the analytic map + `plots/kinematics_C17dp.png`; run it on 14C parameters as a control |
| `dp_summary_C17.C` | full-range resolution, floor, vertex correction, level separation |
| `dp_plots_C17.C` | `plots/C17dp_summary.png` — acceptance, E_x spectra, resolution |

Reco, genfit, acceptance and Ex resolution reuse the 14C macros unchanged — they already take the
target, ejectile, residual and beam masses as arguments.

## Traps hit while building this, worth not repeating

1. **`ex_res_C14_hf.C` hard-codes θ_lab slices 20–90°.** Correct for (p,p'), where 90° *is* the
   physical limit; blind to this channel's entire transfer peak. Its `ALL` row is complete and its
   flat `res` TTree holds every accepted event, so `dp_summary_C17.C` re-bins 0–180° with nothing
   re-run.
2. **`exTrue` in that tree is not the generated E_x** — it is the truth kinematics through the same
   inversion at the same constant beam energy (`ex_res_C14_hf.C:162`). Measuring a residual or
   solving for E_beam(z) against it is self-fulfilling: the first version of `dp_summary_C17.C`
   got a floor of exactly 0.000 and an E_beam(z) slope of exactly 0. **Reference the *generated*
   E_x.** `exTrue` is then free — it already *is* the method-floor residual.
3. **Two campaign drivers on one output tree.** A driver launched with `nohup … &` inside a
   backgrounded call had its children killed but survived; a second was then started. They wrote
   the same genfit files concurrently and each read the other's half-written output
   (`probably not closed / missing cbmsim`), losing two hours. `run_C17dp_campaign.sh` now takes a
   PID-checked driver lock.
4. **A stage is complete only when its product reads back.** The genfit resume check used to be
   "file non-empty AND log mentions CATIMA" — both true the moment genfit *starts* writing. It now
   requires `cbmsim` with the full entry count.
5. **`|| true` on `root -b -q -l -e` is load-bearing.** ROOT exits non-zero (8) even when it printed
   the answer; under `set -eo pipefail` that killed all four samples silently, straight after the
   check, with no error line — having just spent two hours on reco.
