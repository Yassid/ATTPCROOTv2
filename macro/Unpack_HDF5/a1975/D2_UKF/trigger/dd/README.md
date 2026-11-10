# 16C(d,d) elastic scale test: first pass (2026-11-10)

Deuteron-hypothesis genfit (r280, CATIMA) of the D2 runs, `/mnt/f/a1975/gf_dd_r280` (production still running).
First pass on the 11 runs done by 2026-11-10 00:30; run_0023 excluded (IC and scalers truncated), so 10 runs.

    root -l -b -q 'dump_dd.C("run_0016","/mnt/f/a1975/caches/dd_r280/run_0016.root")'      # one per run
    python -I dd_first_pass.py /mnt/f/a1975/caches/dd_r280 92.83 OUT 1000 1350 23

## Inputs
- Luminosity: 316.98 mb^-1 x ic_sca(10 runs)/ic_sca(47) x livetime ratio = 84,192,772/289,561,503 x 0.9746/0.9676
  = 92.83 mb^-1 (carbon, f_beam 0.613). The scalers alone give 96.6.
- Selection: theta_lab 60-89.5, KE 0.3-10 MeV, Ex within +-2.5 sigma of the slice median (sigma 0.23-0.27 MeV),
  background from Ex -4.5..-1.5 (small). No chi2 cut (chi2/ndf ~0.03, the low-energy issue). No PID gate:
  `pid/deuteron_band_d2.json` keeps 4 % of the elastic peak (drawn for higher energies).
- **Beam: 16C selected by IC 1000-1350.** The beam is ~62 % 16C, ~31 % 21O at the same velocity, and d+21O recoils
  fall on the same (theta, KE) locus. IC peaks: C ~1150, O ~2000 (ratio 1.74, Z^2 gives 1.78). Without the gate R
  is ~1.3-1.6.

## Result (no acceptance applied)
R = Y / (L f sigma_ZPL dOmega), f = fraction of the (vz, theta_cm) cell kept by the r_res cut. R = acceptance x scale.
- Plateau R = 0.80-0.92 at theta_cm 16-22, vz 200-900, flat in vz. Potential choice moves it by about 3-4 %.
- R falls above 22 deg (0.5-0.7 at 24-26, 0.2-0.5 at 26-28) and below 16 deg (short tracks, KE < 1.4 MeV).
  The fall above 22 deg is a shape difference with every potential, or a theta_cm offset (sigma falls 25 %/deg).
- **No trigger loss for elastic deuterons**: cells outside r_res >= 28 have the same R as those inside.
- vz 100-200 is low at all angles: elastic recoils start at vz ~100 (none below 75 mm), cause not known.
- Trigger-safe cells only: 0.45-0.86, dominated by the vz 100-200 slice.

## Needed before a scale can be quoted (as of the first pass; see the next section)
- r280 acceptance x efficiency for the elastic (the sim in /mnt/f/a1975_C16_dd_sim is v2 geometry, gain 35000, no
  genfit). If it is ~0.85-0.9 on the plateau, the scale is ~1 within ~10 %.
- The IC-gate efficiency relative to f_beam (5 % of the peak has IC < 1000).
- theta_cm calibration, and the remaining 36 runs.

## Scale with the r280 acceptance (2026-11-10, interim: 2 sim seeds, 10 data runs)
Sim: `macro/Simulation/ATTPC/16C_dt/C16_dd_el_sim_r280.C` through `/mnt/f/a1975/sim_rebuild/r280/run_dd_el.sh`
(gen theta_cm 8-40, Ebeam 184.25 -> reco gain 100000 -> genfit with the run_prod_dd.sh arguments), then
`dump_dd_sim.C` (trees g: generated, d: dump_dd.C columns + truth) and
`python -I dd_scale.py /mnt/f/a1975/caches/dd_r280 92.83 "<dumps>" OUT 23 [NSIG]`.
- Acceptance x efficiency of the sim is ~1 (0.9-1.1 per cell, noise of ~60 generated per cell) at every vz, theta_cm.
- **S = data / (L eps sigma) on the plateau (theta_cm 16-22, vz 250-900): 0.758(12) ZPL; An-Cai .762, Daehnick .754,
  Han-Shi-Shen .797.** Ratio of sums; the error-weighted mean of the cells is biased LOW (a cell whose eps
  fluctuates up gets a low S and a small error: 0.728). Ex window 2.5 / 3.5 / 4.5 sigma moves S by +-2 %.
- theta_cm 12-14 is low at vz < 300: the trigger (deuteron-only events need vz > ~220).
- **Above theta_cm 22 the DATA Ex peak degrades, the sim's does not** (data median -0.12 -> -0.27 MeV and sigma .21 ->
  .30 at 24-26, dissolved, sigma ~1 MeV, at 26+; sim +0.01 +- 0.05 everywhere). The helix there reaches
  r = 2 rho sin(theta_lab) > ~24 cm (KE > 3 MeV): the same class as the (d,t) rmax > 255 KE shifts. Not an
  acceptance; use theta_cm <= 22 only.
- Data resolution is 4x the sim's (Ex sigma .21-.25 vs .05-.06 MeV).

What S is NOT (each measured on the data):
- Spyral validity on the locus (`dump_dd_valid.C`, `dd_valid_frac.py` -> valid_frac.*): data .97-.98 vs sim .97-1.0,
  ratio .98-1.01 for theta_cm 14-28.
- Genfit fit efficiency (`dd_fit_eff.py`, locus from Spyral brho/polar alone; theta_lab = 180 - polar) -> fit_eff.txt:
  data .976-.999, sim 1.000.
- The carbon IC window: 0.1 % of the carbon peak lies at 900-1000 or 1350-1600.
- The IC gate DOES lose the elastic events whose IC pulse is missing: IC < 900 is 4.8 % of all (C+O) locus events,
  i.e. ~5 % of carbon (missing pulses hit both species). Correction x 1/0.95 -> S ~0.80.
- Beam fraction: L uses f_beam 0.613 (measured on H2); the D2 per-track value 0.659 would lower S by 7 %.
Remaining stage: track finding in the real event environment -> `embed_dd.C` (sim deuteron hits overlaid on real
D2 events, chain re-run; same method as embed_vz_dp.C, which gave .68-.91 for the (d,p) protons) and
`embed_dd_ana.py`. Runner `/mnt/f/a1975/caches/dd_r280/embed/run.sh` (runs 16, 31, 79, 102).

## Trace checks (2026-11-10): trigger timing, pile-up, downscaled IC
`trace_pulses.py RUN OUTDIR` (every IC and mesh pulse per event + trigger word) and `trace_checks.py`, runs 16-21, 26, 27,
31 (run 22 left out, see below). Results in `trace_checks.txt/.png`, example traces `trace_ds_noIC_examples_run31.png`.

- **Trigger word** `frib/evt/evt<N>_977`: 1 = IC x mesh coincidence, 2 = downscaled IC, 3 = both (counts match the ic_ds
  scaler). In runs 16-27 the word stored with event i belongs to the waveform of event i+1; in run 31 they pair directly.
  Run 22 fits neither pairing (58839 FRIB events, 39978 reco entries). With the right pairing every coincidence trigger
  has a mesh pulse and 97-98 % of downscaled-IC triggers have none: **the "22 % of triggers without a mesh pulse" are
  the downscaled-IC triggers.**
- **Coincidence window (check 1).** For elastic deuterons the mesh onset follows the drift from the vertex
  (~35 bins earlier than dv 1.104 and ZPP 971.73 predict). No mesh onset is later than ~860 bins (69 us) after the IC
  pulse, while the trace runs to ~917 bins and the full drift is ~1100 bins (88 us). Elastic events whose only mesh
  pulse is the deuteron's disappear below vz ~220 mm; below that they survive only when an early mesh pulse
  (~50 bins after the IC, the residual leaving the beam hole near the pad plane) is present, down to vz ~100. So the
  IC-delay/coincidence explanation is confirmed and quantified: the coincidence covers about the last 69 us of drift.
  The (d,p) adopted result (vz >= 200, r_res >= 28) and the (d,d) plateau (vz 250-900) are outside this loss.
- **Low-IC elastic events (check 2).** 4.8 % of the elastic peak has IC < 1000; 73 % of them have no IC pulse at the
  trigger position, 85 % have a carbon-like pulse (median 1165) elsewhere, before or after at random. Their vz stays
  inside 0-1000 and their Ex is unshifted, so their timing reference is not visibly wrong.
- **Downscaled IC (check 3).** 11.9 % of downscaled-IC triggers have no IC pulse (> 300 ADC) within +-30 bins of the
  trigger position, and 1.2 % have none anywhere in the 164 us trace. Coincidence triggers: 11.3 % and 0.5 %. Some of
  the former have a sub-threshold pulse at the trigger position (23 % above 100 ADC) and a few traces are truncated;
  most have full-size pulses elsewhere, mostly later. A Si pulse at its usual time (26 bins after the IC) is present in
  about half of them. Since the rate is the same for both trigger types, a "true pulse at the trigger time" classifier
  would flag ~11 % of pure-beam triggers too; the downscaled sample measures that false-positive rate.
- Si amplitudes show two groups (~1140 and ~1860), apparently carbon and oxygen.

## H2 (proton-target) trace checks (2026-11-10): `h2_trace_checks.py`
`trace_pulses.py` on runs 106-189 (83 runs; run_0148 has no FRIB waveform group, the known broken run), then
`python -I h2_trace_checks.py /mnt/f/a1975/caches/pd_trig/pulses h2_trace_checks`. Output `h2_trace_checks.{txt,png}`.
- Trigger word pairs DIRECTLY (+0) with the waveform in all 83 runs (downscaled-IC events with no mesh 0.95-0.97).
- Coincidence ceiling: first mesh onset at most 856 bins (68.5 us) after the IC pulse; p99 792 bins. Same as D2.
- IC pulse missing at the trigger position: 8.1 % of coincidence triggers, 8.7 % of downscaled; none anywhere 0.6 % / 1.2 %.
- (p,d) deuterons (adopted selection, joined entry i -> FRIB event i+1, IC correlation 0.83): 99.6 % are word 1, and
  94-98 % have the mesh onset within 150 bins of the IC at every vz (median 32-36 bins): the mesh fires at the
  start of the event, from the residual near the beam hole exit, not from the deuteron drift. The 0-15 mm r_res
  slice is lowest (0.88). These are triggered events only, so this is the timing of the trigger, not its efficiency
  (efficiency vs r_res: project note 2026-11-08).

### Embedding result (2026-11-10 12:05; embed_dd.txt, scale_2seeds_emb.txt)
4 data runs (16, 31, 79, 102) x sim seeds 6001/6002, 2 overlays each: 9678 sim deuterons, reference found alone for
0.912, 17648 overlays.
- **Survival on the plateau (theta_cm 16-22, vz 250-900): 0.794(9) all hosts, 0.814(11) carbon hosts.** The loss is
  the event environment: 0.99 on hosts with < 50 hits, 0.86 at 50-150, 0.74 at 150-300, 0.42-0.50 above 300.
- It falls along the drift: carbon hosts, theta_cm 14-24, vz 100..900 by 100: .905 .929 .880 .821 .859 .785 .769 .695.
  theta_cm: .66 (12-14) rising to .84 (20-24), .74 at 26+.
- `dd_scale.py ... 23 2.5 "/mnt/f/a1975/caches/dd_r280/embed/e_r*.root"` multiplies eps per vz slice by it:
  **S(plateau) = 0.945(16) ZPL** (An-Cai .950, Daehnick .940, Han-Shi-Shen .994); chi2 72 -> 63 / 17.
- With the IC-missing loss (4.8 % of carbon, low-IC events split C:O like the gated ones): **S ~ 0.99**, potential
  spread 0.99-1.04. The D2 luminosity (316.98 scaled, f_beam 0.613) reproduces the 16C(d,d) elastic to a few %.
- **Consequence, not applied:** every D2 channel's acceptance from clean simulation misses this environment loss
  ((d,p) knew it: "22 % embedding loss not applied"; (d,t) acceptance has none). Their cross sections are ~20 % low
  unless corrected. Decision for Yassid. (p,d)/(p,p) partly cancel (same environment, both from clean sims).

### Final with 4 sim seeds (2026-11-10 12:50; scale_4seeds_emb.*): 9991 generated, 8283 selected
- No embedding: S(plateau) 0.765(9) ZPL. **With the embedding survival: 0.954(14) ZPL** (An-Cai .958, Daehnick .949,
  Han-Shi-Shen 1.003). **With the 4.8 % IC-missing loss: 1.002 (ZPL), 0.997-1.054 over the four potentials.**
- The cells scatter beyond their errors (chi2 108/17): read ~10 % per cell as systematic; the plateau mean is not
  affected by it at the stat level quoted.
