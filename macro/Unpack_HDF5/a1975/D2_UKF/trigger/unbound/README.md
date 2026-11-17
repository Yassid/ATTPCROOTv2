# 16C(d,p)17C unbound states on the adopted 2026-11-16 products (NOT ADOPTED)

Live copies run from `/mnt/f/a1975/sim_rebuild/unb_hdbtau/` (outputs there: `unb_vzscan.csv`, `unb_ad_2026-11-16.png`).

- `unb_ad.py` first pass (3.9 MeV sim efficiency for all states; superseded by unb_ad2.py at the theta_cm edges).
- `unb_ad2.py val` validates the state-specific efficiency on the g.s. sim; `unb_ad2.py scan` gives dsigma/dOmega for the
  three states (2.59, 3.44, 4.01 MeV calibrated) in theta_cm 22-38 for vz >= 100 ... 700 mm. **No trigger term.**
- `unb_lfit.py` L = 0-3 one-normalisation fits (FRESCO templates of ~/C17_dp_fits, post coupling fix).
- `fig_unb.py` the figure.
- Truth for the Ex 3.9 sim: `macro/Simulation/ATTPC/16C_dp/dump_truth_dp.C` (reproduces trig/keys.txt for the g.s. sim).

## 2026-11-17: Breit-Wigner fits (current)
- `unb_ad3.py` gaussian model with the acceptance in the fit, common shift per bin, Ex 1.3-5 (positions do not drift with angle; the 3-gaussian
  description is ambiguous, the third peak wanders 3.8-4.5).
- `unb_bw.py` BW states of the C17_dp_fits table (2.763 1/2-, 2.980 3/2+, 3.661, 4.231 3/2+, 4.841 1/2+, 5.91 3/2+; Gamma fixed, penetrabilities),
  shared offset/sigma/tilt by Nelder-Mead, per-bin linear Poisson fits. Free offset -0.21.
- `unb_bw_force.py` / `unb_bw_force32.py` 2.980 forced in (dropped 2.763, or tied at the YSOX/WBP C2S x ADWA ratio), offset free / -0.21 / -0.32.
  -0.32 = minus the bound-blend centroid (5/2+ carries 93-96 % of the bound peak for YSOX, WBP and PLB 811 C2S).
- `unb_ad32.py` ADs, L=0-3 fits and shell-model S at offset -0.32 -> unb_ad32.csv; `fig_*.py` the figures. `unb_sm.py` S vs shell model.
- Shell-model inputs (YSOX, WBP C2S; from Yassid) are NOT in the repo: /mnt/f/a1975/sim_rebuild/unb_hdbtau/shellmodel/.
- 2026-11-17: TBSTEP=2 SUF=_2deg gives 2-deg bins (outputs *_2deg*); fig_ad32b.py sums 2.763+2.980 with the Fisher covariance (split unstable at 2 deg).
- 2026-11-17: THI=137 (theta_lab window; efficiency falls back to the g.s. sim where the 3.9 sim is empty, validated to theta_cm 14) and KEMIN=0.8
  (sim: KE 0.8-1.2 fits 0.87-0.92, bias -0.01, Ex res as above 1.2; below 0.8 degrades). Run: THI=137 KEMIN=0.8 TBLO=10 TBSTEP=4 SUF=_137ke08.
