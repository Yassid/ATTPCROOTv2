# 16C(d,p)17C unbound states on the adopted 2026-11-16 products (NOT ADOPTED)

Live copies run from `/mnt/f/a1975/sim_rebuild/unb_hdbtau/` (outputs there: `unb_vzscan.csv`, `unb_ad_2026-11-16.png`).

- `unb_ad.py` first pass (3.9 MeV sim efficiency for all states; superseded by unb_ad2.py at the theta_cm edges).
- `unb_ad2.py val` validates the state-specific efficiency on the g.s. sim; `unb_ad2.py scan` gives dsigma/dOmega for the
  three states (2.59, 3.44, 4.01 MeV calibrated) in theta_cm 22-38 for vz >= 100 ... 700 mm. **No trigger term.**
- `unb_lfit.py` L = 0-3 one-normalisation fits (FRESCO templates of ~/C17_dp_fits, post coupling fix).
- `fig_unb.py` the figure.
- Truth for the Ex 3.9 sim: `macro/Simulation/ATTPC/16C_dp/dump_truth_dp.C` (reproduces trig/keys.txt for the g.s. sim).
