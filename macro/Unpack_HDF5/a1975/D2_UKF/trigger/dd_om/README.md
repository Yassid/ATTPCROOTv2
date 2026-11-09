# 16C(d,d) elastic reference curves

FRESCO decks and dsigma/dOmega for d+16C at Ed 21.9, 22.4 and 23.1 MeV (vz 900, 550 and 0 of the 16C beam).
Potentials from PtolemyCpp's library: A An-Cai, H Han-Shi-Shen, Z Zhang-Pang-Lou, D Daehnick (outside its
A >= 27 range, comparison only). `spread.txt` is produced by `compare.py`.

Conventions tested against PtolemyCpp on the same potential:
- Surface depth is the peak depth. FRESCO with W_s absorbs 1090.4 mb, PtolemyCpp 1092.9 mb; 4 W_s gives 1368 mb.
- Spin-orbit: FRESCO needs Vso/2 of the quoted value for the deuteron. Without it the curves differ by up to
  36 % at 40-50 deg; with it, or with Vso = 0, they agree to 0.4 %.

Use for a scale test: spread of the three in-range potentials is below 7 % only at theta_cm 16-24 deg. Above
that it is 12-90 %. The energy dependence is steep (about 4-15 % between 21.9 and 23.1 MeV in that window), so
the comparison has to be made per vz slice.
