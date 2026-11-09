# 16C(d,d) elastic reference curves

FRESCO decks and dsigma/dOmega for d+16C at Ed 21.9, 22.4 and 23.1 MeV (vz 900, 550 and 0 of the 16C beam).
Potentials from PtolemyCpp's library: A An-Cai, H Han-Shi-Shen, Z Zhang-Pang-Lou, D Daehnick (outside its
A >= 27 range, comparison only). `spread.txt` is produced by `compare.py`.

Conventions tested against PtolemyCpp on the same potential:
- Surface depth is the peak depth. FRESCO with W_s absorbs 1090.4 mb, PtolemyCpp 1092.9 mb; 4 W_s gives 1368 mb.
- Spin-orbit: FRESCO needs Vso/2 of the quoted value for the deuteron. Without it the curves differ by up to
  36 % at 40-50 deg; with it, or with Vso = 0, they agree to 0.4 %.

Use for a scale test: spread of the three in-range potentials is 3-6 % at theta_cm 18-24 deg and 8 % at 16 deg.
At 26 deg and above it is 10-86 %. The energy dependence is steep (about 4-15 % between 21.9 and 23.1 MeV in that window), so
the comparison has to be made per vz slice.

The usable angles also have to pass the fiducial cut r_res >= 28 mm (`fiducial_window.txt`): theta_cm 18 deg needs
vz <= 193 mm (below the 200 mm window start), so only about 20-24 deg survives, at vz 200-390 mm. The cross section
falls by a factor 3 from 20 to 24 deg, so the result will be sensitive to the theta_cm calibration.
