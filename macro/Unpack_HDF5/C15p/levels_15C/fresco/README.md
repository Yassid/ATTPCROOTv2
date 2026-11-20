# FRESCO / FRESCOX for 15C(d,d') -- a2091

## The binary

Built here from the maintained LLNL fork (the original fresco.org.uk source links are 404):

    git clone --depth 1 https://github.com/LLNL/Frescox.git
    cd Frescox/source && make MACH=gfortran -j8

gives `frescox` and `sfrescox` (version 7.4-2-gd99affb) in
`/home/yassid/fresco_build/Frescox/source/`. Verified on the shipped `b8ex.nin` test, exit 0.

## The decks

    C15_dd_elastic.in            15C(d,d) elastic, Ebeam 191.25 MeV
    C15_dd_inelastic_0740.in     + the 5/2+ at 0.740, rotational coupling

Run with   `frescox < C15_dd_elastic.in > C15_dd_elastic.out`
Outputs: fort.201 = state 1 (g.s.), fort.202 = state 2 (5/2+), saved here as
`C15_dd_elastic_xsec.dat` and `C15_dd_0740_xsec.dat`.

⚠ **ANGLE CONVENTION.** FRESCO's theta is the DEUTERON centre-of-mass angle. The analysis quotes
theta_CM in the PROJECTILE convention. **Compare at 180 - theta**, or flip one of them before
overlaying. Getting this wrong inverts the distribution.

⚠ **THE OPTICAL POTENTIAL IS NOT FITTED.** The deuteron parameters are starting values from a
global parameterisation, and the deformation length is a placeholder 0.90 fm. Nothing here is a
prediction until the elastic has been searched against the measured distribution with `sfrescox`.
That is the next step, not a detail.

## The measured distributions

    measured_gs.dat      ground state (elastic)
    measured_0740.dat    5/2+ at 0.740

Fitted Gaussian integrals from the adopted configuration (fit g.s. to 10 MeV with no gap,
Gaussian peaks, background free per angle), both codes, only bins that reproduce under a 5 %
change in the sample. **NOT acceptance-corrected** -- the efficiency is flat at 0.93-1.00 inside
the vz window, so it is a small correction, but it has not been applied. Yields are counts, not
cross sections, so only the SHAPE is comparable until a normalisation is established.
