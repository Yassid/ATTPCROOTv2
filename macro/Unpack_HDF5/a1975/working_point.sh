#!/usr/bin/env bash
# a1975 DEFAULT WORKING POINT -- the single place these numbers are allowed to live.
#
# Settled 2026-08-17. Source it; do not copy values out of it:
#
#     source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/working_point.sh
#     echo "$A1975_PD_EBEAM"
#
# WHY THIS FILE EXISTS. The (p,d) beam energy alone was carried as FOUR different defaults --
# 185.0 in fit_pd_ps.C and the saved explorer config, 192.0 in pd/ex_pd_a1975.C, 195.5 in
# pp/open_explorer_pd.sh, and 184.17 inherited from (d,t). At dEx/dEbeam ~ 0.035 MeV/MeV that
# spread is ~0.37 MeV of excitation energy, comparable to the level spacing being resolved. Every
# one of those defaults was individually reasonable at the time it was written, which is exactly
# why documentation alone does not hold: the value has to come from ONE place that scripts read.
#
# THE RULE THAT DECIDES EVERYTHING ELSE: CATIMA WITH MATERIAL EFFECTS ON IS THE DEFAULT.
# It is what fixed the energy scale (on (d,t), 2.490 -> 2.822 against a true 3.103, ~8 sigma) and
# it is now viable at full statistics because the shared-TGeo-navigator poisoning was fixed
# (collapse 58.6% -> 2.5%). But matEffects ON is only safe WITH ITS CO-REQUISITES, listed under
# each channel below. Turning it on without them is what produced every earlier disaster:
#   - wrong geometry  -> 2.5x too much material, and it is silent (density cannot reach the fit
#                        while matEffects is off, so the wrong default survived unnoticed)
#   - matFallback ON  -> a throwing track is refitted with setNoEffects(true) and KEPT, blending
#                        two physics models into one spectrum
#   - manual dE/dx    -> AtGenfitter's SetManualELoss block runs unconditionally, so pairing it
#                        with matEffects counts the vertex gap TWICE (+4.65% -> +9.5% on (d,t))
#   - no dE/dx table  -> genfit applies ZERO stopping power below beta*gamma = 0.05

# ---------------------------------------------------------------------------------------------
# 16C(p,d)15C   runs 0106-0189 (84), H2 at 300 torr
# ---------------------------------------------------------------------------------------------
export A1975_PD_RUNS_LO=106
export A1975_PD_RUNS_HI=189
export A1975_PD_EBEAM=185.0                 # DECIDED 2026-08-17. The value the adopted cross
                                            # sections were fitted at (saved explorer config
                                            # 2026-08-12, fit_pd_ps.C:82). 192.0 and 195.5 are
                                            # stale defaults and are to be purged, not carried.
export A1975_PD_GEO=ATTPC_H300torr_RT       # 3.308e-5 g/cm3. NOT ATTPC_H1bar (8.27e-5, 2.5x).
export A1975_PD_RHO=3.308e-5
                                            # Table names are BARE and the two channels resolve
                                            # them from DIFFERENT directories, which is a trap in
                                            # itself: fitGenfitter_a1975.C:51 prepends
                                            # $VMCWORKDIR/resources/energy_loss/ while
                                            # fitGenfitter_a1975_deuterium.C:90 prepends
                                            # $VMCWORKDIR/macro/Unpack_HDF5/a1975/D2_UKF/.
                                            # A name valid for one channel is missing in the other.
export A1975_PD_ELOSS_TABLE=deuteron_H2_catima.txt
export A1975_PD_PDG=1000010020
export A1975_PD_MASS_AMU=2.0135532
export A1975_PD_Z=1

# genfit co-requisites -- these travel together or not at all
export A1975_PD_MATEFFECTS=kTRUE
export A1975_PD_CATIMA_MSC=kTRUE
export A1975_PD_CATIMA_STRAG=kTRUE
export A1975_PD_MATFALLBACK=kFALSE          # never kTRUE for a material production
export A1975_PD_BACKEXTRAP=kTRUE            # recovers the vertex gap; genfit does it ONCE
export A1975_PD_MANUAL_ELOSS=0              # MUST stay 0 while matEffects is on (double count)
export A1975_PD_BFIELD=-2.85
export A1975_PD_MEASSIGMA=4.0
export A1975_PD_MINITER=2
export A1975_PD_MAXITER=5

# PID gates. The FIT gate is deliberately the LOOSE band and the ANALYSIS gate the tight one:
# genfit gates BEFORE fitting, so whatever the production uses is frozen into the cache. Fitting
# the superset keeps the choice open downstream, where PID is fit-independent and can be re-cut
# without refitting. (Re-cutting needs the PID observables in the cache -- see KNOWN GAPS.)
export A1975_PD_GATE_FIT=pid/deuteron_band.json
export A1975_PD_GATE_ANALYSIS=deuteron_tight.json

# selection
export A1975_PD_IC_TB_LO=1000               # cache_pd_run.C ion-chamber TIME window
export A1975_PD_IC_TB_HI=1350
export A1975_PD_IC_LO=950                   # explorer ion-chamber AMPLITUDE gate
export A1975_PD_IC_HI=1350
export A1975_PD_CHI2NDF_MAX=5
export A1975_PD_KE_LO=0
export A1975_PD_KE_HI=30
export A1975_PD_VZ_LO=0                     # mm
export A1975_PD_VZ_HI=500
export A1975_PD_THETA_CORR_DENOM=2950.8     # theta correction 360/denom about the pivot
export A1975_PD_THETA_CORR_PIVOT=27.0
export A1975_PD_LUMINOSITY=168.3            # mb^-1 = 316.4 * 500/940 (the vz cut keeps that
                                            # fraction of the target); from the (p,p) elastic

# products
export A1975_PD_RECO=/mnt/f/a1975/reco/
export A1975_PD_PROD=/mnt/f/a1975/reco_pd_catima_bx/
export A1975_PD_SUFFIX=_pdcatbx
export A1975_PD_CACHE=/mnt/f/a1975/caches/pd_kin_catima_bx.root

# ---------------------------------------------------------------------------------------------
# 16C(d,t)15C   runs 0016-0103 (47), D2 at 300 torr
# ---------------------------------------------------------------------------------------------
export A1975_DT_RUNS_LO=16
export A1975_DT_RUNS_HI=103
export A1975_DT_EBEAM=184.17                # 11.5 MeV/u x 16
export A1975_DT_DRIFTVEL=1.10424
export A1975_DT_PAR=ATTPC.a1975_deuterium_dv1104.par
export A1975_DT_ZPADPLANE=971.7312          # = DriftLength. NOT 1000; over-determined by dv and
                                            # the TB anchors, so it is forced, not chosen.
export A1975_DT_TBENTRANCE=560
export A1975_DT_GEO=ATTPC_D300torr_v2_geomanager.root
export A1975_DT_RHO=6.61e-5                 # Spyral's config declares 302.115 torr to land here,
                                            # because Spyral builds density with INTEGER A
export A1975_DT_PDG=1000010030
export A1975_DT_MASS_AMU=3.01550072
export A1975_DT_ELOSS_TABLE=triton_D2_300torr.txt
export A1975_DT_MATEFFECTS=kTRUE
export A1975_DT_BACKEXTRAP=kTRUE
export A1975_DT_MANUAL_ELOSS=0
export A1975_DT_MATFALLBACK=kFALSE
export A1975_DT_MATA=2
export A1975_DT_IC_LO=900
export A1975_DT_IC_HI=1400
export A1975_DT_CACHE=/mnt/f/a1975/caches/dt_kin_cateloss.root

# NATIVE CATIMA dE/dx — ADOPTED 2026-08-17 as the default for all genfit productions.
# catima::dedx + dedx_n per RK step with the STEP's own material, replacing the ASCII curve x one
# global density. Measured on the full 47 runs against the otherwise-identical table production:
# collapsed fits 2.51 % -> 1.08 %, usable 35601 -> 36119 (+518), while the 3.103 peak
# (2.872 -> 2.866) and its width (0.247 -> 0.252) are unchanged within errors. The low doublet was
# checked by eye on the explorer and is likewise unchanged.
# SO: it buys STATISTICS, not energy scale. It does NOT address the ~0.24 MeV deficit on 3.103.
export A1975_DT_CATIMA_ELOSS=kTRUE
export A1975_DT_CATIMA_ELOSS_FULL=kFALSE   # replaces Bethe-Bloch above bg=0.05 as well; on 8 runs
                                           # it changed nothing (2.918 vs 2.919). A full 47-run
                                           # test is running — update this line with the result.
# The ASCII table is STILL LOADED as a per-step fallback, so A1975_DT_ELOSS_TABLE above is still
# required on disk. CATIMA answers in normal operation and the table is never consulted; dropping
# it is a one-line change once the full-range test settles.

# THE (d,t) PHYSICS RESULT DOES NOT COME FROM GENFIT. The adopted excitation energies are the
# Spyral InterpSolver's (5 levels, constant +0.261 MeV offset, the 4.657 predicted to 27 keV).
# None of the three genfit/UKF energies resolves the two low 15C levels. "CATIMA is the default"
# is a statement about how the genfit productions are configured, NOT a decision to replace the
# solver as the source of the (d,t) numbers.

# ---------------------------------------------------------------------------------------------
# PURGE LIST -- values that exist in the tree and are WRONG. Delete on sight; do not "support".
# ---------------------------------------------------------------------------------------------
#   Ebeam 192.0            pd/ex_pd_a1975.C:50                  -> 185.0
#   Ebeam 195.5            pp/open_explorer_pd.sh               -> 185.0
#   geoName ATTPC_H1bar    fitGenfitter_a1975.C default         -> ATTPC_H300torr_RT
#   gasDensity 9.0e-5      fitUKF_a1975_deuterium.C default     -> 6.61e-5 (9.0e-5 is ~408 torr)
#   chi2max 1e9            every explorer caller, meaning "no cut" -- 1e9 is ALSO the collapsed-fit
#                          sentinel, so `chi2ndf > 1e9` is false and collapsed fits sail through.
#                          Use 1e8. mkexp_pp now drops them unconditionally regardless.
#   [ -s "$f" ]            as a "job finished" test -- a truncated file is non-empty. Use
#                          pd/root_ok.C (opens, closed cleanly, has a TTree).
#
# KNOWN GAPS, open rather than decided:
#   - cache_pd_run.C stores no PID observables (no brho/dedx), so the tight gate CANNOT yet be
#     re-cut at analysis time as the gate policy above intends. Adding them is the enforcement
#     step for that policy.
#   - (p,p) has NOT been re-run under the CATIMA default; its numbers are still matFX-off.
#   - AtGenfitter.cxx:593 runs the manual-eloss block unconditionally. Nothing but convention
#     stops a future caller from re-creating the double count; a guard needs a library rebuild.
