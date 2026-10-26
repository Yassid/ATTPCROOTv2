Dec 2014 AT-TPC  --  4He + 4He elastic cross section
====================================================

You start from the UNPACKED ROOT FILES and do the tracking yourself.

QUICK START
  1. edit config.py -> set DATA_ROOT to the data folder on this machine
     (or:  export DEC2014_DATA=/path/to/dec2014_data )
  2. source <your ROOT>/bin/thisroot.sh        # needed for step 1 only
  3. python3 00_check.py                       # verifies everything, fix what it reports
  4. ./run_all.sh                              # the whole chain

REQUIREMENTS
  ROOT            any recent version, for step 1 only.  ATTPCROOT is NOT needed:
                  the files carry their own StreamerInfo, so plain ROOT reads them.
  python3         numpy scipy matplotlib scikit-learn==1.3.2
                  pip3 install --user numpy scipy matplotlib scikit-learn==1.3.2

THE CHAIN
  00_check.py          verify packages, ROOT, paths, data
  01_dump.sh           ROOT files -> hit CSVs            (needs ROOT)
  02_track.py          TRACKING, option A: standalone clustering (needs scikit-learn)
  02b_track_ransac.sh  TRACKING, option B: ATTPCROOT's own AtSampleConsensus (needs
                       the ATTPCROOT environment).  Writes the same pickle format, so
                       steps 3-6 take either without changes.
  03_kinematics.py     angles and energies, no cuts
  04_selection.py      the elastic constraint, and the energy-scale check
  05_acceptance.py     geometric acceptance eps(E_cm, theta_cm)
  06_cross_section.py  relative dsigma/dOmega at theta_cm = 90 deg

  lib/dump_hits.C      the ROOT macro step 1 runs
  lib/clustering.py    the tracking algorithm
  lib/geometry.py      rebuild an event in 3D, fit lines, find the vertex
  lib/eloss.py         range <-> energy
  lib/loader.py        loop over runs, cache the reconstructed table
  data/                range-energy table
  out/                 plots and intermediate files (created on first run)

DISK USED
  CSVs    ~16 GB   (step 1; delete them once step 2 has run)
  tracks  ~1.9 GB  (step 2; this is what steps 3-6 read)

ON TRACKING -- TWO FINDERS, AND THEY DISAGREE
  Option A, 02_track.py: standalone clustering.  HDBSCAN in [position, orientation]
  space, because plain density clustering merges the two arms exactly at the vertex.

  Option B, 02b_track_ransac.sh: ATTPCROOT's SampleConsensus::AtSampleConsensus, the
  method of Ayyad et al., NIM A 880 (2018) 166.  It was missing from this branch --
  AtRansac/AtPRA sit behind if(PCL_FOUND) and PCL is not installed -- and was ported
  across from develop (commit a2abf613), which needed no PCL.

  ON THE SAME 200 EVENTS OF RUN 38 THEY GIVE DIFFERENT ANSWERS:
      standalone clustering   110/200 events with 3 segments
      AtSampleConsensus        63/200   (iterations 1000, minHits 10, dist 6 mm)
  Loosening the RANSAC settings raises it from 44 to about 70 but it saturates there.
  Of the events each keeps, a similar fraction pass the elastic constraint (13% vs 16%),
  so the difference is in FINDING tracks, not in their quality.

  Which is right is not established.  A plausible cause is that RANSAC fits the dominant
  line first and removes its inliers, so the beam -- which has the most hits -- can take
  hits from the arms near the vertex and leave them below the minimum.  Resolving this is
  real work and it matters: the tracking efficiency is what the cross section is divided
  by.  Run both, compare, and do not assume either is the truth.

THE ONE NUMBER THAT IS NOT SETTLED
  V_DRIFT.  Three methods give 2.142, 2.237, 2.278 cm/us.  Run the chain at all three
  and quote the spread as a systematic:
     DEC2014_VD=2.142 ./run_all.sh
  Each caches separately in out/, so nothing is overwritten.

READ THE FULL GUIDE BEFORE CHANGING ANY CONSTANT IN config.py.
