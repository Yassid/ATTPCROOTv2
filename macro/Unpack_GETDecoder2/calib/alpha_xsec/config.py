"""
EDIT THIS FILE FIRST.  Everything else reads its settings from here.

The only line you normally need to change is DATA_ROOT, which must point at the
folder on the hard drive that contains  nocut/P300_B0/cands_XXXX.pkl
"""
import os

# ---------------------------------------------------------------- where the data is
# Point this at the top of the data folder on your hard drive.
# On Linux a USB drive is usually under /media/<user>/<label>/ ; on macOS /Volumes/<label>/ .
DATA_ROOT = os.environ.get("DEC2014_DATA", "/media/student/alphas/dec2014_data")

ROOT_DIR   = os.path.join(DATA_ROOT, "P300_B0")             # the unpacked ROOT files
CSV_DIR    = os.path.join(DATA_ROOT, "csv")                 # step 1 writes hit CSVs here
TRACKS_DIR = os.path.join(DATA_ROOT, "tracks")              # step 2 writes tracks here
CANDS_DIR  = TRACKS_DIR                                     # what steps 3+ read
ELOSS     = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "data", "eloss_alpha_heco2_11p4.txt")
OUT_DIR   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# ---------------------------------------------------------------- beam and gas
# These are SETTLED.  Do not change them without reading the guide's status table.
E_BEAM   = 11.40     # MeV, 4He entering the gas.  NOT the 1.95 MeV/u in the Excel run log.
PRESSURE = float(os.environ.get("DEC2014_P", 299.5))   # torr; DEC2014_P overrides for the density study
P_TABLE  = 150.0     # the range table is tabulated at this pressure; ranges scale as P/P_TABLE
TILT     = 6.47      # degrees, detector tilt relative to the beam

# ---------------------------------------------------------------- detector geometry
R_PAD        = 275.0   # mm, outer radius of the pad plane
DRIFT_LENGTH = 1000.0  # mm, window (z=1000) to pad plane (z=0)
TB_NS        = 160.0   # ns per time bucket
T0_BUCKETS   = -34.0   # drift-time origin:  z = V_DRIFT * (tb - T0_BUCKETS) * TB_NS * 1e-3 * 10

# ---------------------------------------------------------------- the one unsettled number
# Three methods disagree by 6%.  Run the analysis at all three and quote the spread.
#   2.142  from the beam's Bragg peak
#   2.237  from the arm-arm opening angle
#   2.278  from the theta1+theta2 = 90 deg constraint   <- currently preferred
# THE DRIFT VELOCITY IS SETTLED at 2.265 cm/us.  The "6% spread between three methods"
# that this analysis carried as its dominant systematic was NOT real:
#
#   * The value 2.142 came from assigning 145 time buckets to 497 mm, taken as 90% of the
#     range of a 7.80 MeV alpha.  That beam energy was later corrected to 11.40 MeV, whose
#     range (1071 mm) EXCEEDS the 1006 mm the tilted chamber offers -- so the beam never
#     stops inside and there is no Bragg peak there to measure.  The determination was
#     invalid, not merely imprecise.
#   * 2.237 and 2.278 were the same kinematic constraint read off coarsely.
#
# So there is ONE valid determination, the theta1 + theta2 = 90 degree constraint, and it is
# cleaner than it looks: the sum is insensitive to the beam direction to first order,
# because rotating the beam within the reaction plane raises one arm angle and lowers the
# other.  It therefore measures v_D alone and does NOT trade off against the tilt.
#
# Measured on all 68 runs, peak centroid against v_D:
#     2.245 -> 90.492     2.255 -> 90.244     2.265 -> 89.995     2.275 -> 89.748
# The +-1% band below moves the peak by about +-0.25 degrees, against a 4 degree FWHM.
V_DRIFT      = float(os.environ.get("DEC2014_VD", 2.265))
V_DRIFT_SET  = (2.245, 2.265, 2.285)   # the +-1% band, see below

# ---------------------------------------------------------------- analysis cuts
SUM_WINDOW   = 6.0    # deg; keep events with |theta1 + theta2 - 90| < this
MIN_ARM_MM   = 0.0    # NO arm-length cut.  The clustering already imposes a ~5 mm floor.
