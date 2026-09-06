#!/usr/bin/env bash
# Browser explorer for 16C(d,p)17C at dv 1.10424 — the CATIMA/genfit production of 2026-09-05.
# Sibling of open_explorer_dt_dv1104.sh; same pipeline, same generator.
#
#   ./open_explorer_dp.sh [Ebeam] [icMin] [icMax]
#
# THE IC GATE IS APPLIED HERE, unlike the dt_dv1104 launcher which opens it to 0..1e9. The page
# has no notion of an ion chamber, so a gate not applied in mkexp_pp cannot be applied at all,
# and ungated these caches carry a large beam contaminant. 900-1300 is the window every other
# (d,p) number in this analysis was measured with; pass 0 1e9 to see it ungated.
#
# BEWARE run_0023 and run_0080: their IC files are TRUNCATED (360 entries for 18862 events, and
# 8579 for 13120), so 98% and 35% of their tracks carry ic = -1 and are dropped by this gate.
# That is upstream of this channel — the adopted (d,t) cache has the same damage.
#
# Masses are the NUCLEAR light / ATOMIC heavy convention used throughout this channel; mixing
# them is 0.511 MeV per electron. Checked: Q = -1.498 vs Sn(17C) - B(d) = -1.496.
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.25}"
ICMIN="${2:-900}"
ICMAX="${3:-1300}"
# Opening CUT window written into the page (all still movable in the browser). The point of
# making these arguments is that the interesting (d,p) regions are not the template's defaults:
#   ./open_explorer_dp.sh 184.25 900 1300 90 180 1 backward
# gives theta_lab 90-180 and KE < 1 MeV -- the sub-genfit-threshold backward region.
# NOTE that window CANNOT contain the 17C bound states: at 170 deg the g.s. proton is 1.51 MeV,
# above a 1 MeV ceiling, and the Ex band such a cut shows (2.0-9.3 MeV) is the cut's own
# kinematic window, not structure. Verified by direct kinematics 2026-09-05.
THCUTLO="${4:-0}"
THCUTHI="${5:-180}"
KECUTHI="${6:-1000}"
TAG="${7:-}"
CACHE="${CACHE:-/mnt/f/a1975/caches/dp_kin_findscan.root}"
TMP="${DPTMP:-/mnt/f/a1975/caches/.explorer}"
mkdir -p "$TMP"
OUT="$HOME/a1975_C16_dp_dv1104${TAG:+_$TAG}_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$CACHE" ]] || { echo "ERROR: no cache at $CACHE — run cache_dp_dv1104.sh first"; exit 1; }
[[ -f "$MK" ]]    || { echo "ERROR: generator missing at $MK"; exit 1; }
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

root -b -l -q "$HERE/mkexp_pp.C(\"$CACHE\",\"$TMP/exp_dpfindscan.root\",1e9,$ICMIN,$ICMAX)"
# 17C levels: the three BOUND states (unresolved here at FWHM ~0.85), Sn, then the resonances
# from the C17_dp_fits analysis.
root -b -l -q "$MK(\"$TMP/exp_dpfindscan.root\",\"$OUT\",\"16C(d,p)17C  dv 1.10424  CATIMA\",$EBEAM,16.0147013,2.0135532,1.00727646688,17.0225864,16,\"0:g.s.,0.217:1/2+,0.331:5/2+,0.729:Sn,2.763,3.661,4.231,4.841\",\"\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
THCUTLO="$THCUTLO" THCUTHI="$THCUTHI" KECUTHI="$KECUTHI" python3 - "$OUT" <<'PY'
import sys
p=sys.argv[1]; s=open(p,encoding='utf-8').read()
# --- MAP AXES. The template defaults are the FORWARD-ONLY ones inherited from (p,p')/(d,t):
# theta_lab display 0..95 and KE 0..40. In (d,p) the proton covers the whole angular range and
# the transfer strength is BACKWARD, so 95 hides the physics outright; and the g.s. proton is
# 48.6 MeV at 20 deg, so a 40 MeV ceiling clips the forward end of the locus off the panel.
# Data here spans theta 10-170 (the production window) with a median of 73.5 deg.
# These are DISPLAY ranges only -- they never remove events, that is what the cuts group does.
import os
_cuts = (("thLo","1","0",os.environ.get("THCUTLO","0"),"theta_lab cut lo"),
         ("thHi","1","180",os.environ.get("THCUTHI","180"),"theta_lab cut hi"),
         ("keHi","0.5","1000",os.environ.get("KECUTHI","1000"),"KE cut hi"))
for _id,_st,_o,_n,_what in _cuts:
    if _o == _n: continue
    _a='id="%s" step="%s" value="%s"'%(_id,_st,_o); _b='id="%s" step="%s" value="%s"'%(_id,_st,_n)
    if _a not in s:
        raise SystemExit("open_explorer_dp.sh: %s default is not %s any more"%(_id,_o))
    s=s.replace(_a,_b,1); print("  %-22s %s -> %s"%(_what,_o,_n))
for _id, _old, _new, _what in (("thMax","95","180","theta_lab display max"),
                               ("keMax","40","60","KE display max")):
    _a = 'id="%s" step="%s" value="%s"' % (_id, "5" if _id=="thMax" else "1", _old)
    _b = 'id="%s" step="%s" value="%s"' % (_id, "5" if _id=="thMax" else "1", _new)
    if _a not in s:
        raise SystemExit("open_explorer_dp.sh: %s default is not %s any more -- template changed, "
                         "re-check the axis patch" % (_id, _old))
    s = s.replace(_a, _b, 1)
    print("  %-22s %s -> %s" % (_what, _old, _new))
s=s.replace("const SETS = RAW.ke ? {ukf: RAW} : RAW;",
            "const SETS = RAW.ke ? {ukf: RAW} : RAW;\nconst LBL={ukf:'genfit'};\nconst lbl=f=>LBL[f]||f;")
s=s.replace("b.textContent = f;","b.textContent = lbl(f);")
s=s.replace("${ACTIVE.toUpperCase()}","${lbl(ACTIVE).toUpperCase()}")
open(p,'w',encoding='utf-8').write(s)
PY
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to Desktop"
echo "explorer -> $OUT"
