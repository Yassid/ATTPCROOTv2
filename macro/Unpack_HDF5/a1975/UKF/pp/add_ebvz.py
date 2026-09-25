#!/usr/bin/env python3
"""Add a VERTEX-DEPENDENT BEAM ENERGY control to a generated explorer page.

WHY. The page computes Ex with ONE beam energy for every track. The beam loses energy traversing
the gas, so a reaction deep in the chamber happened at a lower energy than one near the entrance,
and a constant Ebeam therefore returns an Ex that is too high -- by more at the far end.

Measured on the 16C(d,p)17C GROUND-STATE simulation, TRUTH ONLY, no reconstruction anywhere in it
(macro/Simulation/ATTPC/16C_dp, sim_res_dp.C): generated Ex = 0 came back as

    <vz> mm      71    200    297    398    498    600    698    800    877
    Ex  MeV   +0.068 +0.134 +0.189 +0.243 +0.296 +0.356 +0.406 +0.461 +0.514

a +0.45 MeV drift across the chamber that is NOT reconstruction. Scanning the beam energy that
zeroes it gives a straight line to ~3 %:  11.0 MeV per metre, ~8.9 MeV lost end to end.

ONLY THE SLOPE IS FROM THE SIMULATION; THE ANCHOR IS THE DATA'S. The fitted intercept (183.6) is
not the 184.25 the simulation was generated at, so importing the absolute value would re-calibrate
the data's energy scale on a simulation artefact. The control is therefore a PIVOT:

    Ebeam(vz) = ebeam - (slope/1000) * (vz - pivot)          slope in MeV/m, vz and pivot in mm

At vz = pivot nothing changes, so the overall scale stays exactly where the data's calibration put
it and only the DRIFT is removed. Default pivot 500 mm, mid-chamber.

ON THE DATA (find805, the adopted selection) this takes the bound peak's spread across the chamber
from 0.300 MeV to 0.090 -- from swamping the whole 17C bound level scheme (g.s./0.217/0.331) to
below the g.s.-to-1/2+ spacing.

SIGN: applied backwards this DOUBLES the drift. The data's peak rises with vz and the simulation's
truth rises with vz, so the same sign applies -- an empirical check, not an assumption about which
end the beam enters. If a future cache mirrors vz, re-check it before trusting the knob.

  python add_ebvz.py <explorer.html> [slope MeV/m, default 11.0] [note shown under the formula]

The optional slope/note let another experiment set ITS OWN default (a1954 12Be at 600 torr H2:
CATIMA 9.1 MeV/m). With no extra arguments the page is exactly what it always was.
"""
import sys

p = sys.argv[1]
SLOPE = sys.argv[2] if len(sys.argv) > 2 else "11.0"
NOTE = sys.argv[3] if len(sys.argv) > 3 else "sim says 11.0 MeV/m"
s = open(p, encoding="utf-8").read()
orig = len(s)

def must(old, what):
    if old not in s:
        raise SystemExit("add_ebvz.py: anchor missing (%s) -- the template changed, re-check the patch" % what)

# 1. register the two numeric controls so readState() and the params save/load pick them up
must("'kcDenom','kcPivot','keOff'", "IDS list")
s = s.replace("'kcDenom','kcPivot','keOff'", "'kcDenom','kcPivot','keOff','ebzSlope','ebzPivot'", 1)

# 2. register the checkbox as a FLAG, same machinery as kcOn
must("const FLAGS = ['logY','logZ','kinLines','kcOn'];", "FLAGS list")
s = s.replace("const FLAGS = ['logY','logZ','kinLines','kcOn'];",
              "const FLAGS = ['logY','logZ','kinLines','kcOn','ebzOn'];", 1)

# 3. the UI, immediately after the theta-correction checkbox so the two corrections sit together
anchor = '<div class="row"><label class="check"><input type="checkbox" id="kcOn"> apply &theta; correction</label></div>'
must(anchor, "kcOn row")
blk = anchor + '''
      <div class="row"><label class="check"><input type="checkbox" id="ebzOn"> apply E<sub>beam</sub>(v<sub>z</sub>) &mdash; beam energy loss</label></div>
      <div class="row"><label for="ebzSlope">loss [MeV/m]</label><input type="number" id="ebzSlope" step="0.5" value="''' + SLOPE + '''"></div>
      <div class="row"><label for="ebzPivot">pivot v<sub>z</sub> [mm]</label><input type="number" id="ebzPivot" step="25" value="500"></div>
      <div class="row"><span class="mono num" style="font-size:11.5px;color:var(--ink-3)">E<sub>beam</sub> = ebeam &minus; loss&middot;(v<sub>z</sub>&minus;pivot); ''' + NOTE + '''</span></div>'''
s = s.replace(anchor, blk, 1)

# 4. the per-track beam energy, right after thCorr's helpers
must("function kcSlopeDeg(s) {", "kcSlopeDeg")
fn = '''/** Beam energy AT THE VERTEX. A constant Ebeam is what makes the reconstructed Ex drift with
    vz -- see this file's header. Pivoted so the overall energy scale is untouched and only the
    drift is removed. Returns s.ebeam unchanged when the knob is off or vz is unavailable. */
function ebAt(s, vz) {
  if (!s.ebzOn || vz === null || vz === undefined || !isFinite(vz)) return s.ebeam;
  const sl = isFinite(s.ebzSlope) ? s.ebzSlope : 0;
  const pv = isFinite(s.ebzPivot) ? s.ebzPivot : 500;
  return s.ebeam - (sl/1000)*(vz - pv);
}
function kcSlopeDeg(s) {'''
s = s.replace("function kcSlopeDeg(s) {", fn, 1)

# 5. both kinematics call sites. There are exactly two and BOTH must move together, or the
#    comparison overlay silently keeps the uncorrected Ex while the main histogram shifts --
#    the same failure add_keoff.py had to patch separately.
o1 = "const [ex] = kine2b(s.ebeam, thCorr(s, d.th[i], _k)*Math.PI/180, _k);"
must(o1, "exOnly kine2b")
s = s.replace(o1, "const [ex] = kine2b(ebAt(s, d.vz ? d.vz[i] : null), thCorr(s, d.th[i], _k)*Math.PI/180, _k);", 1)

o2 = "const [ex, tcm] = kine2b(s.ebeam, th*Math.PI/180, ke);"
must(o2, "compute kine2b")
s = s.replace(o2, "const [ex, tcm] = kine2b(ebAt(s, d.vz ? d.vz[i] : null), th*Math.PI/180, ke);", 1)

# 6. say so in the status line, beside the theta-correction note, so a reader of a screenshot can
#    tell whether the correction was on
o3 = "+ (s.kcOn ? ' &middot; <b>propagated to E<sub>x</sub></b>' : ' &middot; not propagated');"
if o3 in s:
    s = s.replace(o3,
        "+ (s.kcOn ? ' &middot; <b>propagated to E<sub>x</sub></b>' : ' &middot; not propagated')"
        "+ (s.ebzOn ? ' &middot; <b>E<sub>beam</sub>(v<sub>z</sub>) ON, ' + (isFinite(s.ebzSlope)?s.ebzSlope:0)"
        " + ' MeV/m</b>' : ' &middot; E<sub>beam</sub> constant');", 1)
    print("  status line: annotated")
else:
    print("  status line: anchor not found, skipped (cosmetic only)")

open(p, "w", encoding="utf-8").write(s)
print("  patched %s: +%d bytes" % (p, len(s) - orig))
print("  controls: ebzOn (checkbox), ebzSlope [MeV/m, default " + SLOPE + "], ebzPivot [mm, default 500]")
print("  DEFAULT IS OFF -- the page opens exactly as before until the box is ticked.")
