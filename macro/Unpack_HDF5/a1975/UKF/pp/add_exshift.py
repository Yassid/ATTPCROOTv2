#!/usr/bin/env python3
"""Add an additive Ex SHIFT control (exShift, MeV) to a generated explorer page.

WHY. The adopted 16C(p,d) chain (a1975_analysis/pd/run_pd_chain.sh) applies EXSHIFT = -0.6309 MeV
after the two-body kinematics: the g.s. sat at +0.6309(76) on the adopted configuration and the
shift has to be additive (dEx/dEbeam is only 0.005 here; keOff would need ~-2 MeV and wrecks the
width). The explorer had no such control, so the page showed the g.s. at +0.63 while the chain
showed it at 0.

It is added inside kine2b(), so EVERY Ex the page computes (main histogram, comparison overlay,
Ex-vs-theta_cm, slices) gets it. theta_cm is untouched. The KE-vs-theta reference loci are physical
and are not shifted. Registered in IDS, so it saves and loads with the params JSON.

    python3 add_exshift.py page.html [default_MeV]
"""
import sys

path = sys.argv[1]
default = sys.argv[2] if len(sys.argv) > 2 else "0"
s = open(path, encoding="utf-8").read()
if 'id="exShift"' in s:
    print("  exShift already present -- not patched again"); sys.exit(0)

row_anchor = '<div class="row"><label for="keOff">KE offset [MeV]</label>'
i = s.find(row_anchor)
if i < 0: sys.exit("ERROR: keOff row not found")
j = s.find('</div>', i) + len('</div>')
row = ('\n      <div class="row" title="additive Ex shift after the kinematics (run_pd_chain.sh EXSHIFT)">'
       f'<label for="exShift">Ex shift [MeV]</label><input type="number" id="exShift" step="0.01" value="{default}"></div>')
s = s[:j] + row + s[j:]

old = "const m4x = Math.sqrt(a), ex = m4x - M4;"
if s.count(old) != 1: sys.exit("ERROR: kine2b Ex line not found exactly once")
s = s.replace(old, "const m4x = Math.sqrt(a), ex = m4x - M4 + EXSHIFT();")
s = s.replace("function kine2b(Kp, thl, Ke) {",
              "let _exsEl = null;\nfunction EXSHIFT() { _exsEl = _exsEl || document.getElementById('exShift');"
              " const v = _exsEl ? parseFloat(_exsEl.value) : 0; return isFinite(v) ? v : 0; }\n"
              "function kine2b(Kp, thl, Ke) {", 1)

ids_old = "'kcDenom','kcPivot','keOff',"
if ids_old not in s: sys.exit("ERROR: IDS list not found")
s = s.replace(ids_old, "'kcDenom','kcPivot','keOff','exShift',", 1)

open(path, "w", encoding="utf-8").write(s)
print(f"  patched {path}: exShift control, default {default} MeV, applied inside kine2b")
