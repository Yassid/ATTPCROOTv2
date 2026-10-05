#!/usr/bin/env python3
"""Make the per-run Ex map honour the KE-offset and Ebeam(vz) knobs.

    python3 pp/fix_runpanel_corr_C15p.py <explorer.html>

inject_run_panel.py adds its own build() with its own kine2b call, and add_keoff / add_ebvz patch
by exact string match against the TEMPLATE's call sites -- so the run panel's copy is invisible to
both. The result is a panel that keeps computing the uncorrected Ex while the main histogram
shifts under the knobs: the same class of divergence add_ebvz's own comments warn about for the
comparison overlay, and just as silent, because a per-run map that is subtly wrong still looks
entirely reasonable.

Run AFTER inject_run_panel, add_keoff and add_ebvz -- it needs ebAt() to exist and it patches the
panel that inject_run_panel wrote.
"""
import sys

p = sys.argv[1]
s = open(p, encoding="utf-8").read()

OLD = "const [ex] = kine2b(sAll.ebeam, thCorr(sAll, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);"
NEW = ("const _kr = d.ke[i] + (sAll.keOff||0);\n"
       "      const [ex] = kine2b(typeof ebAt === 'function' ? ebAt(sAll, d.vz ? d.vz[i] : null) : sAll.ebeam,"
       " thCorr(sAll, d.th[i], _kr)*Math.PI/180, _kr);")

if OLD not in s:
    if "_kr" in s:
        print("  run panel already corrected -- nothing to do")
        sys.exit(0)
    raise SystemExit("fix_runpanel_corr: anchor missing -- inject_run_panel's build() changed, "
                     "re-check that the per-run map still honours the knobs")

n = s.count(OLD)
s = s.replace(OLD, NEW)
open(p, "w", encoding="utf-8").write(s)
print("  per-run Ex map now honours keOff and Ebeam(vz)  (%d call site%s)" % (n, "" if n == 1 else "s"))
