#!/usr/bin/env python3
"""Set the generated explorer's DEFAULT control values.

    python3 pp/set_defaults_C15p.py <explorer.html> [key=value ...]

Why this exists: the page ships with icLo/icHi/npLo/npHi controls whose defaults silently discard
most of the sample. npulse has three states and only two are physics --

    npulse == 1   single pulse
    npulse >  1   genuine pile-up
    npulse == 0   NO IC VALUE WAS JOINED, not a pile-up verdict

so the shipped npLo=npHi=1 throws away every track the IC join could not reach. On the (p,d)
cache that is 11,269 of 18,264 tracks, and the page looks four times emptier than an older one
that has no IC controls at all.

It is worse than a display choice here: the kin cache's ic values come from a join on
(run, event, trackID), and gate_events_C15p.C RENUMBERS events from zero when it reduces a run to
its gated subset. The key means different things on the two sides, so the ic values that did join
are matched to the wrong events. Filtering on them is filtering on noise. The beam gate has
ALREADY been applied upstream -- gate_events only kept IC-passing events before fitting -- so
opening these controls wide loses no beam selection.

Defaults written here are therefore the honest ones: no IC filtering, and exBins matched to the
older pages so spectra can be compared by eye.
"""
import re, sys

DEFAULTS = {"npLo": "0", "npHi": "9999", "icLo": "-1", "icHi": "99999", "exBins": "200"}

path = sys.argv[1]
for kv in sys.argv[2:]:
    k, _, v = kv.partition("=")
    DEFAULTS[k] = v

s = open(path, encoding="utf-8").read()
for k, v in DEFAULTS.items():
    pat = re.compile(r'(id="%s"[^>]*?value=")[^"]*(")' % re.escape(k))
    s, n = pat.subn(lambda m: m.group(1) + v + m.group(2), s, count=1)
    if not n:
        pat = re.compile(r'(value=")[^"]*("[^>]*?id="%s")' % re.escape(k))
        s, n = pat.subn(lambda m: m.group(1) + v + m.group(2), s, count=1)
    print("  %-8s -> %-8s %s" % (k, v, "ok" if n else "CONTROL NOT FOUND"))
open(path, "w", encoding="utf-8").write(s)
print("defaults written to %s" % path)
