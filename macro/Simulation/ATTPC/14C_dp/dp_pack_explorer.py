#!/usr/bin/env python3
"""Pack the dp_export_explorer.C dumps into one self-contained HTML explorer.

    root -b -q 'dp_export_explorer.C("/mnt/f/a1954_C14dp","<bindir>")'
    ./dp_pack_explorer.py <bindir> [out.html]

Paths are arguments on purpose: hard-coding a session scratchpad is what silently broke several
a1975 explorer launchers once those sessions ended.
"""
import base64, json, os, sys
d = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dpexp"
samples = []
for line in open(os.path.join(d, "index.txt")):
    cfg, lev, base, n = line.split()
    raw = open(os.path.join(d, base + ".bin"), "rb").read()
    samples.append({"cfg": cfg, "lev": lev, "n": int(n),
                    "b64": base64.b64encode(raw).decode()})
tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "dp_explorer_template.html")).read()
out = tpl.replace("__DATA__", json.dumps(samples, separators=(",", ":")))
p = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/a1954_C14dp_explorer.html")
open(p, "w").write(out)
print("wrote", p, round(len(out)/1e6, 2), "MB,", sum(s["n"] for s in samples), "events")
