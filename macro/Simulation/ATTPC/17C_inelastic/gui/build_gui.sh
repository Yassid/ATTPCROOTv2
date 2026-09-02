#!/usr/bin/env bash
# Build the standalone 17C kinematics explorer: one self-contained HTML file, no network, no server.
#
#   ./build_gui.sh [output.html]
#
# Two steps, because the data and the page have very different lifetimes. export_gui_C17.C reads the
# campaign products (which live outside the repo, on INEL_ROOT) and writes C17_gui_data.json; this
# script inlines that JSON into template.html. Edit the template freely and re-run -- the export only
# has to be redone when the campaign changes.
#
# The JSON is committed so the page can be rebuilt by anyone with the repo, since the ~40 GB of
# campaign products it came from cannot be.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUT=${1:-$HOME/C17_kinematics_explorer.html}
[ -s "$HERE/C17_gui_data.json" ] || {
   echo "missing C17_gui_data.json -- run:  root -b -q 'export_gui_C17.C'  from the parent directory"
   exit 2
}
python3 - "$HERE" "$OUT" <<'PY'
import sys
here, out = sys.argv[1], sys.argv[2]
tpl = open(f"{here}/template.html").read()
data = open(f"{here}/C17_gui_data.json").read().strip()
if "__DATA__" not in tpl:
    raise SystemExit("template.html has no __DATA__ placeholder")
open(out, "w").write(tpl.replace("__DATA__", data))
print(f"wrote {out}  ({len(tpl)+len(data):,} bytes)")
PY
