#!/bin/bash
# Move the a2091 15C data products onto the Seagate and leave symlinks behind.
#
# WHY SYMLINKS instead of rewriting paths: /home/yassid/a2091_C15_* is hard-coded in 77
# places across 38 files (macros and shell scripts), and the external mount point contains
# a SPACE ("/media/yassid/Seagate Hub/..."). Rewriting all of them would risk breaking every
# unquoted shell expansion. Replacing each directory with a symlink keeps all existing paths
# working verbatim while the bytes live on the external drive. The symlinks themselves sit on
# ext4, so exfat's lack of symlink support does not matter.
#
# Verifies every source file exists at the destination with the same size BEFORE deleting
# any original.
#
# NB: run this with its log OUTSIDE the directories being migrated. Logging into
# /home/yassid/a2091_C15_ic/ once cost the second half of the output: rsync copied the log,
# rm -rf removed the original, and the still-open fd then pointed at a deleted inode, so
# every later line went nowhere.
set -u
SG="/media/yassid/Seagate Hub/ATTPC"
say(){ echo "[$(date '+%H:%M:%S')] $*"; }

mkdir -p "$SG" || { say "cannot write to $SG"; exit 1; }

migrate(){                      # $1 = directory name under /home/yassid
  local name="$1" src="/home/yassid/$1" dst="$SG/$1"
  if [ -L "$src" ]; then say "$name already a symlink -> $(readlink "$src")"; return 0; fi
  if [ ! -d "$src" ]; then say "$name: no local directory, skipping"; return 0; fi

  local nsrc ssrc
  nsrc=$(find "$src" -type f | wc -l)
  ssrc=$(du -sb "$src" | cut -f1)
  say "$name: $nsrc files, $(numfmt --to=iec "$ssrc") -> external"
  mkdir -p "$dst"
  rsync -a "$src"/ "$dst"/ || { say "$name: RSYNC FAILED, originals kept"; return 1; }

  # Verify PER FILE: every source file must exist at the destination with the same size.
  # Do NOT compare file COUNTS -- the destination may legitimately hold more (e.g. the fit
  # directory already contained fits migrated earlier, 150 vs 31, which made an equality
  # check reject a perfectly good copy).
  local miss=0 rel s1 s2
  while IFS= read -r f; do
    rel="${f#"$src"/}"
    s1=$(stat -c%s "$f")
    s2=$(stat -c%s "$dst/$rel" 2>/dev/null || echo -1)
    if [ "$s1" != "$s2" ]; then say "  MISSING/DIFFERS: $rel ($s1 vs $s2)"; miss=$((miss+1)); fi
  done < <(find "$src" -type f)
  if [ "$miss" -ne 0 ]; then say "$name: $miss file(s) failed verification -- originals kept"; return 1; fi

  say "$name: verified all $nsrc files present at destination; replacing with symlink"
  rm -rf "$src" && ln -s "$dst" "$src" || { say "$name: symlink step FAILED"; return 1; }
  say "$name: done -> $(readlink "$src")"
}

# a2091_C15_fit was already migrated by hand earlier; just link it if the local dir lingers.
for d in a2091_C15_reco a2091_C15_ic a2091_C15_fit; do migrate "$d"; done

echo
say "RESULT:"
for d in a2091_C15_reco a2091_C15_ic a2091_C15_fit; do
  printf "  %-20s -> %s\n" "$d" "$(readlink "/home/yassid/$d" 2>/dev/null || echo '(still local)')"
done
df -h /home | tail -1
df -h "$SG" | tail -1
