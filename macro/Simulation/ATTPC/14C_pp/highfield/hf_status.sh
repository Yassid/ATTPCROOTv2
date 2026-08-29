#!/usr/bin/env bash
# What the campaign has actually finished, counted from the COMPLETED markers rather than from
# files existing -- a killed job leaves files behind.
ROOTDIR=${HF_ROOT:-/mnt/f/a1954_C14_hf}
printf "%-14s %8s  %s\n" config done levels
for cfg in b285_attpc b285_2mm b400_attpc b400_2mm b700_attpc b700_2mm; do
   d="$ROOTDIR/$cfg"
   n=$(ls "$d"/*.marker 2>/dev/null | wc -l)
   lv=$(ls "$d"/*.marker 2>/dev/null | sed 's#.*/##; s/_s[0-9]*_.*//' | sort -u | tr '\n' ' ')
   printf "%-14s %5s/5  %s\n" "$cfg" "$n" "${lv:--}"
done
echo
du -sh "$ROOTDIR" 2>/dev/null
