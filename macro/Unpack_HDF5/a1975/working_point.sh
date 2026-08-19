#!/usr/bin/env bash
# MOVED 2026-08-19. This file is no longer the working point; it is a pointer to it.
#
# The a1975 analysis now lives in one repository, with the data it is measured from:
#
#     ~/a1975_analysis        (p,p) (p,d) (d,t) (p,t), one cache per channel
#
# The constants are in a1975_analysis/working_point.sh. Source THAT. Keeping a second copy here
# is exactly the failure this file was created to prevent -- and it had already happened: the copy
# here carried A1975_DT_EBEAM=184.17 while every adopted (d,t) producer used 184.25.
#
#     source ~/a1975_analysis/working_point.sh
#     ~/a1975_analysis/common/check_consistency.sh
#
# This shim forwards, so an old call still works.
A1975_REPO="${A1975_REPO:-$HOME/a1975_analysis}"
if [ -r "$A1975_REPO/working_point.sh" ]; then
   # shellcheck disable=SC1091
   . "$A1975_REPO/working_point.sh"
else
   echo "a1975 working point not found at $A1975_REPO/working_point.sh" >&2
   echo "clone https://github.com/Yassid/a1975_analysis or set A1975_REPO" >&2
   return 1 2>/dev/null || exit 1
fi
