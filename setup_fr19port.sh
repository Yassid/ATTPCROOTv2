#!/bin/bash
# Environment for the FairRootv18.00-fairroot19-port worktree.
#
#   source setup_fr19port.sh
#
# This branch's build/config.sh predates the modern FairRoot config template: it
# leaves ROOT_INCLUDE_PATH empty, so cling cannot find any FairRoot or At* header
# and every class in a macro comes out as "unknown type name". Source this instead
# of config.sh on its own.

ATTPCROOT_DIR=/home/yassid/fair_install/ATTPCROOTv2_fr19port
SIMPATH_DIR=/home/yassid/fair_install/FairSoft/install
FAIRROOT_DIR=/home/yassid/fair_install/FairRoot_18.6   # the only FairRoot built against ROOT 6.26

source ${SIMPATH_DIR}/bin/thisroot.sh
export SIMPATH=${SIMPATH_DIR}
export FAIRROOTPATH=${FAIRROOT_DIR}
source ${ATTPCROOT_DIR}/build/config.sh > /dev/null 2>&1

export ROOT_INCLUDE_PATH=${ATTPCROOT_DIR}/install/include:${FAIRROOT_DIR}/include:${SIMPATH_DIR}/include/vmc:${SIMPATH_DIR}/include/root:${SIMPATH_DIR}/include

echo "ATTPCROOT (fr19port) : ${ATTPCROOT_DIR}"
echo "FairRoot             : ${FAIRROOTPATH}"
echo "ROOT                 : $(root-config --version 2>/dev/null)"
