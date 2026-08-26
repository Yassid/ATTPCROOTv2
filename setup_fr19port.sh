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


# Geant4 data. build/config.sh on this branch sets none of these, so any simulation dies
# with "G4ENSDFSTATEDATA environment variable must be set" before the first event.
G4DATA=${SIMPATH_DIR}/share/Geant4-11.0.1/data
export G4ENSDFSTATEDATA=${G4DATA}/G4ENSDFSTATE2.3
export G4LEVELGAMMADATA=${G4DATA}/PhotonEvaporation5.7
export G4RADIOACTIVEDATA=${G4DATA}/RadioactiveDecay5.6
export G4LEDATA=${G4DATA}/G4EMLOW8.0
export G4NEUTRONHPDATA=${G4DATA}/G4NDL4.6
export G4PARTICLEXSDATA=${G4DATA}/G4PARTICLEXS4.0
export G4PIIDATA=${G4DATA}/G4PII1.3
export G4SAIDXSDATA=${G4DATA}/G4SAIDDATA2.0
export G4ABLADATA=${G4DATA}/G4ABLA3.1
export G4INCLDATA=${G4DATA}/G4INCL1.0
export G4REALSURFACEDATA=${G4DATA}/RealSurface2.2

echo "ATTPCROOT (fr19port) : ${ATTPCROOT_DIR}"
echo "FairRoot             : ${FAIRROOTPATH}"
echo "ROOT                 : $(root-config --version 2>/dev/null)"
