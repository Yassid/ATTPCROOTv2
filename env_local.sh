#!/bin/bash
# Local WSL build environment for ATTPCROOTv2 FairRootv18.00 branch
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh

export SIMPATH=/home/yassid/fair_install/FairSoft/install
export FAIRROOTPATH=/home/yassid/fair_install/FairRootInstall
export GENFIT=/home/yassid/fair_install/GenFit
export HDF5_ROOT=/home/yassid/fair_install/hdf5

export CMAKE_PREFIX_PATH=$FAIRROOTPATH:$SIMPATH:$HDF5_ROOT:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=$GENFIT/lib:$HDF5_ROOT/lib:$LD_LIBRARY_PATH
