#!/usr/bin/env bash
REPO=~/fair_install/ATTPCROOTv2-OpenKF; FAIR=~/fair_install/FairRootInstall; SIM=~/fair_install/FairSoft/install
g++ -std=c++17 -O2 match_diag_prog.cpp -o match_diag $(root-config --cflags --libs) \
  -I$REPO/build/include -I$FAIR/include -I$SIM/include -L$REPO/build/lib -lAtData -Wl,-rpath,$REPO/build/lib
