#!/bin/bash
cd /mnt/f/a1975/sim_rebuild/trig
ls /mnt/f/a1975/gf_dp_find805/ | sed -nE 's/^run_0*([0-9]+)_multifit_genfitter_p\.root$/\1/p' | xargs -P 6 -I{} ~/Spyral/venv/bin/python -I mesh_features.py {}
echo FEATURES_ALLDONE
