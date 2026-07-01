"""Spyral pipeline config for a1975 16C(d,p)17C — DEUTERIUM target, PID up to Estimation.

For the Spyral<->ATTPCROOT clustering comparison on BACKWARD (d,p) protons.
Parameters from the original C16_dp Spyral config (programs/Spyral-1.0/main.py),
which was written for exactly these D2 runs. The ONE D2-critical difference vs the
proton-target Spyral config is window_time_bucket=560 (slower D2 drift, 1.15 cm/us
vs 1.30; encodes the z-calibration that decides forward/backward).

Runs PointcloudPhase -> ClusterPhase -> EstimationPhase (no solver). Spyral
parallelizes across runs, so feed it chunked "runs" (see chunk_run.py).

Run:  /home/yassid/Spyral/venv/bin/python <this file>
Override run range / processes / trace / workspace via env vars:
  SPY_RUN_MIN, SPY_RUN_MAX, SPY_NPROC, SPY_TRACE, SPY_WS, SPY_ACTIVE (e.g. 011 to reuse pointcloud)
"""
import os
from pathlib import Path
import multiprocessing

from spyral import (
    Pipeline,
    start_pipeline,
    PointcloudPhase,
    ClusterPhase,
    EstimationPhase,
    PadParameters,
    GetParameters,
    FribParameters,
    DetectorParameters,
    ClusterParameters,
    OverlapJoinParameters,
    EstimateParameters,
    DEFAULT_MAP,
)
from spyral.core.config import HdbscanParameters

trace_path = Path(os.environ.get("SPY_TRACE", "/home/yassid/spyral_d2/h5"))
workspace_path = Path(os.environ.get("SPY_WS", "/home/yassid/spyral_d2/workspace"))

run_min = int(os.environ.get("SPY_RUN_MIN", "300"))
run_max = int(os.environ.get("SPY_RUN_MAX", "305"))
n_processes = int(os.environ.get("SPY_NPROC", "6"))

# phase activity mask: "111" all phases, "011" reuse existing Pointcloud
_active = os.environ.get("SPY_ACTIVE", "111")
active = [c == "1" for c in _active]

pad_params = PadParameters(
    pad_geometry_path=DEFAULT_MAP,
    pad_time_path=DEFAULT_MAP,
    pad_scale_path=DEFAULT_MAP,
)

get_params = GetParameters(
    baseline_window_scale=20.0,
    peak_separation=50.0,
    peak_prominence=20.0,
    peak_max_width=50.0,
    peak_threshold=40.0,
)

frib_params = FribParameters(
    baseline_window_scale=100.0,
    peak_separation=50.0,
    peak_prominence=20.0,
    peak_max_width=500.0,
    peak_threshold=100.0,
    ic_delay_time_bucket=1100,
    ic_multiplicity=1,
)

det_params = DetectorParameters(
    magnetic_field=2.85,
    electric_field=45000.0,
    detector_length=1000.0,
    beam_region_radius=25.0,
    micromegas_time_bucket=10.0,
    window_time_bucket=560.0,   # D2-critical: slower drift than proton target (490)
    get_frequency=6.25,
    garfield_file_path=Path("invalid.txt"),
    do_garfield_correction=False,
)

cluster_params = ClusterParameters(
    min_cloud_size=50,
    hdbscan_parameters=HdbscanParameters(
        min_points=3,
        min_size_scale_factor=0.05,
        min_size_lower_cutoff=10,
        cluster_selection_epsilon=10.0,
    ),
    tripclust_parameters=None,
    overlap_join=OverlapJoinParameters(
        circle_overlap_ratio=0.25,
        min_cluster_size_join=15,
    ),
    continuity_join=None,
    direction_threshold=0.5,
    outlier_scale_factor=0.05,
)

estimate_params = EstimateParameters(
    min_total_trajectory_points=30, smoothing_factor=100.0
)

pipe = Pipeline(
    [
        PointcloudPhase(get_params, frib_params, det_params, pad_params),
        ClusterPhase(cluster_params, det_params),
        EstimationPhase(estimate_params, det_params),
    ],
    active,
    workspace_path,
    trace_path,
)


def main():
    start_pipeline(pipe, run_min, run_max, n_processes)


if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    main()
