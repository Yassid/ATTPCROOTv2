"""Export GNN arc-walk clusters to CSV for C++ UKF integration.

Reads the cached graph data, runs GNN inference + arc-walk clustering,
and writes cluster centers with covariance to CSV files that the
run_ukf_gnn.C macro can read.

Usage:
    source ~/gnn_env/bin/activate
    python gnn_cluster_export.py [--hits-per-cluster 8] [--n-events 10000]
"""

import argparse
import numpy as np
import torch

from dataset import build_dataset, train_val_split
from model import EdgeClassifier
from cluster import cluster_arc_walk, cluster_positions


# Diffusion parameters — must match AtTrackTransformer.cxx
COEF_T = 0.00009      # cm^2/us  (transverse diffusion)
COEF_L = 0.0000009    # cm^2/us  (longitudinal diffusion)
DRIFT_VEL = 1.0       # cm/us
TB_TIME = 0.320        # us
PAD_RES_XY = 2.3       # mm
PAD_RES_Z = 2.3 * 1.5  # mm


def load_model(path="best_model.pt", device="cuda"):
    model = EdgeClassifier(in_channels=6, hidden=64, n_gravnet=4).to(device)
    model.load_state_dict(torch.load(path, weights_only=True, map_location=device))
    model.eval()
    return model


def fit_circle_xy(x, y):
    """Algebraic circle fit (Kasa method). Returns (cx, cy, R)."""
    A = np.column_stack([2 * x, 2 * y, np.ones_like(x)])
    b = x**2 + y**2
    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = result[0], result[1]
    R = np.sqrt(result[2] + cx**2 + cy**2)
    return cx, cy, R


def compute_cluster_covariance(hit_positions, hit_charges):
    """Compute cluster covariance matching AtTrackTransformer.cxx lines 110-127.

    Returns (cov_xx, cov_yy, cov_zz) — diagonal elements of the 3x3 covariance.
    """
    var_x, var_y, var_z = 0.0, 0.0, 0.0
    total_q = 0.0

    for pos, q in zip(hit_positions, hit_charges):
        drift_time = pos[2] / (10.0 * DRIFT_VEL)  # us
        varT = 100.0 * COEF_T * 2.0 * drift_time   # mm^2
        varL = 100.0 * COEF_L * 2.0 * drift_time   # mm^2
        tbRes_mm = DRIFT_VEL * TB_TIME * 10.0       # mm per TB
        varTB = tbRes_mm * tbRes_mm / 12.0           # mm^2

        var_x += q * q * (PAD_RES_XY**2 + varT)
        var_y += q * q * (PAD_RES_XY**2 + varT)
        var_z += q * q * (PAD_RES_Z**2 + varTB + varL)
        total_q += q

    if total_q < 1e-9:
        return PAD_RES_XY**2, PAD_RES_XY**2, PAD_RES_Z**2

    sigma_x = np.sqrt(var_x) / total_q
    sigma_y = np.sqrt(var_y) / total_q
    sigma_z = np.sqrt(var_z) / total_q
    return sigma_x**2, sigma_y**2, sigma_z**2


def compute_charge_weighted_center(hit_positions, hit_charges):
    """Charge-weighted mean position (X,Y charge-weighted, Z simple mean)."""
    total_q = sum(hit_charges)
    if total_q < 1e-9:
        return np.mean(hit_positions, axis=0)
    x = sum(p[0] * q for p, q in zip(hit_positions, hit_charges)) / total_q
    y = sum(p[1] * q for p, q in zip(hit_positions, hit_charges)) / total_q
    z = np.mean([p[2] for p in hit_positions])
    return np.array([x, y, z])


def export_clusters(data_list, model, device, output_clusters, output_hits,
                    output_tracks, hits_per_cluster=8):
    """Run GNN arc-walk on all events and export to CSV."""

    f_cl = open(output_clusters, 'w')
    f_cl.write("event_id,cluster_id,track_pred,n_hits,x,y,z,charge,"
               "cov_xx,cov_yy,cov_zz\n")

    f_ht = open(output_hits, 'w')
    f_ht.write("event_id,hit_id,cluster_id,track_pred\n")

    f_tr = open(output_tracks, 'w')
    f_tr.write("event_id,track_pred,n_hits,n_clusters,geo_radius,geo_theta,geo_phi,"
               "geo_cx,geo_cy\n")

    for ev_idx in range(len(data_list)):
        data = data_list[ev_idx]
        pos = data.pos.cpu().numpy()
        n = data.num_nodes

        # Get raw charges from un-normalized features
        # x_norm = (x_raw - means) / stds, but we need raw charge
        # Feature order: x, y, z, charge, r, phi
        # We can get charge from pos and the raw data, but the normalized
        # features don't directly give raw charge. Use a default charge of 1.0
        # for covariance computation (the relative weighting is what matters).
        # Actually, we have x_raw stored in data.pos for x,y,z.
        # For charge, use uniform weighting since we don't have raw charge in
        # the graph data structure.
        charges = np.ones(n)  # uniform charge weighting

        # Run arc-walk clustering
        r = cluster_arc_walk(model, data, device, hits_per_cluster=hits_per_cluster)
        labels = r["labels"]

        # Identify tracks (connected components from arc-walk)
        track_labels = set()
        # arc-walk assigns cluster IDs sequentially per-track, so we need to
        # map cluster IDs back to track IDs. Use the track structure from
        # the arc-walk internals: clusters from the same track have consecutive IDs.
        # We can recover track membership from which nodes share clusters.

        # Build track_pred from the arc-walk's connected components
        # (the arc-walk already did union-find internally; we reconstruct from labels)
        # Nodes with same cluster belong to same track. Nodes with label=-1 are noise.
        # Group clusters that belong to the same connected component.

        # Simpler: re-run union-find on the predicted edges to get track IDs
        edge_index = data.edge_index.cpu().numpy()
        with torch.no_grad():
            scores = torch.sigmoid(model(data.to(device))).cpu().numpy()

        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        threshold = 0.5
        for j in range(edge_index.shape[1]):
            if scores[j] > threshold:
                a, b = edge_index[0, j], edge_index[1, j]
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb

        # Map each node to its track_pred (normalized to sequential IDs)
        root_to_track = {}
        track_counter = 0
        node_track = np.full(n, -1, dtype=int)
        for node in range(n):
            root = find(node)
            if root not in root_to_track:
                root_to_track[root] = track_counter
                track_counter += 1
            node_track[node] = root_to_track[root]

        # For each cluster, determine its track_pred (majority vote)
        cluster_track = {}
        for cl in range(r["n_clusters"]):
            cl_mask = labels == cl
            if cl_mask.sum() == 0:
                continue
            tracks_in_cl = node_track[cl_mask]
            cluster_track[cl] = np.bincount(tracks_in_cl[tracks_in_cl >= 0]).argmax()

        # Compute all cluster centers and covariances
        cluster_data = []
        for cl in range(r["n_clusters"]):
            cl_mask = labels == cl
            if cl_mask.sum() == 0:
                continue

            hit_indices = np.where(cl_mask)[0]
            hit_pos = pos[hit_indices]
            hit_q = charges[hit_indices]
            tp = cluster_track.get(cl, -1)

            center = compute_charge_weighted_center(hit_pos, hit_q)
            cov_xx, cov_yy, cov_zz = compute_cluster_covariance(hit_pos, hit_q)

            cluster_data.append({
                "cl": cl, "tp": tp, "n_hits": len(hit_indices),
                "center": center, "charge": hit_q.sum(),
                "cov_xx": cov_xx, "cov_yy": cov_yy, "cov_zz": cov_zz,
            })

        # Write per-cluster data, ordered per track with vertex (high Z) first.
        # The UKF expects clusters ordered from vertex (highest Z in digi frame)
        # toward the Bragg peak (lowest Z in digi frame).
        from itertools import groupby
        cluster_data.sort(key=lambda c: c["tp"])
        new_cl_id = 0
        for tp, group in groupby(cluster_data, key=lambda c: c["tp"]):
            track_clusters = sorted(group, key=lambda c: -c["center"][2])  # descending Z
            for cd in track_clusters:
                f_cl.write(f"{ev_idx},{new_cl_id},{cd['tp']},{cd['n_hits']},"
                           f"{cd['center'][0]:.6f},{cd['center'][1]:.6f},{cd['center'][2]:.6f},"
                           f"{cd['charge']:.2f},{cd['cov_xx']:.6f},{cd['cov_yy']:.6f},{cd['cov_zz']:.6f}\n")
                new_cl_id += 1

        # Write per-hit mapping
        for hit_id in range(n):
            cl = labels[hit_id]
            tp = node_track[hit_id]
            f_ht.write(f"{ev_idx},{hit_id},{cl},{tp}\n")

        # Write per-track geometry (circle fit for Brho seeding)
        tracks_in_event = set(node_track[node_track >= 0])
        for tp in sorted(tracks_in_event):
            tp_mask = node_track == tp
            tp_pos = pos[tp_mask]

            if len(tp_pos) < 5:
                continue

            # Circle fit in XY
            try:
                cx, cy, R = fit_circle_xy(tp_pos[:, 0], tp_pos[:, 1])
            except Exception:
                cx, cy, R = 0, 0, 0

            # Initial direction from first few hits (sorted by Z descending)
            z_order = np.argsort(-tp_pos[:, 2])
            tp_sorted = tp_pos[z_order]
            n_use = min(5, len(tp_sorted))
            if n_use >= 2:
                dr = tp_sorted[n_use - 1] - tp_sorted[0]
                norm = np.linalg.norm(dr)
                if norm > 1e-6:
                    dr /= norm
                    if dr[2] > 0:
                        dr = -dr
                    r_xy = np.sqrt(dr[0]**2 + dr[1]**2)
                    theta = np.degrees(np.arctan2(r_xy, dr[2]))
                    phi = np.degrees(np.arctan2(dr[1], dr[0]))
                else:
                    theta, phi = 0, 0
            else:
                theta, phi = 0, 0

            # Count clusters for this track
            n_cl = sum(1 for c, t in cluster_track.items() if t == tp)

            f_tr.write(f"{ev_idx},{tp},{tp_mask.sum()},{n_cl},"
                       f"{R:.4f},{theta:.4f},{phi:.4f},{cx:.4f},{cy:.4f}\n")

        if (ev_idx + 1) % 500 == 0:
            print(f"  Processed {ev_idx + 1}/{len(data_list)} events...")

    f_cl.close()
    f_ht.close()
    f_tr.close()


def main():
    parser = argparse.ArgumentParser(description="Export GNN arc-walk clusters to CSV")
    parser.add_argument("--hits-per-cluster", type=int, default=8)
    parser.add_argument("--n-events", type=int, default=10000,
                        help="Max events to process (0 = all)")
    parser.add_argument("--model", default="best_model.pt")
    parser.add_argument("--cache", default="../data/gnn_training/graphs_k10.pt")
    parser.add_argument("--csv", default="../data/gnn_training/all_hits.csv")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load data — use ALL events (not just val split) for full export
    data_list = build_dataset(args.csv, k=10, cache_path=args.cache)
    if args.n_events > 0:
        data_list = data_list[:args.n_events]
    print(f"Exporting clusters for {len(data_list)} events...")

    model = load_model(args.model, device)

    out_dir = "../data/gnn_training"
    export_clusters(
        data_list, model, device,
        output_clusters=f"{out_dir}/gnn_clusters.csv",
        output_hits=f"{out_dir}/gnn_hit_to_cluster.csv",
        output_tracks=f"{out_dir}/gnn_tracks.csv",
        hits_per_cluster=args.hits_per_cluster,
    )

    print(f"\nExported to:")
    print(f"  {out_dir}/gnn_clusters.csv     (cluster centers + covariance)")
    print(f"  {out_dir}/gnn_hit_to_cluster.csv (hit-to-cluster mapping)")
    print(f"  {out_dir}/gnn_tracks.csv        (per-track geometry for Brho seeding)")


if __name__ == "__main__":
    main()
