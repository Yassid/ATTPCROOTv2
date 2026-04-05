"""Export GNN-separated raw hits per track for C++ ClusterizeSmooth3D.

Instead of doing arc-walk clustering in Python, this script only uses the
GNN for track separation (edge classification + connected components) and
exports the raw hit positions per predicted track. The C++ side then applies
ClusterizeSmooth3D to each pure-track hit set.

Usage:
    source ~/gnn_env/bin/activate
    python gnn_hits_export.py [--n-events 5002]
"""

import argparse
import numpy as np
import pandas as pd
import torch

from dataset import build_dataset, train_val_split
from model import EdgeClassifier


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


def export_hits(data_list, model, device, output_hits, output_tracks,
                csv_path="../data/gnn_training/all_hits.csv"):
    """Run GNN track separation and export raw hits per track."""

    # Load real charges from CSV
    print("Loading hit charges from CSV...")
    df = pd.read_csv(csv_path)
    # Build per-event charge lookup: event_charges[event_id] = array of charges
    event_charges = {}
    for eid, group in df.groupby("event_id"):
        event_charges[eid] = group["charge"].values

    f_ht = open(output_hits, 'w')
    f_ht.write("event_id,hit_id,track_pred,x,y,z,charge\n")

    f_tr = open(output_tracks, 'w')
    f_tr.write("event_id,track_pred,n_hits,geo_radius,geo_theta,geo_phi,"
               "geo_cx,geo_cy\n")

    for ev_idx in range(len(data_list)):
        data = data_list[ev_idx]
        pos = data.pos.cpu().numpy()
        n = data.num_nodes

        # Use original digi event_id (not sequential graph index)
        orig_event_id = int(data.event_id) if hasattr(data, 'event_id') else ev_idx

        # Get real charges for this event (keyed by original event_id)
        charges = event_charges.get(orig_event_id, np.ones(n))

        # GNN inference: predict edge scores
        edge_index = data.edge_index.cpu().numpy()
        with torch.no_grad():
            scores = torch.sigmoid(model(data.to(device))).cpu().numpy()

        # Union-find on edges with score > 0.5 → track separation
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for j in range(edge_index.shape[1]):
            if scores[j] > 0.5:
                a, b = edge_index[0, j], edge_index[1, j]
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb

        # Map nodes to sequential track IDs
        root_to_track = {}
        track_counter = 0
        node_track = np.full(n, -1, dtype=int)
        for node in range(n):
            root = find(node)
            if root not in root_to_track:
                root_to_track[root] = track_counter
                track_counter += 1
            node_track[node] = root_to_track[root]

        # Write raw hits per track with real charges
        for hit_id in range(n):
            tp = node_track[hit_id]
            q = charges[hit_id] if hit_id < len(charges) else 1.0
            f_ht.write(f"{orig_event_id},{hit_id},{tp},"
                       f"{pos[hit_id, 0]:.6f},{pos[hit_id, 1]:.6f},{pos[hit_id, 2]:.6f},"
                       f"{q:.4f}\n")

        # Write per-track geometry (circle fit for Brho seeding)
        for tp in range(track_counter):
            tp_mask = node_track == tp
            tp_pos = pos[tp_mask]

            if len(tp_pos) < 5:
                continue

            # Circle fit in XY
            try:
                cx, cy, R = fit_circle_xy(tp_pos[:, 0], tp_pos[:, 1])
            except Exception:
                cx, cy, R = 0, 0, 0

            # Initial direction from first few hits (Z-sorted)
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

            f_tr.write(f"{orig_event_id},{tp},{tp_mask.sum()},"
                       f"{R:.4f},{theta:.4f},{phi:.4f},{cx:.4f},{cy:.4f}\n")

        if (ev_idx + 1) % 500 == 0:
            print(f"  Processed {ev_idx + 1}/{len(data_list)} events...")

    f_ht.close()
    f_tr.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-events", type=int, default=5002)
    parser.add_argument("--model", default="best_model.pt")
    parser.add_argument("--cache", default="../data/gnn_training/graphs_k10.pt")
    parser.add_argument("--csv", default="../data/gnn_training/all_hits.csv")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    data_list = build_dataset(args.csv, k=10, cache_path=args.cache)
    if args.n_events > 0:
        data_list = data_list[:args.n_events]
    print(f"Exporting hits for {len(data_list)} events...")

    model = load_model(args.model, device)

    out_dir = "../data/gnn_training"
    export_hits(
        data_list, model, device,
        output_hits=f"{out_dir}/gnn_separated_hits.csv",
        output_tracks=f"{out_dir}/gnn_separated_tracks.csv",
        csv_path=args.csv,
    )

    print(f"\nExported to:")
    print(f"  {out_dir}/gnn_separated_hits.csv   (raw hits per GNN track)")
    print(f"  {out_dir}/gnn_separated_tracks.csv  (per-track geometry)")


if __name__ == "__main__":
    main()
