"""Evaluate arc-walk clustering quality for UKF kinematics.

Compares cluster centers from arc-walk against truth to check:
  1. Track range preservation (sum of inter-cluster distances)
  2. Initial direction (theta, phi from first clusters)
  3. Circle fit in XY → radius → Brho → KE estimate
  4. Spatial residuals of cluster centers vs truth hits
"""

import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from dataset import build_dataset, train_val_split
from model import EdgeClassifier
from cluster import cluster_arc_walk, cluster_positions


def load_model(path="best_model.pt", device="cuda"):
    model = EdgeClassifier(in_channels=6, hidden=64, n_gravnet=4).to(device)
    model.load_state_dict(torch.load(path, weights_only=True, map_location=device))
    model.eval()
    return model


def get_truth_kinematics(csv_path, event_ids):
    """Extract truth KE and initial direction per proton track from CSV."""
    df = pd.read_csv(csv_path)
    truth = {}
    for eid in event_ids:
        ev = df[df.event_id == eid]
        proton = ev[ev.track_id == 1]
        if len(proton) == 0:
            continue
        truth[eid] = {
            "KE": proton["energy"].iloc[0],
            "A": proton["A"].iloc[0],
            "Z_particle": proton["Z"].iloc[0],
        }
    return truth


def fit_circle_xy(x, y):
    """Algebraic circle fit in XY plane. Returns (cx, cy, R)."""
    # Kasa method
    A = np.column_stack([2 * x, 2 * y, np.ones_like(x)])
    b = x**2 + y**2
    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy = result[0], result[1]
    R = np.sqrt(result[2] + cx**2 + cy**2)
    return cx, cy, R


def range_from_centers(centers_ordered):
    """Total arc length from ordered cluster centers."""
    if len(centers_ordered) < 2:
        return 0.0
    diffs = np.diff(centers_ordered, axis=0)
    return np.sum(np.linalg.norm(diffs, axis=1))


def direction_from_hits(hits, n_use=5):
    """Estimate initial direction from the first few hits.

    Uses the direction vector from hit[0] to hit[n_use-1], which matches
    how the UKF seeds the initial state. Hits must be pre-sorted by arc
    order (vertex end first).
    """
    if len(hits) < 3:
        return None, None
    n = min(n_use, len(hits))

    # Direction from first to nth hit
    dr = hits[n - 1] - hits[0]
    norm = np.linalg.norm(dr)
    if norm < 1e-6:
        return None, None
    dr /= norm

    # Orient: travel direction is from high Z to low Z in digi frame
    if dr[2] > 0:
        dr = -dr

    # theta from Z axis, phi in XY
    r_xy = np.sqrt(dr[0]**2 + dr[1]**2)
    theta = np.arctan2(r_xy, dr[2])
    phi = np.arctan2(dr[1], dr[0])
    return np.degrees(theta), np.degrees(phi)


def initial_direction_from_clusters(pos, labels, proton_clusters, track_ids,
                                     n_hits=10):
    """Extract initial direction using raw hits from the first clusters.

    1. Collect hits from proton clusters near the vertex (high Z end)
    2. Sort those hits by Z (descending)
    3. Fit a line through the first n_hits to get direction
    """
    # Collect all proton hits that belong to proton clusters
    proton_hit_indices = []
    for cl in proton_clusters:
        cl_mask = labels == cl
        hit_indices = np.where(cl_mask)[0]
        # Only keep hits that are actually proton
        for idx in hit_indices:
            if track_ids[idx] == 1:
                proton_hit_indices.append(idx)

    if len(proton_hit_indices) < 3:
        return None, None

    proton_hit_indices = np.array(proton_hit_indices)
    proton_pos = pos[proton_hit_indices]

    # Sort by Z descending (vertex = high Z in digi frame)
    z_order = np.argsort(-proton_pos[:, 2])
    proton_pos = proton_pos[z_order]

    return direction_from_hits(proton_pos, n_use=n_hits)


def truth_initial_direction(pos, track_ids, tid=1, n_hits=10):
    """Estimate truth initial direction using PCA on first hits."""
    mask = track_ids == tid
    t_pos = pos[mask]
    if len(t_pos) < 3:
        return None, None
    # Sort by Z descending (vertex = high Z in digi frame)
    z_order = np.argsort(-t_pos[:, 2])
    t_pos = t_pos[z_order]
    return direction_from_hits(t_pos, n_use=n_hits)


def truth_range(pos, track_ids, tid=1):
    """Compute truth track range from sorted hits."""
    mask = track_ids == tid
    t_pos = pos[mask]
    if len(t_pos) < 2:
        return 0.0
    # Sort by Z
    z_order = np.argsort(-t_pos[:, 2])
    t_pos = t_pos[z_order]
    diffs = np.diff(t_pos, axis=0)
    return np.sum(np.linalg.norm(diffs, axis=1))


def evaluate_arc_walk_kinematics(model, val_data, device, csv_path,
                                  n_events=200, B_field=2.85,
                                  hits_per_cluster=8):
    """Run arc-walk clustering and compare kinematics with truth.

    B_field: magnetic field in Tesla (AT-TPC nominal ~2.85 T)
    """
    # Get truth KE from CSV
    # We need to map val_data indices back to event IDs
    # Since val_data is shuffled, we need the CSV for truth KE
    # but the graph data has track positions

    results = {
        "truth_KE": [], "range_KE": [], "circle_KE": [],
        "truth_theta": [], "cluster_theta": [],
        "truth_range": [], "cluster_range": [],
        "truth_radius": [], "cluster_radius": [],
        "n_proton_clusters": [], "proton_coverage": [],
    }

    proton_mass = 938.272  # MeV/c^2

    for i in range(min(n_events, len(val_data))):
        data = val_data[i]
        track_ids = data.track_id.cpu().numpy()
        pos = data.pos.cpu().numpy()

        # Skip events without proton track
        if 1 not in track_ids:
            continue

        # Run arc-walk clustering
        r = cluster_arc_walk(model, data, device, hits_per_cluster=hits_per_cluster)
        labels = r["labels"]

        # Find proton clusters (clusters where majority of hits are track 1)
        proton_clusters = []
        for cl in range(r["n_clusters"]):
            cl_mask = labels == cl
            if cl_mask.sum() == 0:
                continue
            cl_tracks = track_ids[cl_mask]
            if np.bincount(cl_tracks).argmax() == 1:  # majority proton
                proton_clusters.append(cl)

        if len(proton_clusters) < 3:
            continue

        # Get ordered cluster centers for proton track
        centers = cluster_positions(pos, labels)
        proton_centers = np.array([centers[cl] for cl in sorted(proton_clusters)])

        # Order by Z (descending = vertex end first in digi frame)
        z_order = np.argsort(-proton_centers[:, 2])
        proton_centers = proton_centers[z_order]

        # 1. Track range
        c_range = range_from_centers(proton_centers)
        t_range = truth_range(pos, track_ids, tid=1)
        results["cluster_range"].append(c_range)
        results["truth_range"].append(t_range)

        # 2. Initial direction (from raw hits, not cluster centers)
        c_theta, _ = initial_direction_from_clusters(
            pos, labels, proton_clusters, track_ids, n_hits=10)
        t_theta, _ = truth_initial_direction(pos, track_ids, tid=1, n_hits=10)
        if c_theta is not None and t_theta is not None:
            results["cluster_theta"].append(c_theta)
            results["truth_theta"].append(t_theta)

        # 3. Circle fit in XY → radius → KE
        if len(proton_centers) >= 5:
            try:
                cx, cy, R_mm = fit_circle_xy(proton_centers[:, 0],
                                              proton_centers[:, 1])
                R_m = R_mm / 1000.0  # convert to meters
                # p = qBR in natural units: p[MeV/c] = 0.3 * Z * B[T] * R[m] * 1000
                # Actually p[MeV/c] = 0.3 * Z * B[T] * R[mm]
                p_MeV = 0.3 * 1 * B_field * R_mm  # Z=1 for proton
                KE_circle = np.sqrt(p_MeV**2 + proton_mass**2) - proton_mass
                results["circle_KE"].append(KE_circle)
                results["cluster_radius"].append(R_mm)

                # Truth circle fit
                proton_mask = track_ids == 1
                t_pos = pos[proton_mask]
                _, _, R_truth = fit_circle_xy(t_pos[:, 0], t_pos[:, 1])
                results["truth_radius"].append(R_truth)
                p_truth = 0.3 * 1 * B_field * R_truth
                KE_truth_circle = np.sqrt(p_truth**2 + proton_mass**2) - proton_mass
                results["truth_KE"].append(KE_truth_circle)
            except Exception:
                pass

        # 4. Coverage
        proton_mask = track_ids == 1
        proton_in_clusters = 0
        for cl in proton_clusters:
            cl_mask = labels == cl
            proton_in_clusters += (track_ids[cl_mask] == 1).sum()
        coverage = proton_in_clusters / proton_mask.sum()
        results["proton_coverage"].append(coverage)
        results["n_proton_clusters"].append(len(proton_clusters))

    return results


def plot_kinematics_eval(results, save_path=None):
    """Plot kinematic comparison: arc-walk clusters vs truth."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 1. Range comparison
    ax = axes[0, 0]
    t_range = np.array(results["truth_range"])
    c_range = np.array(results["cluster_range"])
    ax.scatter(t_range, c_range, s=5, alpha=0.5)
    lim = max(t_range.max(), c_range.max()) * 1.05
    ax.plot([0, lim], [0, lim], 'r--', lw=1)
    ax.set_xlabel("Truth range [mm]")
    ax.set_ylabel("Cluster range [mm]")
    ax.set_title(f"Track Range (N={len(t_range)})")
    ratio = c_range / np.clip(t_range, 1, None)
    ax.text(0.05, 0.9, f"ratio: {np.mean(ratio):.3f}±{np.std(ratio):.3f}",
            transform=ax.transAxes, fontsize=9)

    # 2. Theta comparison
    ax = axes[0, 1]
    t_theta = np.array(results["truth_theta"])
    c_theta = np.array(results["cluster_theta"])
    ax.scatter(t_theta, c_theta, s=5, alpha=0.5)
    lim_t = [min(t_theta.min(), c_theta.min()) - 5,
             max(t_theta.max(), c_theta.max()) + 5]
    ax.plot(lim_t, lim_t, 'r--', lw=1)
    ax.set_xlabel("Truth θ [deg]")
    ax.set_ylabel("Cluster θ [deg]")
    dtheta = c_theta - t_theta
    ax.set_title(f"Initial θ (Δ={np.mean(dtheta):.2f}±{np.std(dtheta):.2f}°)")

    # 3. Circle radius comparison
    ax = axes[0, 2]
    if len(results["truth_radius"]) > 0:
        t_R = np.array(results["truth_radius"])
        c_R = np.array(results["cluster_radius"])
        ax.scatter(t_R, c_R, s=5, alpha=0.5)
        lim = max(t_R.max(), c_R.max()) * 1.05
        ax.plot([0, lim], [0, lim], 'r--', lw=1)
        ax.set_xlabel("Truth R [mm]")
        ax.set_ylabel("Cluster R [mm]")
        dR = (c_R - t_R) / t_R * 100
        ax.set_title(f"Circle Radius (ΔR/R={np.mean(dR):.1f}±{np.std(dR):.1f}%)")

    # 4. KE from circle fit
    ax = axes[1, 0]
    if len(results["truth_KE"]) > 0:
        t_KE = np.array(results["truth_KE"])
        c_KE = np.array(results["circle_KE"])
        ax.scatter(t_KE, c_KE, s=5, alpha=0.5)
        lim = max(t_KE.max(), c_KE.max()) * 1.05
        ax.plot([0, lim], [0, lim], 'r--', lw=1)
        ax.set_xlabel("Truth KE (circle fit) [MeV]")
        ax.set_ylabel("Cluster KE (circle fit) [MeV]")
        dKE = (c_KE - t_KE) / np.clip(t_KE, 0.1, None) * 100
        ax.set_title(f"KE from Circle Fit (Δ={np.mean(dKE):.1f}±{np.std(dKE):.1f}%)")

    # 5. Range ratio distribution
    ax = axes[1, 1]
    ax.hist(ratio, bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(1.0, color='r', ls='--')
    ax.set_xlabel("Cluster range / Truth range")
    ax.set_ylabel("Events")
    ax.set_title(f"Range Ratio Distribution")

    # 6. Coverage and cluster count
    ax = axes[1, 2]
    cov = np.array(results["proton_coverage"])
    n_cl = np.array(results["n_proton_clusters"])
    ax.hist(cov, bins=30, edgecolor='black', alpha=0.7, label=f"coverage")
    ax.set_xlabel("Proton hit coverage fraction")
    ax.set_ylabel("Events")
    ax.set_title(f"Coverage: {np.mean(cov):.3f}±{np.std(cov):.3f}, "
                 f"clusters/ev: {np.mean(n_cl):.1f}")

    fig.suptitle("Arc-Walk Clustering: Kinematic Quality Check", fontsize=14)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def plot_theta_residuals(results, save_path=None):
    """Detailed theta residual analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    t_theta = np.array(results["truth_theta"])
    c_theta = np.array(results["cluster_theta"])
    dtheta = c_theta - t_theta

    # Residual distribution
    ax = axes[0]
    ax.hist(dtheta, bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(0, color='r', ls='--')
    ax.set_xlabel("Δθ (cluster - truth) [deg]")
    ax.set_ylabel("Events")
    ax.set_title(f"θ Residual: μ={np.mean(dtheta):.2f}°, σ={np.std(dtheta):.2f}°")

    # Residual vs truth theta
    ax = axes[1]
    ax.scatter(t_theta, dtheta, s=5, alpha=0.5)
    ax.axhline(0, color='r', ls='--')
    ax.set_xlabel("Truth θ [deg]")
    ax.set_ylabel("Δθ [deg]")
    ax.set_title("θ Residual vs Truth θ")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    csv_path = "../data/gnn_training/all_hits.csv"
    data_list = build_dataset(
        csv_path, k=10,
        cache_path="../data/gnn_training/graphs_k10.pt"
    )
    _, val_data = train_val_split(data_list)
    model = load_model("best_model.pt", device)

    print(f"Evaluating arc-walk kinematics on {min(200, len(val_data))} events...")
    results = evaluate_arc_walk_kinematics(
        model, val_data, device, csv_path,
        n_events=200, B_field=2.85, hits_per_cluster=8
    )

    print(f"\n{'='*60}")
    print(f"Arc-Walk Kinematic Quality Summary")
    print(f"{'='*60}")

    t_range = np.array(results["truth_range"])
    c_range = np.array(results["cluster_range"])
    ratio = c_range / np.clip(t_range, 1, None)
    print(f"Range ratio (cluster/truth): {np.mean(ratio):.3f} ± {np.std(ratio):.3f}")

    cov = np.array(results["proton_coverage"])
    print(f"Proton coverage: {np.mean(cov):.3f} ± {np.std(cov):.3f}")

    n_cl = np.array(results["n_proton_clusters"])
    print(f"Proton clusters/event: {np.mean(n_cl):.1f} ± {np.std(n_cl):.1f}")

    if len(results["truth_theta"]) > 0:
        dtheta = np.array(results["cluster_theta"]) - np.array(results["truth_theta"])
        print(f"θ residual: {np.mean(dtheta):.2f} ± {np.std(dtheta):.2f} deg")

    if len(results["truth_radius"]) > 0:
        t_R = np.array(results["truth_radius"])
        c_R = np.array(results["cluster_radius"])
        dR = (c_R - t_R) / t_R * 100
        print(f"Radius residual (ΔR/R): {np.mean(dR):.1f} ± {np.std(dR):.1f} %")

    if len(results["truth_KE"]) > 0:
        t_KE = np.array(results["truth_KE"])
        c_KE = np.array(results["circle_KE"])
        dKE = (c_KE - t_KE) / np.clip(t_KE, 0.1, None) * 100
        print(f"KE residual (circle fit): {np.mean(dKE):.1f} ± {np.std(dKE):.1f} %")

    # Plots
    out = "../data/gnn_training"
    plot_kinematics_eval(results, save_path=f"{out}/arcwalk_kinematics.png")
    plot_theta_residuals(results, save_path=f"{out}/arcwalk_theta_residuals.png")

    # Also try different stride values
    print(f"\n{'='*60}")
    print("Stride scan...")
    for stride in [4, 8, 12, 16]:
        res = evaluate_arc_walk_kinematics(
            model, val_data, device, csv_path,
            n_events=200, hits_per_cluster=stride
        )
        r = np.array(res["cluster_range"]) / np.clip(np.array(res["truth_range"]), 1, None)
        cov = np.array(res["proton_coverage"])
        n_cl = np.array(res["n_proton_clusters"])
        if len(res["truth_radius"]) > 0:
            dR = (np.array(res["cluster_radius"]) - np.array(res["truth_radius"])) / np.array(res["truth_radius"]) * 100
            print(f"  stride={stride:2d}: range_ratio={np.mean(r):.3f}, "
                  f"coverage={np.mean(cov):.3f}, "
                  f"clusters/ev={np.mean(n_cl):.1f}, "
                  f"ΔR/R={np.mean(dR):.1f}±{np.std(dR):.1f}%")


if __name__ == "__main__":
    main()
