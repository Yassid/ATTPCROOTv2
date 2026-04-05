"""Analyze GNN tracking performance and visualize results."""

import numpy as np
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

from dataset import build_dataset, train_val_split
from model import EdgeClassifier


def load_model(path="best_model.pt", device="cuda"):
    model = EdgeClassifier(in_channels=6, hidden=64, n_gravnet=4).to(device)
    model.load_state_dict(torch.load(path, weights_only=True, map_location=device))
    model.eval()
    return model


def predict_event(model, data, device, threshold=0.5):
    """Run inference on a single event, return clusters via connected components."""
    data = data.to(device)
    with torch.no_grad():
        scores = torch.sigmoid(model(data)).cpu().numpy()

    edge_index = data.edge_index.cpu().numpy()
    track_ids = data.track_id.cpu().numpy()
    pos = data.pos.cpu().numpy()
    n = data.num_nodes

    # Connected components on predicted edges
    mask = scores > threshold
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for j in range(edge_index.shape[1]):
        if mask[j]:
            a, b = edge_index[0, j], edge_index[1, j]
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

    pred_labels = np.array([find(i) for i in range(n)])
    # Relabel to 0, 1, 2, ...
    unique_labels = np.unique(pred_labels)
    label_map = {old: new for new, old in enumerate(unique_labels)}
    pred_labels = np.array([label_map[l] for l in pred_labels])

    return {
        "pos": pos,
        "track_ids": track_ids,
        "pred_labels": pred_labels,
        "edge_index": edge_index,
        "edge_scores": scores,
        "n_truth_tracks": len(np.unique(track_ids)),
        "n_pred_clusters": len(np.unique(pred_labels)),
    }


def compute_per_track_metrics(result):
    """Compute efficiency and purity per truth track."""
    track_ids = result["track_ids"]
    pred_labels = result["pred_labels"]
    metrics = []

    for tid in np.unique(track_ids):
        truth_mask = track_ids == tid
        truth_nodes = set(np.where(truth_mask)[0])

        # Find best matching cluster
        best_overlap, best_cluster, best_size = 0, -1, 0
        for cl in np.unique(pred_labels):
            cl_nodes = set(np.where(pred_labels == cl)[0])
            overlap = len(truth_nodes & cl_nodes)
            if overlap > best_overlap:
                best_overlap = overlap
                best_cluster = cl
                best_size = len(cl_nodes)

        eff = best_overlap / len(truth_nodes) if len(truth_nodes) > 0 else 0
        pur = best_overlap / best_size if best_size > 0 else 0

        metrics.append({
            "track_id": tid,
            "n_hits": len(truth_nodes),
            "efficiency": eff,
            "purity": pur,
            "best_cluster": best_cluster,
            "best_cluster_size": best_size,
            "n_missed": len(truth_nodes) - best_overlap,
        })

    return metrics


def investigate_efficiency(model, val_data, device, n_events=500):
    """Investigate where efficiency is lost."""
    all_track_metrics = []
    missed_positions = []  # positions of missed hits

    for i, data in enumerate(val_data[:n_events]):
        result = predict_event(model, data, device)
        track_metrics = compute_per_track_metrics(result)

        for tm in track_metrics:
            tm["event_idx"] = i
            all_track_metrics.append(tm)

            # Collect missed hit positions
            if tm["n_missed"] > 0:
                truth_mask = result["track_ids"] == tm["track_id"]
                pred_mask = result["pred_labels"] == tm["best_cluster"]
                missed_mask = truth_mask & ~pred_mask
                missed_positions.append(result["pos"][missed_mask])

    # Analyze
    beam_metrics = [m for m in all_track_metrics if m["track_id"] == 0]
    proton_metrics = [m for m in all_track_metrics if m["track_id"] != 0]

    print("=== Efficiency Analysis ===")
    print(f"Total tracks analyzed: {len(all_track_metrics)}")
    print(f"\nBeam tracks ({len(beam_metrics)}):")
    beam_eff = [m["efficiency"] for m in beam_metrics]
    print(f"  Efficiency: {np.mean(beam_eff):.4f} ± {np.std(beam_eff):.4f}")
    print(f"  Perfect (eff=1.0): {sum(1 for e in beam_eff if e > 0.99)}/{len(beam_eff)}")

    print(f"\nProton tracks ({len(proton_metrics)}):")
    proton_eff = [m["efficiency"] for m in proton_metrics]
    print(f"  Efficiency: {np.mean(proton_eff):.4f} ± {np.std(proton_eff):.4f}")
    print(f"  Perfect (eff=1.0): {sum(1 for e in proton_eff if e > 0.99)}/{len(proton_eff)}")

    # Efficiency vs number of hits
    print(f"\nEfficiency vs track size:")
    for lo, hi in [(0, 20), (20, 50), (50, 100), (100, 200), (200, 500), (500, 9999)]:
        subset = [m for m in all_track_metrics if lo <= m["n_hits"] < hi]
        if subset:
            eff = np.mean([m["efficiency"] for m in subset])
            print(f"  {lo:4d}-{hi:4d} hits: eff={eff:.4f} (N={len(subset)})")

    # Where are missed hits?
    if missed_positions:
        all_missed = np.vstack(missed_positions)
        print(f"\nMissed hit positions ({len(all_missed)} total):")
        print(f"  X: {all_missed[:,0].mean():.1f} ± {all_missed[:,0].std():.1f}")
        print(f"  Y: {all_missed[:,1].mean():.1f} ± {all_missed[:,1].std():.1f}")
        print(f"  Z: {all_missed[:,2].mean():.1f} ± {all_missed[:,2].std():.1f}")
        print(f"  R: {np.sqrt(all_missed[:,0]**2 + all_missed[:,1]**2).mean():.1f}")

    return all_track_metrics, missed_positions


def plot_event_3d(result, save_path=None, title=None):
    """3D event display: truth vs predicted clustering."""
    pos = result["pos"]
    truth = result["track_ids"]
    pred = result["pred_labels"]

    fig = plt.figure(figsize=(16, 6))

    # Truth
    ax1 = fig.add_subplot(131, projection='3d')
    for tid in np.unique(truth):
        mask = truth == tid
        label = f"Track {tid} ({mask.sum()} hits)"
        ax1.scatter(pos[mask, 0], pos[mask, 1], pos[mask, 2],
                    s=3, alpha=0.6, label=label)
    ax1.set_xlabel("X [mm]")
    ax1.set_ylabel("Y [mm]")
    ax1.set_zlabel("Z [mm]")
    ax1.set_title("MC Truth")
    ax1.legend(fontsize=7, loc='upper left')

    # Predicted
    ax2 = fig.add_subplot(132, projection='3d')
    for cl in np.unique(pred):
        mask = pred == cl
        ax2.scatter(pos[mask, 0], pos[mask, 1], pos[mask, 2],
                    s=3, alpha=0.6, label=f"Cluster {cl} ({mask.sum()})")
    ax2.set_xlabel("X [mm]")
    ax2.set_ylabel("Y [mm]")
    ax2.set_zlabel("Z [mm]")
    ax2.set_title(f"GNN Predicted ({result['n_pred_clusters']} clusters)")
    ax2.legend(fontsize=7, loc='upper left')

    # Errors: highlight misassigned hits
    ax3 = fig.add_subplot(133, projection='3d')
    # Map each truth track to its best cluster
    track_to_cluster = {}
    for tid in np.unique(truth):
        truth_mask = truth == tid
        best_cl, best_overlap = -1, 0
        for cl in np.unique(pred):
            overlap = (truth_mask & (pred == cl)).sum()
            if overlap > best_overlap:
                best_overlap = overlap
                best_cl = cl
        track_to_cluster[tid] = best_cl

    correct_mask = np.array([pred[i] == track_to_cluster[truth[i]] for i in range(len(truth))])
    ax3.scatter(pos[correct_mask, 0], pos[correct_mask, 1], pos[correct_mask, 2],
                s=3, alpha=0.3, c='gray', label=f"Correct ({correct_mask.sum()})")
    if (~correct_mask).sum() > 0:
        ax3.scatter(pos[~correct_mask, 0], pos[~correct_mask, 1], pos[~correct_mask, 2],
                    s=15, alpha=1.0, c='red', marker='x',
                    label=f"Misassigned ({(~correct_mask).sum()})")
    ax3.set_xlabel("X [mm]")
    ax3.set_ylabel("Y [mm]")
    ax3.set_zlabel("Z [mm]")
    ax3.set_title("Errors")
    ax3.legend(fontsize=7, loc='upper left')

    if title:
        fig.suptitle(title, fontsize=12)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def plot_event_projections(result, save_path=None, title=None):
    """2D projections: XY, XZ, YZ for truth vs predicted."""
    pos = result["pos"]
    truth = result["track_ids"]
    pred = result["pred_labels"]

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    projections = [
        (0, 1, "X [mm]", "Y [mm]", "XY"),
        (2, 0, "Z [mm]", "X [mm]", "XZ"),
        (2, 1, "Z [mm]", "Y [mm]", "YZ"),
    ]

    for col, (xi, yi, xlabel, ylabel, pname) in enumerate(projections):
        # Truth
        ax = axes[0, col]
        for tid in np.unique(truth):
            mask = truth == tid
            ax.scatter(pos[mask, xi], pos[mask, yi], s=2, alpha=0.6,
                       label=f"Track {tid}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"Truth {pname}")
        ax.legend(fontsize=7)

        # Predicted
        ax = axes[1, col]
        for cl in np.unique(pred):
            mask = pred == cl
            ax.scatter(pos[mask, xi], pos[mask, yi], s=2, alpha=0.6,
                       label=f"Cl {cl}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"Predicted {pname}")
        ax.legend(fontsize=7)

    if title:
        fig.suptitle(title, fontsize=13)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def plot_efficiency_summary(all_track_metrics, missed_positions, save_path=None):
    """Summary plots for efficiency analysis."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # 1. Efficiency distribution
    ax = axes[0, 0]
    beam = [m["efficiency"] for m in all_track_metrics if m["track_id"] == 0]
    proton = [m["efficiency"] for m in all_track_metrics if m["track_id"] != 0]
    ax.hist(beam, bins=50, range=(0, 1.05), alpha=0.7, label=f"Beam (N={len(beam)})")
    ax.hist(proton, bins=50, range=(0, 1.05), alpha=0.7, label=f"Proton (N={len(proton)})")
    ax.set_xlabel("Track efficiency")
    ax.set_ylabel("Tracks")
    ax.set_title("Efficiency distribution")
    ax.legend()

    # 2. Efficiency vs n_hits
    ax = axes[0, 1]
    nhits = [m["n_hits"] for m in all_track_metrics]
    effs = [m["efficiency"] for m in all_track_metrics]
    ax.scatter(nhits, effs, s=3, alpha=0.3)
    ax.set_xlabel("Number of hits")
    ax.set_ylabel("Efficiency")
    ax.set_title("Efficiency vs track size")

    # 3. Purity distribution
    ax = axes[0, 2]
    purs = [m["purity"] for m in all_track_metrics]
    ax.hist(purs, bins=50, range=(0, 1.05), color='green', alpha=0.7)
    ax.set_xlabel("Track purity")
    ax.set_ylabel("Tracks")
    ax.set_title(f"Purity (mean={np.mean(purs):.4f})")

    # 4. Missed hit Z distribution
    ax = axes[1, 0]
    if missed_positions:
        all_missed = np.vstack(missed_positions)
        ax.hist(all_missed[:, 2], bins=50, color='red', alpha=0.7)
        ax.set_xlabel("Z [mm]")
        ax.set_ylabel("Missed hits")
        ax.set_title(f"Missed hits Z distribution (N={len(all_missed)})")
    else:
        ax.text(0.5, 0.5, "No missed hits", ha='center', va='center', transform=ax.transAxes)

    # 5. Missed hit R distribution
    ax = axes[1, 1]
    if missed_positions:
        R = np.sqrt(all_missed[:, 0]**2 + all_missed[:, 1]**2)
        ax.hist(R, bins=50, color='red', alpha=0.7)
        ax.set_xlabel("R [mm]")
        ax.set_ylabel("Missed hits")
        ax.set_title("Missed hits R distribution")

    # 6. Number of missed hits per event
    ax = axes[1, 2]
    missed_per_event = {}
    for m in all_track_metrics:
        eidx = m["event_idx"]
        missed_per_event[eidx] = missed_per_event.get(eidx, 0) + m["n_missed"]
    vals = list(missed_per_event.values())
    ax.hist(vals, bins=50, color='orange', alpha=0.7)
    ax.set_xlabel("Missed hits per event")
    ax.set_ylabel("Events")
    ax.set_title(f"Missed hits/event (mean={np.mean(vals):.1f})")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load data and model
    data_list = build_dataset(
        "../data/gnn_training/all_hits.csv", k=10,
        cache_path="../data/gnn_training/graphs_k10.pt"
    )
    _, val_data = train_val_split(data_list)
    model = load_model("best_model.pt", device)
    print(f"Val events: {len(val_data)}")

    # 1. Investigate efficiency
    print("\n" + "=" * 60)
    track_metrics, missed_positions = investigate_efficiency(
        model, val_data, device, n_events=len(val_data)
    )

    # 2. Summary plots
    plot_efficiency_summary(track_metrics, missed_positions,
                            save_path="../data/gnn_training/efficiency_analysis.png")

    # 3. Visualize individual events
    # Good events
    for i in range(min(5, len(val_data))):
        result = predict_event(model, val_data[i], device)
        tm = compute_per_track_metrics(result)
        eff_str = " ".join(f"T{m['track_id']}:{m['efficiency']:.2f}" for m in tm)

        plot_event_projections(
            result,
            save_path=f"../data/gnn_training/event_{i}_projections.png",
            title=f"Event {i} | {result['n_pred_clusters']} clusters | {eff_str}"
        )

    # Find worst events (lowest proton efficiency)
    proton_by_event = {}
    for m in track_metrics:
        if m["track_id"] != 0:
            proton_by_event[m["event_idx"]] = m["efficiency"]

    worst_events = sorted(proton_by_event.items(), key=lambda x: x[1])[:5]
    print(f"\nWorst proton events: {worst_events}")

    for rank, (eidx, eff) in enumerate(worst_events):
        result = predict_event(model, val_data[eidx], device)
        tm = compute_per_track_metrics(result)
        eff_str = " ".join(f"T{m['track_id']}:{m['efficiency']:.2f}" for m in tm)

        plot_event_projections(
            result,
            save_path=f"../data/gnn_training/worst_{rank}_ev{eidx}_projections.png",
            title=f"Worst #{rank+1} (event {eidx}) | eff={eff:.2f} | {eff_str}"
        )

    print("\nDone!")


if __name__ == "__main__":
    main()
