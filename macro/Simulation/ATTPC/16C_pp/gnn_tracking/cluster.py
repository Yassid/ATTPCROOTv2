"""GNN-based clustering for UKF track fitting.

Three clustering strategies using the trained GNN:
  Option 1: HDBSCAN in learned embedding space
  Option 2: Soft clustering head (future — requires retraining)
  Option 3: Arc-walk on weighted graph, fixed-stride grouping
"""

import numpy as np
import torch
import hdbscan
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dataset import build_dataset, train_val_split
from model import EdgeClassifier


def load_model(path="best_model.pt", device="cuda"):
    model = EdgeClassifier(in_channels=6, hidden=64, n_gravnet=4).to(device)
    model.load_state_dict(torch.load(path, weights_only=True, map_location=device))
    model.eval()
    return model


# ============================================================
# Option 1: HDBSCAN in embedding space
# ============================================================

def cluster_embedding_hdbscan(model, data, device, min_cluster_size=5, min_samples=3):
    """Cluster hits using HDBSCAN on GNN-learned embeddings.

    Returns dict with cluster labels, embeddings, and positions.
    """
    data_dev = data.to(device)
    with torch.no_grad():
        _, embeddings = model(data_dev, return_embeddings=True)

    emb = embeddings.cpu().numpy()
    pos = data.pos.cpu().numpy()
    track_ids = data.track_id.cpu().numpy()

    # HDBSCAN on learned embeddings
    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="euclidean",
    )
    labels = clusterer.fit_predict(emb)

    return {
        "labels": labels,
        "embeddings": emb,
        "pos": pos,
        "track_ids": track_ids,
        "n_clusters": len(set(labels)) - (1 if -1 in labels else 0),
        "n_noise": (labels == -1).sum(),
        "probabilities": clusterer.probabilities_,
    }


# ============================================================
# Option 3: Arc-walk on weighted graph + fixed-stride grouping
# ============================================================

def cluster_arc_walk(model, data, device, hits_per_cluster=10):
    """Cluster by walking along the track arc using edge scores as weights.

    1. Predict edge scores
    2. For each track (from connected components on high-score edges):
       a. Build subgraph of same-track edges
       b. Find maximum spanning tree → natural arc ordering
       c. Walk the tree from one endpoint to the other
       d. Group every N consecutive hits into one cluster
    """
    data_dev = data.to(device)
    with torch.no_grad():
        scores = torch.sigmoid(model(data_dev)).cpu().numpy()

    edge_index = data.edge_index.cpu().numpy()
    pos = data.pos.cpu().numpy()
    track_ids = data.track_id.cpu().numpy()
    n = data.num_nodes

    # Step 1: connected components on edges with score > 0.5
    threshold = 0.5
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for j in range(edge_index.shape[1]):
        if scores[j] > threshold:
            a, b = edge_index[0, j], edge_index[1, j]
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

    # Group nodes by track
    tracks = {}
    for node in range(n):
        tracks.setdefault(find(node), []).append(node)

    # Step 2: for each track, build weighted adjacency and find arc order
    cluster_labels = np.full(n, -1, dtype=int)
    cluster_id = 0

    for track_root, track_nodes in tracks.items():
        if len(track_nodes) < 3:
            continue

        node_set = set(track_nodes)
        node_to_local = {node: i for i, node in enumerate(track_nodes)}
        tn = len(track_nodes)

        # Build adjacency with edge scores as weights
        adj = [[] for _ in range(tn)]
        for j in range(edge_index.shape[1]):
            a, b = edge_index[0, j], edge_index[1, j]
            if a in node_set and b in node_set and scores[j] > threshold:
                la, lb = node_to_local[a], node_to_local[b]
                adj[la].append((lb, scores[j]))
                adj[lb].append((la, scores[j]))

        # Maximum spanning tree via Prim's algorithm
        visited = [False] * tn
        tree_adj = [[] for _ in range(tn)]

        # Start from the node with highest Z (vertex end in digi frame)
        start = max(range(tn), key=lambda i: pos[track_nodes[i], 2])
        visited[start] = True
        heap = []

        import heapq
        for nb, w in adj[start]:
            heapq.heappush(heap, (-w, start, nb))  # max heap via negation

        while heap:
            neg_w, u, v = heapq.heappop(heap)
            if visited[v]:
                continue
            visited[v] = True
            tree_adj[u].append(v)
            tree_adj[v].append(u)
            for nb, w in adj[v]:
                if not visited[nb]:
                    heapq.heappush(heap, (-w, v, nb))

        # Walk the tree from start (DFS to find longest path = arc)
        def dfs_farthest(root):
            dist = [-1] * tn
            dist[root] = 0
            stack = [root]
            order = []
            while stack:
                node = stack.pop()
                order.append(node)
                for nb in tree_adj[node]:
                    if dist[nb] == -1:
                        dist[nb] = dist[node] + 1
                        stack.append(nb)
            farthest = max(range(tn), key=lambda i: dist[i])
            return farthest, order

        # Find endpoints of the longest path
        end1, _ = dfs_farthest(start)
        end2, arc_order = dfs_farthest(end1)

        # Group consecutive hits into clusters
        for i in range(0, len(arc_order), hits_per_cluster):
            group = arc_order[i:i + hits_per_cluster]
            for local_idx in group:
                cluster_labels[track_nodes[local_idx]] = cluster_id
            cluster_id += 1

    return {
        "labels": cluster_labels,
        "pos": pos,
        "track_ids": track_ids,
        "n_clusters": cluster_id,
        "n_noise": (cluster_labels == -1).sum(),
    }


# ============================================================
# Evaluation: compare clustering with truth and ClusterizeSmooth3D
# ============================================================

def cluster_positions(pos, labels):
    """Compute charge-weighted cluster centers (equal weight here)."""
    clusters = {}
    for i, cl in enumerate(labels):
        if cl < 0:
            continue
        clusters.setdefault(cl, []).append(pos[i])
    centers = {}
    for cl, pts in clusters.items():
        centers[cl] = np.mean(pts, axis=0)
    return centers


def evaluate_clustering(result, method_name=""):
    """Print clustering summary."""
    labels = result["labels"]
    track_ids = result["track_ids"]
    pos = result["pos"]

    n_clusters = result["n_clusters"]
    n_noise = result["n_noise"]
    n_hits = len(labels)

    print(f"\n{'='*50}")
    print(f"Method: {method_name}")
    print(f"  Hits: {n_hits}, Clusters: {n_clusters}, Noise: {n_noise}")

    # Per truth-track breakdown
    for tid in np.unique(track_ids):
        mask = track_ids == tid
        t_labels = labels[mask]
        t_clusters = set(t_labels[t_labels >= 0])
        n_in_clusters = (t_labels >= 0).sum()
        print(f"  Track {tid}: {mask.sum()} hits → {len(t_clusters)} clusters, "
              f"{n_in_clusters} clustered, {(t_labels == -1).sum()} noise")

    # Cluster size distribution
    sizes = []
    for cl in range(n_clusters):
        s = (labels == cl).sum()
        sizes.append(s)
    if sizes:
        print(f"  Cluster sizes: min={min(sizes)}, max={max(sizes)}, "
              f"mean={np.mean(sizes):.1f}, median={np.median(sizes):.0f}")

    # Track purity per cluster
    purities = []
    for cl in range(n_clusters):
        cl_mask = labels == cl
        if cl_mask.sum() == 0:
            continue
        cl_tracks = track_ids[cl_mask]
        dominant = np.bincount(cl_tracks).max()
        purities.append(dominant / cl_mask.sum())
    if purities:
        print(f"  Cluster purity: {np.mean(purities):.4f}")


def plot_clustering_comparison(results, event_idx, save_path=None):
    """Side-by-side 2D projections for different clustering methods."""
    n_methods = len(results)
    fig, axes = plt.subplots(n_methods + 1, 3, figsize=(18, 5 * (n_methods + 1)))

    pos = results[0][1]["pos"]
    track_ids = results[0][1]["track_ids"]

    projections = [
        (0, 1, "X", "Y"),
        (2, 0, "Z", "X"),
        (2, 1, "Z", "Y"),
    ]

    # Row 0: truth
    for col, (xi, yi, xl, yl) in enumerate(projections):
        ax = axes[0, col]
        for tid in np.unique(track_ids):
            mask = track_ids == tid
            ax.scatter(pos[mask, xi], pos[mask, yi], s=2, alpha=0.6,
                       label=f"Track {tid}")
        ax.set_xlabel(f"{xl} [mm]")
        ax.set_ylabel(f"{yl} [mm]")
        ax.set_title(f"Truth {xl}{yl}")
        ax.legend(fontsize=7)

    # Subsequent rows: each method
    for row, (name, result) in enumerate(results, 1):
        labels = result["labels"]
        for col, (xi, yi, xl, yl) in enumerate(projections):
            ax = axes[row, col]

            # Noise
            noise_mask = labels == -1
            if noise_mask.sum() > 0:
                ax.scatter(pos[noise_mask, xi], pos[noise_mask, yi],
                           s=2, alpha=0.3, c='gray', label="noise")

            # Clusters
            unique_labels = sorted(set(labels[labels >= 0]))
            colors = plt.cm.tab20(np.linspace(0, 1, max(len(unique_labels), 1)))
            for ci, cl in enumerate(unique_labels):
                mask = labels == cl
                ax.scatter(pos[mask, xi], pos[mask, yi],
                           s=4, alpha=0.7, color=colors[ci % len(colors)])

            # Draw cluster centers
            centers = cluster_positions(pos, labels)
            cx = [centers[cl][xi] for cl in sorted(centers)]
            cy = [centers[cl][yi] for cl in sorted(centers)]
            ax.scatter(cx, cy, s=30, marker='x', c='black', linewidths=1, zorder=5)

            ax.set_xlabel(f"{xl} [mm]")
            ax.set_ylabel(f"{yl} [mm]")
            ax.set_title(f"{name} {xl}{yl} ({result['n_clusters']} cl)")

    fig.suptitle(f"Event {event_idx} — Clustering Comparison", fontsize=14)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


def plot_embeddings_2d(result, save_path=None):
    """Visualize HDBSCAN clustering in 2D PCA of embedding space."""
    from sklearn.decomposition import PCA

    emb = result["embeddings"]
    labels = result["labels"]
    track_ids = result["track_ids"]

    pca = PCA(n_components=2)
    emb_2d = pca.fit_transform(emb)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Color by truth
    ax = axes[0]
    for tid in np.unique(track_ids):
        mask = track_ids == tid
        ax.scatter(emb_2d[mask, 0], emb_2d[mask, 1], s=3, alpha=0.6,
                   label=f"Track {tid}")
    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    ax.set_title("Embedding space (colored by truth)")
    ax.legend(fontsize=8)

    # Color by HDBSCAN cluster
    ax = axes[1]
    noise_mask = labels == -1
    if noise_mask.sum() > 0:
        ax.scatter(emb_2d[noise_mask, 0], emb_2d[noise_mask, 1],
                   s=3, alpha=0.3, c='gray', label="noise")
    for cl in sorted(set(labels[labels >= 0])):
        mask = labels == cl
        ax.scatter(emb_2d[mask, 0], emb_2d[mask, 1], s=3, alpha=0.6,
                   label=f"Cl {cl}")
    ax.set_xlabel("PCA 1")
    ax.set_ylabel("PCA 2")
    ax.set_title("Embedding space (colored by HDBSCAN)")
    ax.legend(fontsize=7)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {save_path}")
    plt.close()


# ============================================================
# Main: compare all methods
# ============================================================

def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    data_list = build_dataset(
        "../data/gnn_training/all_hits.csv", k=10,
        cache_path="../data/gnn_training/graphs_k10.pt"
    )
    _, val_data = train_val_split(data_list)
    model = load_model("best_model.pt", device)

    # Run on several events
    for ev_idx in [0, 1, 5, 10, 50]:
        if ev_idx >= len(val_data):
            break

        data = val_data[ev_idx]
        print(f"\n{'#'*60}")
        print(f"Event {ev_idx}: {data.num_nodes} hits, "
              f"tracks: {data.track_id.unique().tolist()}")

        # Option 1: HDBSCAN in embedding space
        r1 = cluster_embedding_hdbscan(model, data, device,
                                        min_cluster_size=5, min_samples=3)
        evaluate_clustering(r1, "Option 1: HDBSCAN (embedding)")

        # Option 3: Arc-walk + fixed stride
        r3 = cluster_arc_walk(model, data, device, hits_per_cluster=8)
        evaluate_clustering(r3, "Option 3: Arc-walk (stride=8)")

        # Comparison plot
        results = [
            ("HDBSCAN (embedding)", r1),
            ("Arc-walk (stride=8)", r3),
        ]
        plot_clustering_comparison(
            results, ev_idx,
            save_path=f"../data/gnn_training/clustering_ev{ev_idx}.png"
        )

        # Embedding visualization for first event
        if ev_idx == 0:
            plot_embeddings_2d(
                r1, save_path="../data/gnn_training/embeddings_pca.png"
            )

    # Aggregate statistics over many events
    print(f"\n{'#'*60}")
    print("Aggregate over 200 events...")
    hdb_clusters, arc_clusters = [], []
    hdb_purity, arc_purity = [], []

    for i in range(min(200, len(val_data))):
        data = val_data[i]

        r1 = cluster_embedding_hdbscan(model, data, device)
        r3 = cluster_arc_walk(model, data, device, hits_per_cluster=8)

        # Proton track stats only
        proton_mask = data.track_id.cpu().numpy() != 0
        for r, cl_list, pur_list in [(r1, hdb_clusters, hdb_purity),
                                       (r3, arc_clusters, arc_purity)]:
            proton_labels = r["labels"][proton_mask]
            n_cl = len(set(proton_labels[proton_labels >= 0]))
            cl_list.append(n_cl)

            # Purity of proton clusters
            for cl in set(proton_labels[proton_labels >= 0]):
                cl_mask = r["labels"] == cl
                dominant = np.bincount(r["track_ids"][cl_mask]).max()
                pur_list.append(dominant / cl_mask.sum())

    print(f"\nProton clusters per event:")
    print(f"  HDBSCAN:  {np.mean(hdb_clusters):.1f} ± {np.std(hdb_clusters):.1f}")
    print(f"  Arc-walk: {np.mean(arc_clusters):.1f} ± {np.std(arc_clusters):.1f}")
    print(f"\nProton cluster purity:")
    print(f"  HDBSCAN:  {np.mean(hdb_purity):.4f}")
    print(f"  Arc-walk: {np.mean(arc_purity):.4f}")


if __name__ == "__main__":
    main()
