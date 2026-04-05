"""Training script for AT-TPC GNN edge classifier."""

import time
import numpy as np
import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

from dataset import build_dataset, train_val_split
from model import EdgeClassifier


def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    n_batches = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        scores = model(batch)
        loss = criterion(scores, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
        n_batches += 1
    return total_loss / n_batches


@torch.no_grad()
def evaluate(model, loader, device, threshold=0.5):
    model.eval()
    all_scores, all_labels = [], []
    for batch in loader:
        batch = batch.to(device)
        scores = model(batch)
        all_scores.append(torch.sigmoid(scores).cpu())
        all_labels.append(batch.y.cpu())

    scores = torch.cat(all_scores).numpy()
    labels = torch.cat(all_labels).numpy()
    preds = (scores > threshold).astype(int)

    return {
        "accuracy": accuracy_score(labels, preds),
        "precision": precision_score(labels, preds, zero_division=0),
        "recall": recall_score(labels, preds, zero_division=0),
        "f1": f1_score(labels, preds, zero_division=0),
    }


@torch.no_grad()
def evaluate_tracking(model, data_list, device, threshold=0.5):
    """Tracking metrics: efficiency and purity via connected components."""
    model.eval()
    efficiencies, purities = [], []

    for data in data_list:
        data = data.to(device)
        scores = torch.sigmoid(model(data)).cpu()
        edge_index = data.edge_index.cpu()
        track_ids = data.track_id.cpu()

        mask = scores > threshold
        kept = edge_index[:, mask]
        if kept.shape[1] == 0:
            continue

        # Union-find
        n = data.num_nodes
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for j in range(kept.shape[1]):
            a, b = kept[0, j].item(), kept[1, j].item()
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        clusters = {}
        for node in range(n):
            clusters.setdefault(find(node), []).append(node)

        truth_tracks = {}
        for node in range(n):
            truth_tracks.setdefault(track_ids[node].item(), set()).add(node)

        for tid, truth_nodes in truth_tracks.items():
            best_overlap, best_size = 0, 0
            for cl_nodes in clusters.values():
                overlap = len(truth_nodes & set(cl_nodes))
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_size = len(cl_nodes)
            if len(truth_nodes) > 0:
                efficiencies.append(best_overlap / len(truth_nodes))
                purities.append(best_overlap / best_size if best_size > 0 else 0)

    return {
        "track_efficiency": np.mean(efficiencies) if efficiencies else 0,
        "track_purity": np.mean(purities) if purities else 0,
    }


def main():
    csv_path = "../data/gnn_training/all_hits.csv"
    cache_path = "../data/gnn_training/graphs_k10.pt"
    k = 10
    hidden = 64
    n_gravnet = 4
    lr = 1e-3
    epochs = 30
    batch_size = 32
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load data
    data_list = build_dataset(csv_path, k=k, cache_path=cache_path)
    train_data, val_data = train_val_split(data_list)
    print(f"Train: {len(train_data)}, Val: {len(val_data)}")

    train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=batch_size, shuffle=False)

    # Model
    model = EdgeClassifier(in_channels=6, hidden=hidden, n_gravnet=n_gravnet).to(device)
    params = sum(p.numel() for p in model.parameters())
    print(f"Parameters: {params:,}")

    # Class weighting
    pos_frac = data_list[0].y.mean().item()
    pos_weight = torch.tensor([(1 - pos_frac) / pos_frac], device=device)
    print(f"Edge balance: {pos_frac:.3f} positive, pos_weight={pos_weight.item():.2f}")

    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)

    best_f1 = 0
    for epoch in range(1, epochs + 1):
        t0 = time.time()
        loss = train_epoch(model, train_loader, optimizer, criterion, device)
        metrics = evaluate(model, val_loader, device)
        dt = time.time() - t0
        scheduler.step(loss)

        print(f"Ep {epoch:2d}/{epochs} | loss={loss:.4f} | "
              f"acc={metrics['accuracy']:.4f} prec={metrics['precision']:.4f} "
              f"rec={metrics['recall']:.4f} f1={metrics['f1']:.4f} | {dt:.1f}s",
              flush=True)

        if metrics["f1"] > best_f1:
            best_f1 = metrics["f1"]
            torch.save(model.state_dict(), "best_model.pt")

    # Final
    print("\n=== Final Evaluation ===")
    model.load_state_dict(torch.load("best_model.pt", weights_only=True))
    metrics = evaluate(model, val_loader, device)
    print(f"Edge: acc={metrics['accuracy']:.4f} prec={metrics['precision']:.4f} "
          f"rec={metrics['recall']:.4f} f1={metrics['f1']:.4f}")

    print("\nTracking evaluation (200 events)...")
    track = evaluate_tracking(model, val_data[:200], device)
    print(f"Track efficiency: {track['track_efficiency']:.4f}")
    print(f"Track purity: {track['track_purity']:.4f}")


if __name__ == "__main__":
    main()
