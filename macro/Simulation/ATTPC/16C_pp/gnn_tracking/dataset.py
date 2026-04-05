"""AT-TPC GNN Dataset: loads hit CSV and builds per-event PyG Data objects."""

import os
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data
from torch_geometric.nn import knn_graph


def build_dataset(csv_path, k=10, reaction_only=True, cache_path=None):
    """Build list of PyG Data objects from CSV, with disk caching.

    Returns a plain list[Data] — use with DataLoader directly.
    """
    # Check cache
    if cache_path and os.path.exists(cache_path):
        print(f"Loading cached graphs from {cache_path}...")
        return torch.load(cache_path, weights_only=False)

    print(f"Loading {csv_path}...")
    df = pd.read_csv(csv_path)

    # Group by event, filter — keep original event_id
    events = []
    event_ids = []
    for eid, group in df.groupby("event_id"):
        if reaction_only and group["track_id"].nunique() < 2:
            continue
        events.append(group.reset_index(drop=True))
        event_ids.append(int(eid))

    print(f"Building graphs for {len(events)} events...")

    # Normalization stats
    all_data = pd.concat(events, ignore_index=True)
    cols = ["x", "y", "z", "charge", "r", "phi"]
    means = torch.tensor([all_data[c].mean() for c in cols], dtype=torch.float)
    stds = torch.tensor([max(all_data[c].std(), 1e-9) for c in cols], dtype=torch.float)

    data_list = []
    for i, group in enumerate(events):
        if (i + 1) % 1000 == 0:
            print(f"  {i+1}/{len(events)}...")

        x_raw = torch.tensor(group[cols].values, dtype=torch.float)
        x_norm = (x_raw - means) / stds
        pos = x_raw[:, :3]

        edge_index = knn_graph(pos, k=k, loop=False)

        track_ids = torch.tensor(group["track_id"].values, dtype=torch.long)
        src, dst = edge_index
        edge_labels = (track_ids[src] == track_ids[dst]).float()

        data = Data(
            x=x_norm, pos=pos, edge_index=edge_index, y=edge_labels,
            track_id=track_ids, num_nodes=len(group),
            event_id=event_ids[i],
        )
        data_list.append(data)

    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        print(f"Caching to {cache_path}...")
        torch.save(data_list, cache_path)

    return data_list


def train_val_split(data_list, val_fraction=0.15, seed=42):
    """Split into train and val lists."""
    n = len(data_list)
    indices = np.arange(n)
    np.random.RandomState(seed).shuffle(indices)
    n_val = int(n * val_fraction)
    val_list = [data_list[i] for i in indices[:n_val]]
    train_list = [data_list[i] for i in indices[n_val:]]
    return train_list, val_list


if __name__ == "__main__":
    data_list = build_dataset(
        "../data/gnn_training/all_hits.csv", k=10,
        cache_path="../data/gnn_training/graphs_k10.pt"
    )
    print(f"\nDataset: {len(data_list)} events")

    data = data_list[0]
    print(f"\nEvent 0:")
    print(f"  Nodes: {data.num_nodes}")
    print(f"  Edges: {data.edge_index.shape[1]}")
    print(f"  Features: {data.x.shape}")
    print(f"  True edge fraction: {data.y.mean().item():.3f}")
    print(f"  Tracks: {data.track_id.unique().tolist()}")
