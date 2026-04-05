"""GNN edge classifier for AT-TPC track separation."""

import torch
import torch.nn as nn
from torch_geometric.nn import GravNetConv


class EdgeClassifier(nn.Module):
    """GNN-based edge classifier using GravNet layers.

    Architecture:
        node features → encoder MLP → N × GravNet layers (with residual) →
        concatenate edge endpoints → MLP → sigmoid → edge score
    """

    def __init__(self, in_channels=6, hidden=64, n_gravnet=4, gravnet_k=20,
                 gravnet_space=4):
        super().__init__()

        # Node encoder
        self.encoder = nn.Sequential(
            nn.Linear(in_channels, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
        )

        # GravNet layers with residual connections
        self.gravnet_layers = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        for _ in range(n_gravnet):
            self.gravnet_layers.append(
                GravNetConv(hidden, hidden, space_dimensions=gravnet_space,
                            propagate_dimensions=hidden, k=gravnet_k)
            )
            self.layer_norms.append(nn.LayerNorm(hidden))

        # Edge classifier: concatenate src + dst node features → score
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden // 2),
            nn.ReLU(),
            nn.Linear(hidden // 2, 1),
        )

    def forward(self, data, return_embeddings=False):
        x, edge_index = data.x, data.edge_index
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None else None

        # Encode
        h = self.encoder(x)

        # Message passing with residuals
        for gravnet, norm in zip(self.gravnet_layers, self.layer_norms):
            h_new = gravnet(h, batch)
            h = norm(h + h_new)  # residual + layer norm

        # Edge classification
        src, dst = edge_index
        edge_features = torch.cat([h[src], h[dst]], dim=-1)
        edge_scores = self.edge_mlp(edge_features).squeeze(-1)

        if return_embeddings:
            return edge_scores, h
        return edge_scores


if __name__ == "__main__":
    from dataset import ATTPCDataset

    ds = ATTPCDataset("../data/gnn_training/all_hits.csv", k=10)
    data = ds[0]

    model = EdgeClassifier(in_channels=6, hidden=64, n_gravnet=4)
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Forward pass
    with torch.no_grad():
        scores = model(data)
    print(f"Edge scores shape: {scores.shape}")
    print(f"Score range: [{scores.min():.3f}, {scores.max():.3f}]")
