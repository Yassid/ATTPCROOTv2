# AT-TPC Track Clustering with Physics-Informed GNNs — Implementation Plan

## 1. Project Overview

Build a GNN-based pipeline to cluster (separate) particle tracks in the Active Target Time Projection Chamber (AT-TPC) using PyTorch Geometric. The pipeline takes 3D point clouds of ionization hits and assigns each hit to a track.

### Why GNNs for AT-TPC?
- AT-TPC data is naturally a 3D point cloud → graphs are the ideal representation
- GNNs handle variable-size inputs, irregular geometry, and curved tracks
- They can be made physics-aware (energy loss, conservation laws, B-field geometry)
- Proven in similar detectors: DUNE LAr-TPC, GlueX drift chamber, Belle II CDC, LHC trackers

---

## 2. Pipeline Architecture (4 Stages)

The pipeline follows the Exa.TrkX / ACORN pattern:

```
Raw Hits → [Stage 1: Data Prep] → [Stage 2: Graph Construction] → [Stage 3: GNN Edge Classification] → [Stage 4: Track Building]
```

### Stage 1 — Data Preparation
- **Input**: AT-TPC spacepoints (x, y from pad plane; z from drift time; charge amplitude)
- **Training data source**: Geant4/Garfield++ simulation with ground-truth track-to-hit assignments
- **Output**: PyG `Data` objects with:
  - `x`: Node feature matrix [N_hits, N_features]
  - `edge_index`: Edge connectivity [2, N_edges]
  - `edge_attr`: Edge features (optional) [N_edges, N_edge_features]
  - `y`: Edge labels (1 = same track, 0 = different track) [N_edges]
  - `track_id`: Per-node ground-truth track assignment [N_hits]

#### Suggested Node Features (per hit):
| Feature | Description | Physics motivation |
|---------|-------------|-------------------|
| x, y, z | 3D position | Spatial location |
| charge | Pad amplitude | Proportional to dE/dx |
| r | Distance to beam axis | Separates beam from reaction products |
| phi | Azimuthal angle | Useful in cylindrical geometry |
| local_charge_density | Charge in neighborhood | Proxy for local dE/dx |
| local_curvature | Curvature from neighbor triplets | Encodes momentum (in B-field) |
| n_neighbors | Number of nearby hits | Track density indicator |

#### Suggested Edge Features (per edge):
| Feature | Description |
|---------|-------------|
| dr | Euclidean distance between hits |
| dz | Difference in z (drift direction) |
| dphi | Difference in azimuthal angle |
| d_charge | Difference in charge |
| direction_cosines | Unit vector between hits (3 components) |

### Stage 2 — Graph Construction

**Option A — Geometric k-NN graph (start here):**
```python
from torch_geometric.nn import knn_graph
edge_index = knn_graph(pos, k=15, loop=False)
```
- Use k = 10–20 as starting point
- Consider using cylindrical coordinates (r, phi, z) instead of (x, y, z) if detector has cylindrical symmetry
- If B-field is present, consider helical distance metric

**Option B — Radius graph:**
```python
from torch_geometric.nn import radius_graph
edge_index = radius_graph(pos, r=max_radius, loop=False)
```
- Better control over local connectivity
- Need to tune radius per detector geometry

**Option C — Metric learning (advanced, Exa.TrkX style):**
- Train an MLP to embed hits into latent space
- Build k-NN graph in latent space
- Hits from same track should be close in learned space

### Stage 3 — GNN Edge Classification

**Recommended architecture for first prototype:**

```
Input features → [Encoder MLP] → [N × Message Passing Layers with residual connections] → [Edge Classifier MLP] → edge scores
```

**GNN layer options (all available in PyG):**

| Layer | PyG Class | Pros | Best for |
|-------|-----------|------|----------|
| GravNet | `GravNetConv` | Dynamic graph, distance-weighted, designed for particle detectors | First choice for AT-TPC |
| Interaction Network | Custom via `MessagePassing` | Explicit edge+node updates, used by Exa.TrkX | High accuracy |
| EdgeConv (DGCNN) | `DynamicEdgeConv` | Recomputes neighbors per layer | Multi-scale tracks |
| GAT | `GATConv` | Attention-weighted messages | Variable importance neighbors |
| GraphSAGE | `SAGEConv` | Scalable sampling-based | Very large events |

**Architecture details:**
- Encoder MLP: maps raw features to hidden dimension (e.g., 64 or 128)
- Message passing: 4–8 layers with skip/residual connections
- Edge classifier: concatenate endpoint node features → MLP → sigmoid → score ∈ [0, 1]
- Hidden dimension: 64–128 (AT-TPC events are smaller than LHC, so you don't need huge models)

**Loss function:**
```python
loss = BCE_loss(edge_scores, edge_labels)
```
- Use class weighting or focal loss since fake edges >> true edges
- Optionally add physics-informed loss terms (see Section 3 below)

### Stage 4 — Track Building (Post-Processing)
1. Threshold edge scores (e.g., score > 0.5)
2. Extract connected components → each component = one track candidate
3. Filter: remove clusters with < min_hits
4. Optional: fit each cluster to physical model (line/helix), reject poor fits
5. Optional: split clusters with multiple tracks using fit residuals

---

## 3. Physics-Informed Enhancements

### 3.1 Physics in Graph Construction
- **Helical distance metric**: If B-field is present, replace Euclidean distance with helical arc length when building k-NN graph
- **Beam-direction bias**: Weight edges along z (beam axis) more heavily
- **Energy-loss aware edges**: Filter edges where charge gradient is unphysical (e.g., charge increases along forward direction for a stopping particle)

### 3.2 Physics in Node/Edge Features
- **dE/dx features**: Local energy loss is characteristic of particle species (Bethe-Bloch)
- **Curvature features**: Local curvature from hit triplets encodes momentum in B-field
- **Bragg peak indicators**: Charge profile near track endpoints differs from mid-track

### 3.3 Physics-Informed Loss Functions
```python
L_total = L_BCE + λ₁ * L_smoothness + λ₂ * L_momentum + λ₃ * L_energy + λ₄ * L_vertex
```

| Loss term | What it penalizes |
|-----------|-------------------|
| `L_smoothness` | Abrupt direction changes along reconstructed tracks |
| `L_momentum` | Violations of transverse momentum conservation at vertex |
| `L_energy` | Violations of energy conservation given known beam energy |
| `L_vertex` | Extrapolated tracks that don't converge to a common vertex |
| `L_bethe_bloch` | Energy loss profiles inconsistent with any known particle species |

### 3.4 Physics in Architecture
- **E(n)-Equivariant GNNs (EGNN)**: Guarantee rotational/translational symmetry → better data efficiency
- **Attention kernels**: Weight messages by physics-motivated functions (e.g., Gaussian of helical distance)
- **Hamiltonian / Lagrangian Neural Networks**: Learn dynamics that conserve energy by construction

### 3.5 Hybrid GNN + Physics Solver
1. GNN predicts edge scores
2. Extract clusters
3. Fit each cluster with differentiable ODE solver (`torchdiffeq`) for charged particle in gas + B-field
4. Backpropagate through physics simulation into GNN weights
5. End-to-end training: GNN learns to produce clusters that are physically consistent

---

## 4. Required Software Stack

### Core dependencies:
```
torch >= 2.0
torch-geometric >= 2.4
torch-cluster
torch-scatter
torch-sparse
pytorch-lightning >= 2.0
numpy
scipy
h5py or uproot (depending on data format)
scikit-learn (for evaluation metrics)
matplotlib (for visualization)
```

### Optional (for advanced features):
```
e3nn                  # E(3)-equivariant neural networks
torchdiffeq           # Differentiable ODE solvers (for Hamiltonian/Neural ODE approaches)
torch-geometric-temporal  # If temporal features are relevant
wandb or tensorboard  # Experiment tracking
optuna                # Hyperparameter optimization
```

---

## 5. Suggested Code Structure

```
attpc-gnn-tracking/
├── config/
│   ├── default.yaml          # Default hyperparameters
│   └── experiment/           # Experiment-specific configs
├── data/
│   ├── dataset.py            # PyG Dataset class for AT-TPC events
│   ├── preprocessing.py      # Hit cleaning, feature engineering
│   ├── graph_construction.py # k-NN, radius, or metric-learning graph builders
│   └── transforms.py         # PyG transforms (normalization, augmentation)
├── models/
│   ├── gnn_edge_classifier.py    # Main GNN model
│   ├── layers.py                 # Custom message-passing layers
│   ├── encoder.py                # Feature encoder MLPs
│   ├── edge_classifier.py        # Edge classification head
│   └── metric_learning.py        # Optional: embedding network for graph construction
├── losses/
│   ├── bce_loss.py               # Weighted BCE / focal loss
│   └── physics_losses.py         # Smoothness, conservation, Bethe-Bloch losses
├── tracking/
│   ├── track_builder.py          # Connected components + post-processing
│   ├── track_fitting.py          # Physical track fitting (line/helix)
│   └── vertex_reconstruction.py  # Vertex finding from fitted tracks
├── evaluation/
│   ├── metrics.py                # Tracking efficiency, purity, fake rate
│   └── visualization.py         # 3D event displays with colored tracks
├── training/
│   ├── lightning_module.py       # PyTorch Lightning training module
│   └── callbacks.py              # Custom callbacks (LR scheduling, checkpointing)
├── scripts/
│   ├── train.py                  # Training entry point
│   ├── evaluate.py               # Evaluation script
│   └── inference.py              # Run on new data
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_graph_visualization.ipynb
│   └── 03_results_analysis.ipynb
├── requirements.txt
└── README.md
```

---

## 6. Implementation Roadmap

### Phase 1 — Minimal Viable Pipeline (start here)
- [ ] Data loading: read AT-TPC simulation, extract hits with track labels
- [ ] Simple feature engineering: (x, y, z, charge, r, phi)
- [ ] k-NN graph construction in Euclidean space
- [ ] GNN model with GravNetConv or InteractionNetwork layers
- [ ] BCE loss with class weighting
- [ ] Connected components track building
- [ ] Basic evaluation: tracking efficiency, purity

### Phase 2 — Physics Enhancements
- [ ] Add physics-motivated features (local dE/dx, curvature, Bragg peak)
- [ ] Implement physics-informed loss terms (smoothness, conservation)
- [ ] Switch to cylindrical or helical distance for graph construction
- [ ] Add edge features to message passing

### Phase 3 — Advanced Methods
- [ ] Metric learning for graph construction
- [ ] E(n)-equivariant architecture
- [ ] Differentiable track fitting in the loss
- [ ] Hyperparameter optimization (Optuna)
- [ ] Uncertainty quantification

### Phase 4 — Production
- [ ] Training on full simulation dataset
- [ ] Validation on real AT-TPC data
- [ ] Integration with ATTPCROOT or experiment framework
- [ ] Performance optimization (ONNX export, batched inference)

---

## 7. Key References

### Frameworks:
- **Exa.TrkX**: https://exatrkx.github.io/ (GNN tracking pipeline, demonstrated on DUNE LAr-TPC)
- **ACORN**: Evolution of Exa.TrkX with uncertainty quantification
- **Ariadne**: https://github.com/t3hseus/ariadne (PyTorch library for DL tracking)
- **ACTS**: https://acts.readthedocs.io/ (common tracking software with Exa.TrkX plugin)
- **PyTorch Geometric**: https://pytorch-geometric.readthedocs.io/

### Key papers:
- Exa.TrkX pipeline performance: Eur. Phys. J. C 81, 876 (2021), arXiv:2103.06995
- GravNet: arXiv:1902.07987 (distance-weighted GNNs for particle detectors)
- GlueX GNN tracking: arXiv:2505.22504 (very relevant — solenoid + drift chamber)
- Belle II GNN tracking: arXiv (2025, end-to-end multi-track with GNN)
- AT-TPC ML classification: Kuchera et al., arXiv:1810.10350
- AT-TPC deep learning inpainting: APS DNP 2019
- PEAR (PointNet for TPC vertex reconstruction): ALPHA-g experiment
- Ariadne library: arXiv:2109.08982
- Hierarchical GNN for tracking: arXiv:2303.01640
- EGNN (E(n)-equivariant): Satorras et al., arXiv:2102.09844
- Hamiltonian Neural Networks: Greydanus et al., arXiv:1906.01563

### AT-TPC specific:
- ATTPCROOT: https://github.com/ATTPC (analysis framework from FRIB/MSU)
- Simulation tools: Geant4 + Garfield++ for detector response

---

## 8. Questions to Resolve Before Coding

1. **Data format**: What format are your AT-TPC hits in? (HDF5, ROOT TTree, numpy, CSV?)
2. **Magnetic field**: Does your AT-TPC operate inside a solenoid? If so, what field strength?
3. **Event complexity**: How many hits per event? How many tracks per event?
4. **Simulation availability**: Do you have Geant4/Garfield++ simulation with truth labels?
5. **Physics reactions**: What reactions are you studying? (affects track topologies)
6. **Computing resources**: GPU availability? (single GPU is fine for prototyping)
7. **Integration target**: Do you need to integrate with ATTPCROOT or another framework?
