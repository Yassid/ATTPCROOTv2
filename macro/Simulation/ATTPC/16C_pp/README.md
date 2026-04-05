# 16C(p,p) UKF Track Reconstruction

Two track reconstruction pipelines for the 16C(p,p) elastic scattering reaction in the AT-TPC.
Both share the same simulation and digitization steps, then diverge at track finding and clustering.

## Common Setup

```bash
# Load environment (required for all steps)
source /path/to/ATTPCROOTv2-OpenKF/build/config.sh

# Working directory
cd macro/Simulation/ATTPC/16C_pp/
```

## Step 0: Simulation and Digitization (shared)

These steps produce the common input file `data/output_digi.root`.

```bash
# 1. Simulate 10000 events (beam + reaction alternating)
root -b -q 'C16_pp_sim.C(10000)'

# 2. Digitize + PSA + Pattern Recognition (RANSAC)
root -b -q run_digi_attpc.C
```

Output: `data/output_digi.root` (contains AtPatternEvent with RANSAC tracks)

---

## Method 1: RANSAC + ClusterizeSmooth3D + UKF

The standard pipeline. RANSAC separates tracks, ClusterizeSmooth3D clusters
hits within each track (adaptive radius based on KE estimate), UKF fits.

### Run

```bash
root -b -q 'run_ukf_only.C(10000)'
```

### Output

- `data/output_ukf_only.root` (AtTrackingEvent with fitted tracks)
- Copy to: `data/smooth3d/output_ukf_only.root`

### Visualize

```bash
# Interactive display with control panel, live refitting, diagnostics
root -l run_ukf_display.C

# Navigate: DrawEvent(n), DrawEvent(-1) next, DrawEvent(-2) prev
```

### Analyze

```bash
# Kinematic curves (KE vs theta)
root -b -q show_kinematics.C

# Chi2 quality cuts
root -b -q chi2_cuts_kinematics.C

# Chi2 vs resolution correlation
root -b -q chi2_vs_resolution.C

# Vertex distance vs error
root -b -q vtxR_vs_error2.C
```

### Performance (10000 events)

- Efficiency: 95.6%
- KE bias: -0.2%, RMS: 5.1%
- Theta bias: 0.21 deg, RMS: 0.73 deg
- Tracks fitted: ~4863 (per 5002 reaction events)

---

## Method 2: GNN + ClusterizeSmooth3D + UKF

A GNN (GravNet edge classifier) replaces RANSAC for track separation.
The GNN achieves 99.98% track purity vs RANSAC. ClusterizeSmooth3D
still handles clustering within each separated track.

### Prerequisites

```bash
# Python environment (one-time setup)
python3 -m venv ~/gnn_env
source ~/gnn_env/bin/activate
pip install torch torch_geometric pandas scikit-learn hdbscan matplotlib
```

### Step 1: Extract training data (one-time)

```bash
root -b -q 'extract_gnn_data.C(10000)'
```

Output: `data/gnn_training/all_hits.csv`

### Step 2: Train the GNN (one-time)

```bash
source ~/gnn_env/bin/activate
cd gnn_tracking
python train.py
cd ..
```

Output: `gnn_tracking/best_model.pt`

### Step 3: Export GNN-separated hits

```bash
source ~/gnn_env/bin/activate
cd gnn_tracking
python gnn_hits_export.py --n-events 5002
cd ..
```

Output:
- `data/gnn_training/gnn_separated_hits.csv` (raw hits per GNN track)
- `data/gnn_training/gnn_separated_tracks.csv` (per-track geometry)

### Step 4: Build AtPatternEvent ROOT file

```bash
root -b -q 'write_gnn_digi.C(5002)'
```

Output: `data/output_gnn_digi.root`

This applies observable-based beam rejection (theta_lab 20-90 deg, radius > 30 mm)
and pre-clusters with ClusterizeSmooth3D for display compatibility.

### Step 5: Run UKF fitter

```bash
root -b -q 'run_ukf_gnn_smooth.C(5002)'
```

Output:
- `data/output_ukf_gnn_smooth.root` (AtTrackingEvent with fitted tracks)
- Copy to: `data/gnn_smooth/output_ukf_gnn_smooth.root`

### Visualize

```bash
# Interactive display (same GUI as Method 1)
root -l run_ukf_display_gnn.C

# Navigate: DrawEvent(n), DrawEvent(-1) next, DrawEvent(-2) prev
```

### Performance (5002 events)

- Tracks fitted: ~4997 (+2.8% vs Method 1)
- KE distribution matches Method 1
- Theta: 59.9 deg (vs 59.5 deg for Method 1)
- 5% more proton tracks recovered

---

## Comparison

```bash
# Side-by-side kinematic plots (KE vs theta, distributions)
# Reads both output_ukf_only.root and output_ukf_gnn_smooth.root
root -b -q compare_gnn_vs_smooth.C
```

Output: `data/gnn_vs_smooth_comparison.png`

A detailed kinematic comparison is also at: `data/gnn_smooth/kinematics_comparison.png`

---

## Directory Structure

```
16C_pp/
├── C16_pp_sim.C                 # Simulation macro
├── run_digi_attpc.C             # Digitization + PSA + RANSAC
│
├── run_ukf_only.C               # Method 1: RANSAC + Smooth3D + UKF
├── run_ukf_display.C            # Method 1: interactive display
├── display_ukf.C                # Method 1: simple display
│
├── extract_gnn_data.C           # Extract hits to CSV for GNN training
├── write_gnn_digi.C             # Build AtPatternEvent from GNN output
├── run_ukf_gnn_smooth.C         # Method 2: GNN + Smooth3D + UKF
├── run_ukf_display_gnn.C        # Method 2: interactive display
│
├── compare_gnn_vs_smooth.C      # Compare both methods
├── show_kinematics.C            # Kinematic analysis
├── chi2_cuts_kinematics.C       # Chi2 quality cuts
├── chi2_vs_resolution.C         # Chi2 vs resolution
│
├── gnn_tracking/                # GNN Python code
│   ├── model.py                 # GravNet edge classifier (82K params)
│   ├── dataset.py               # PyG data loading + k-NN graph
│   ├── train.py                 # Training script (30 epochs)
│   ├── analyze.py               # GNN performance analysis
│   ├── cluster.py               # Arc-walk clustering (experimental)
│   ├── gnn_hits_export.py       # Export GNN-separated hits to CSV
│   ├── gnn_cluster_export.py    # Export arc-walk clusters (experimental)
│   ├── eval_kinematics.py       # Kinematic quality evaluation
│   └── best_model.pt            # Trained model weights
│
└── data/
    ├── output_digi.root         # Shared: digitized events
    ├── attpcsim.root            # Shared: MC truth
    ├── output_gnn_digi.root     # GNN-separated tracks as AtPatternEvent
    │
    ├── smooth3d/                # Method 1 results
    │   ├── output_ukf_only.root
    │   └── *.png                # Analysis plots
    │
    ├── gnn_smooth/              # Method 2 results
    │   ├── output_ukf_gnn_smooth.root
    │   ├── output_gnn_digi.root
    │   └── kinematics_comparison.png
    │
    ├── gnn_arcwalk/             # Experimental: GNN + arc-walk (no Smooth3D)
    │   └── output_ukf_gnn.root
    │
    └── gnn_training/            # GNN training data and exports
        ├── all_hits.csv         # Hit data for training
        ├── graphs_k10.pt       # Cached k-NN graphs
        ├── gnn_separated_hits.csv
        ├── gnn_separated_tracks.csv
        └── *.png                # Clustering/embedding plots
```

## GNN Architecture

- Model: GravNet edge classifier with 4 GravNet layers, 64-dim hidden, 82K parameters
- Input: 6 features per hit (x, y, z, charge, r, phi), normalized
- Graph: k-NN with k=10 neighbors
- Output: edge score (same-track probability)
- Training: BCE loss with class weighting, 30 epochs, Adam optimizer
- Track separation: threshold edge scores > 0.5, connected components via union-find
- Performance: 98.7% edge accuracy, 99.99% precision, 98.7% recall
- Proton efficiency: 98.8%, purity: 99.98%
