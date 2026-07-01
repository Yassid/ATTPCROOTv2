# Clustering / track-finding algorithms for AT-TPC — CTD & TrackML survey
Deep-research result (2026-06-27). 108 claims, all survived 3-vote adversarial verification (0 refuted).

## 1. Learned-embedding + density clustering ("embed then cluster")
- **Object Condensation** (Kieseler, EPJC 2020, arXiv:2002.03605) — grid-free, one-stage
  multi-object reconstruction. **No assumptions on object size/density/number.** Learns a
  latent space where same-object points condense around a high-confidence "condensation
  point"; clustering is learnable + local. Works natively on **graphs and point clouds**.
  Code: **cms-pepr/HGCalML** (GravNet + OC, TensorFlow). → strongest conceptual fit for
  vertex-sharing/overlapping tracks.
- **TrackML accuracy-phase clustering** — engineer/learn per-hit coordinates, then **DBSCAN**
  (CPMP "DBSCAN-for-ever", outrunner, top-quarks). Direct upgrade path from Spyral's HDBSCAN.
  Code: github liaopeiyuan/TrackML-Particle-Tracking-Challenge, nicolefinnie/kaggle-trackml.
- **Exa.TrkX graph construction** uses a **metric-learning embedding** to build the graph, then
  GNN, then DBSCAN/connected-components (see #2).

## 2. GNN track finding (Exa.TrkX / GNN4ITk lineage)
- **Exa.TrkX** (EPJC 81, 2021, arXiv:2103.06995): metric-learning embedding → GNN edge
  classification → track building. Quality comparable to production trackers; **compute scales
  ~linearly with #hits** (good for AT-TPC's hundreds–thousands of hits). Validated on TrackML,
  **generalizes to DUNE & CMS geometries**. Code: **github.com/exatrkx** (acorn/commonframework,
  C++ inference). Pipeline repo: HSF-reco-and-software-triggers/Tracking-ML-Exa.TrkX.
- **GNN4ITk** (ATLAS): 3 stages — graph construction (metric learning OR Module Map) → GNN edge
  labeling → **connected components + walk-through**. Fake & duplicate rates O(10⁻³). Adding
  cluster-position features raised edge purity 70%→>90%. Code: **GNN4ITkTeam CommonFramework**.
- **Edge-classifying Interaction Networks** (Comput. Softw. Big Sci. 5, 2021).
- **Hierarchical GNN** (arXiv:2303.01640); **EggNet** evolving graph-attention (arXiv:2407.13925);
  **Transformers for tracking** (arXiv:2411.07149); **HEPT** — LSH, near-linear, point-cloud
  (CTD2025). **Large-radius / non-pointing tracks** with Exa.TrkX (arXiv:2203.08800).
- **TPC-specific GNN precedents:** DUNE LArTPC end-to-end DL chain (arXiv:2102.01033); CTD2023
  "Multipurpose GNN for LArTPC" (Cerati, Fermilab) — closest detector analogue; CTD2023
  "Combined track finding with GNN & CKF" (Huth) — hybrid ML+Kalman.
- **CTD2025 directly on our problems:** "double metric learning / spacepoint-doublet embedding
  that explicitly handles **shared spacepoints between tracks**" (Jay Chan, LBNL); "GNN4ITk
  adapted to **non-helical** trajectories" (Condren, UC Irvine); GNN for multi-TPC TDIS@JLab
  (Shujie Li, LBNL).

## 3. Classical methods for curved / looping tracks
- **ConformalTracking** (iLCSoft) — **conformal mapping + cellular automaton**. C++, GPLv3,
  RELEASED. Tunable MaxCellAngle / MaxCellAngleRZ / Chi2Cut; vertex-constrained or not.
  Geometry-agnostic (built for CLIC; needs porting). Best ready-to-try classical curved-track code.
- **4D Cellular Automaton** track finder, CBM (CTD2016, EPJ Web Conf.).
- **Belle II** finder: CA + Legendre transform + **CKF**, efficient down to **50 MeV/c** strongly
  curved tracks; FastBDT for background (arXiv:2003.12466). Proof CA+CKF handles low-p loopers.
- **ACTS** Combinatorial Kalman Filter (released, well-documented). **PANDA TPC** TDR
  (arXiv:1207.0013) — TPC-specific pattern recognition.

## 4. TrackML challenge — what actually won
- Task = cluster ~100k 3D points/event into ~10k helical tracks; score = fraction of points
  correctly associated. Accuracy phase (Kaggle, 2018) + Throughput phase (Codalab, 2-core/4GB).
  Papers: arXiv:1904.06778 (accuracy), EPJ 10.1007/s41781-023-00094-w (throughput). Dataset:
  Kaggle, Codalab, CERN Open Data. Site: sites.google.com/site/trackmlparticle.
- **Throughput winners were CLASSICAL combinatorial, not DL:** 1st "Mikado" (Gorbunov, 94.4%,
  0.56 s/evt) — **find easiest high-p (nearly straight) tracks first, remove their hits, recurse**
  (directly portable iterative peeling). 3rd (Kunze) — ANN pattern-rec + triplet prolongation +
  outlier density. All winners released code.
- **Accuracy winners** leaned on embed-then-DBSCAN (see #1).

## 5. AT-TPC / nuclear-physics ML (the gap = opportunity)
- Existing AT-TPC ML is mostly **classification / vertex regression, NOT track-finding/clustering**:
  - Kuchera et al. MSU/NSCL, "ML Methods for Track Classification in the AT-TPC" (arXiv:1810.10350)
    — CNN image-based.
  - "ML for ¹²C event classification & reconstruction in AT-TPC" (arXiv:2304.13233, NIM A 2023) —
    VGG/ResNet; ResNet-34 0.99 class; regression σ_E<77 keV, σ_θ<0.1 rad.
  - PointNet++ rare-event classification in AT-TPC (ResearchGate 385334464).
  - MATE-TPC ¹²C+¹²C event-class + vertex (arXiv:2605.28296).
- **PEAR** (arXiv:2502.12169): **PointNet ensemble** vertex reconstruction directly on TPC
  spacepoint clouds (ALPHA-g radial TPC) — **outperforms trajectory fitting on all metrics,
  recovers vertices where standard fails.** Direct precedent that point-cloud DL beats fitting
  on a TPC.
- **Sim-to-data caveat:** "Unpaired Translation of Point Clouds for Modeling Detector Response"
  (AT-TPC, NeurIPS ML4PS 2024) — addresses exactly the sim→data domain gap we must watch.

## 6. CTD proceedings — where to browse
- Landing: **epj-conferences.org/ctd** (EPJ Web of Conferences).
- **CTD 2023** Toulouse (10–13 Oct 2023): indico.cern.ch/event/1252748 (timetable + proceedings).
- **CTD 2025** Univ. of Tokyo (10–14 Nov 2025): indico.cern.ch/event/1499357.

## Shortlist for AT-TPC (ranked)
1. **Object Condensation (GravNet, HGCalML)** — point-cloud native, grid-free, handles
   overlapping/vertex-sharing & variable #tracks by construction. Train on sim. Best "different
   clustering algorithm" bet.
2. **Exa.TrkX / GNN4ITk (embed → GNN → connected-components)** — most mature, released code,
   proven on LArTPC/DUNE; CTD2025 double-metric-learning variant explicitly handles shared
   spacepoints; non-helical variant exists. Heavier infra.
3. **Embed-then-DBSCAN (TrackML accuracy)** — lightest lift, direct upgrade from Spyral HDBSCAN:
   learn an embedding on sim, cluster in it.
4. **ConformalTracking (conformal+CA)** — released C++, no training, strong on curved tracks;
   good classical baseline / fallback.
All four train or run on labeled simulations (we can generate per-hit truth). Watch the sim→data
domain gap (point-cloud translation paper). PEAR is the precedent that this works on a TPC.
