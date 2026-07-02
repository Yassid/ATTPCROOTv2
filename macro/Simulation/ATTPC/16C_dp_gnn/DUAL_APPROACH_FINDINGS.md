# Track separation — dual-approach exploration (autonomous session)

Goal: separate overlapping AT-TPC tracks robustly, favouring **over-segmentation** (stitchable)
over merges. Explored several "connect-the-dots" algorithms and a dual splitter+refiner.

## Methods built (all offline Python, no training, ~0.3 s/event)
- `dircluster.py` — **local direction-continuity clustering** (splitter). Per hit: tangent from
  PCA of k=12 3D neighbours. Connect neighbour pairs only if both tangents align with the segment
  AND with each other AND charge (dE/dx) ratio > QRATIO. Connected components → tracks.
- `trackfollow.py` — greedy directional track-following (cellular-automaton style).
- `stitch.py` — `stitch()` merges fragment endpoints on a common trajectory + dE/dx band;
  `attach()` grows track cores by absorbing gray hits.
- `quant_eval.py`, `quant_dircluster.py`, `diag_particle.py` — sim-truth quantifiers.

## Key results (sim truth, 250 events, merge = 2 real tracks fused in one cluster)

**dE/dx is the decisive lever** (direction alone merges 22% of events; charge continuity fixes it):
| QRATIO (dE/dx cut) | merge % | proton recovery |
|---|---|---|
| 0.0 (off) | 22.3% | — |
| 0.4 (recommended) | **0.4%** | **92.8%** |

**Per-particle recovery at QRATIO=0.4** (this is the headline):
| particle | hits | clustered % |
|---|---|---|
| proton | 96,554 | **92.8%** |
| noise | 36,642 | 0.1% (correctly rejected) |
| 17C recoil | 4,789 | 23.4% (hard: 23 hits, dense, near central hole) |

## Conclusions
1. **The direction+dE/dx splitter (`dircluster` @ QRATIO=0.4) is the method to use.** It recovers the
   physics-relevant **proton at 93%** with **0.4% merges** and **~0% noise contamination**, captures
   spirals whole, separates vertex-sharing crossings (evt 1198/1718 that the GNN+DBSCAN merged), and
   over-segments only the dense vertex/recoil region.
2. **Greedy track-following is worse** — it fragments smooth curves (connected-components keeps them whole).
3. **The dual splitter→(stitch/attach) did NOT help**: the residual loss is not fragmented protons but
   the intrinsically hard **17C recoil** (short, dense, near the hole). `stitch()` no-ops; `attach()`
   trades merges (0.4%→11.6%) for +3% recovery. So a second *general* algorithm isn't the fix.
4. If 17C recovery matters, it needs a **dedicated dense-vertex/short-recoil handler**, not a
   general refiner — but for (d,p) physics the proton is the ejectile of interest, so 93% is strong.

## Recommended pipeline
`dircluster.cluster(P, Q, qratio=0.4, min_hits=4)` → per-event track labels → downstream stitch/fit.
Galleries: `diagnostics/dircluster_gallery.png`, `dircluster_events.png`.

## Compared against the earlier GNN
Sim-trained DGCNN embedding + DBSCAN (`embed_sim_noisy.pt`) MERGED vertex-sharing tracks and needed
training + a domain-gap-closing noise model. The geometric direction+dE/dx clusterer needs neither and
beats it on the merge failure mode. The noise model work (gain fix + AtPulse additive noise) remains
valuable for making *simulated* training data realistic, if a learned method is revisited.
