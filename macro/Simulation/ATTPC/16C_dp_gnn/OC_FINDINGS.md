# Object Condensation for AT-TPC track separation — exploration

Ref: Lieret & DeZoort, "An Object Condensation Pipeline for Charged Particle Tracking at the
HL-LHC", arXiv:2309.16754 (based on Kieseler 2020).

## Algorithm
Per hit i, a GNN predicts a condensation strength beta_i in (0,1) and clustering coordinates
c_i in a learned latent space. Charge q_i = arctanh^2(beta_i) + q_min. Per track t, the
condensation point (CP) is the max-beta hit; charge q^(t), position c^(t).
- Potential loss: L_V = (1/N) sum_i q_i sum_t [ delta(l_i=t) q^(t) ||c^(t)-c_i||^2
  + s_rep (1-delta) q^(t) max(0, 1-||c^(t)-c_i||) ]  (attract own hits to CP, repel others).
- Beta loss: L_beta = (1/Nt) sum_t (1-beta^(t)) + s_B * <beta of noise hits>.
- L = L_V + s_beta L_beta.  Paper: s_rep=0.6, s_beta=0.004, q_min=0.34, s_B=0.09.
- Inference: DBSCAN on the coords c_i (paper eps=0.279, min_samples=1). (Greedy CP-seeded
  assignment is the alternative.)

## Implementation (this repo)
- train_oc.py : DGCNN (DynamicEdgeConv) backbone -> coord head (latent d=3) + beta head. Full
  OC loss. Trained on noisy-sim truth (sim_noisy.parquet: proton/17C tracks + noise=-1). GPU.
- eval_compare.py : head-to-head OC vs dircluster on the SAME held-out sim events.
- oc_gallery.py : apply sim-trained OC to real events, render DBSCAN-on-coords clusters.

## Results (held-out sim, 131 events)
1. BETA COLLAPSE with paper's s_beta=0.004: beta stayed ~0.05 (signal) / 0.016 (noise), below
   any useful seeding threshold -> greedy-beta inference found ~0 clusters (5% recovery). The
   s_beta=0.004 was tuned for their interaction-network/trackML setup, too weak here.
2. But the LATENT COORDS learned excellent separation: intra-track spread 0.14 vs inter-track
   center distance 0.94 (~7x). DBSCAN on coords recovers well.
3. Retrain with s_beta=0.1 (beta STILL collapses - stubborn: raising beta raises q which raises
   L_V, so pressure keeps beta low; DBSCAN-on-coords doesn't need beta anyway) IMPROVED the
   COORDS. OC DBSCAN eps-sweep vs dircluster(q=0.65,no override), same 131 held-out sim events:
     eps    OC:merge% recov% proton% 17C%    (dircluster: merge1.5 recov39.3 proton95.3 17C48.8)
     0.05     1.5      65.3   96.5   51.3
     0.06     3.8      68.9   97.1   60.5
     0.07     7.6      74.4   97.6   68.5
   => OC DOMINATES: at eps=0.05, SAME merge as dircluster (1.5%) but recovery 65 vs 39 and 17C
   51 vs 49; at eps=0.07, 17C reaches 68.5% (near dircluster+charge-override's 84%, no charge trick).
   OC makes WHOLE tracks (2.1 clus/ev) vs dircluster fragments (4.4). RECOMMEND OC eps=0.05-0.06.
4. *** SIM->REAL TRANSFER FAILS for separation (the key negative result) ***
   The oc_real.png gallery LOOKED ok (cherry-picked), but QUANTITATIVE gold-set validation
   (oc_gold.py, 73 real hand-labeled 2-track events) shows OC MERGES the two tracks badly:
     OC:  eps0.01 merge 33% (ARI0.10, shatters), eps0.03 53%, eps0.05 74%  -- NO good eps.
     dircluster (geometric): merge 1.4%, ARI 0.449, homog 0.734.
   Denoise-first + OC does NOT rescue it (merge 49-68%) -> the gap is in the track SEPARATION,
   not the noise. OC's learned latent space separates SIM tracks (1.5% merge in-domain) but real
   tracks map too close in latent space -> merge. Same failure mode as the earlier embedding GNN.

5. GRID (latent {4,6} x s_beta {0.1,1.0}, 30ep each) CONFIRMS it is a domain gap, not a
   hyperparameter artifact. In-domain sim recovery vs real-gold merge:
     lat4 sb0.1: sim 68.9% (merge1.5) | gold merge 63.0%
     lat6 sb0.1: sim 79.5% (merge2.3) | gold merge 64.4%   <- best in-domain, still 64% real merge
     lat4 sb1.0: sim 57.1%           | gold merge 63.0%
     lat6 sb1.0: sim 56.2%           | gold merge 57.5%
   Every config excels on sim and FAILS on real (57-64% merge). Higher latent helps in-domain,
   never transfers. (s_beta=1.0 hurt coords: beta term dominates -> lower recovery.)

## Verdict
OC is elegant and DOMINATES on SIM (in-domain): recovery 65-74% vs dircluster 39%, 17C 51-68% vs
49%, whole tracks. BUT it FAILS to transfer to real data for separation (33-74% merge on the gold
set) - a learned-latent-space domain gap that denoising/eps-tuning don't fix. LESSON (repeated):
for the sim->real gap in THIS problem, GEOMETRIC methods (dircluster: direction+dE/dx continuity,
1.4% merge on real gold) beat LEARNED methods (OC, embedding GNN) which overfit sim morphology.
=> Keep dircluster/AtDirDeDxCleaner for production. OC would need real training labels or much
better sim track-morphology realism (beyond noise) to be viable on real data.
