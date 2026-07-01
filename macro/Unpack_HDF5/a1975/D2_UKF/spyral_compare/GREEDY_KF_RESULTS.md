# Greedy KF spiral-recovery follower — results (a1975 16C(d,p), run_0016)

Goal: recover the multi-turn-spiral track portions that ATTPCROOT's PRA
(`AtTrackFinderTC`) drops (~8% severe + ~20% partial of substantial events), as a
**post-PRA** step that reuses the existing PRA arc as a seed.

## What was built
`scripts/follower.py` — **grow + bridge + grow** recovery (the "greedy KF" slice):
1. Seed = PRA's largest track.
2. **Region-grow** by connectivity (`dmerge=10 mm`), absorbing only *dense* hits
   (`min_dense=3` neighbours) so isolated noise is rejected.
3. **Bridge**: from each PCA endpoint, jump to the nearest *forward, dense*
   unassigned hit within `dbridge=45 mm` (crosses inter-turn gaps without grabbing noise).
4. Repeat until stable.

Inputs: `scripts/extract_follower.C` -> `follower_hits.csv` (corrected hits + PRA-seed
label per hit, via `AtHit::GetHitID()`). Truth proxy = Spyral's largest cluster size.
Runner: `scripts/follow_validate.py` (334-event validation set, runs in <90 s).

## Results (validation set, 332 events, Spyral main cluster >=30 hits)
Recovered fraction = track size / Spyral-main-cluster (1.0 = matches Spyral).

| category      |  n  | PRA frac | KF frac | over-grow >1.3 |
|---------------|-----|----------|---------|----------------|
| severe (<30%) | 149 |   0.19   | **0.28**|       3%       |
| partial 30-70%|  60 |   0.52   | **0.70**|       5%       |
| match 70-130% |  59 |   1.05   |   1.08  |       8%       |
| clean single  |  60 |   1.04   |   1.08  |       3%       |

Case studies: ev57 0.43->**0.89** (recovers the dropped START turns), ev11090
0.35->0.92 (stays on the thin track, ignores the noise blob), ev15 1.47->1.55
(harmless), **ev49 0.18->0.18 (FAILS)**.

## Verdict
- **Works** when the PRA seed is *connected* to the dropped body (ev57, partial pop):
  recovers a lot, safely (controls barely move, noise rejected).
- **Hard failure** when the seed is a small arc *separated* from the track body by a
  gap and lying "ahead" of it (ev49): the bridge can't re-acquire the body. ~half the
  severe population is this type, so severe median only reaches 0.28.
- Compared to the one-line **loose-TC retune** (`Atriplet a 0.03->0.20, Tcluster t
  4->12`), which got ev49 0.88 / ev57 0.75 / ev36 0.85 but *destroys* messy multi-track
  events (ev11090 -> 0.04): the follower recovers **less** on spirals but is **gentler**
  (never over-merges crossing tracks; small, guarded noise pickup).

## Recommendation
Neither is a clean global win. For maximum spiral recovery, **adaptive loose-TC**
(detect high-curvature events, re-run `AtTrackFinderTC` with loose a/t on just those)
beats the geometric follower and reuses the validated C++ algorithm. The greedy
follower is worth keeping as a **safe top-up** for the connected-tail cases where
loose-TC's over-merge risk is unacceptable. Open improvement for the follower:
cluster-level (not hit-level) bridging to fix the isolated-seed case (ev49).

Next decision point: re-reco a subset with adaptive loose-TC and measure the actual
fit-yield / Ex-resolution impact before investing further.

---

# Ex-impact experiment — VERDICT (adaptive loose-TC, run_0016)

Re-ran PRA-only with loose TC (a=0.20,t=12), fit, built 17C Ex spectra, compared
default vs global-loose vs adaptive (loose only on flagged multi-turn spirals).

| variant  | proton candidates | near g.s. (|Ex|<1.5) |
|----------|-------------------|----------------------|
| default  |       2129        |         53           |
| loose    |       1178 (-45%) |         24           |
| adaptive |       2042 (-87)  |         49           |

**Recovering the spirals does NOT improve the Ex spectrum — it slightly degrades it.**
- The 5419 flagged spirals (median 3.1 turns) yield a proton candidate only **2.1%** of
  the time -> ~98% are NOT (d,p) protons, they are low-energy stopped background
  (tight spiral = low momentum).
- For the few that are protons, fully clustering the spiral makes them fit WORSE
  (115 -> 28 proton candidates): a multi-turn helix is harder for genfit than the clean
  partial arc PRA already keeps.
- Global loose-TC is actively harmful (-45% yield: over-merges normal events).

**Conclusion: PRA's spiral-dropping is a NON-ISSUE for 17C(d,p). Do NOT deploy adaptive
loose-TC or the greedy follower in production — they add complexity for zero (or slightly
negative) physics gain. Default `AtTrackFinderTC` is the correct choice for this channel.**
The whole "PRA drops spiral portions" thread is closed: real, but irrelevant to the
proton physics. (Would only matter for an analysis whose signal IS the low-energy
spiraling species.)
