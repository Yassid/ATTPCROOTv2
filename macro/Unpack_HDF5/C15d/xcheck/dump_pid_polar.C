// Dump AtSpyralPID's own polar angle (the SIMPLE estimator, not the GENFIT fit) for the old
// production reco, so it can be compared against Spyral's Estimation-phase polar -- also a simple
// estimator. If the two simple estimators agree, the 0.45 deg offset lives in the FITTERS; if they
// already disagree, it is upstream in the hits or the geometry.
void pidcsv(){
  TChain c("pid");
  int n=c.Add("/home/yassid/C15d_reco/run_*_pid.root");
  printf("files %d entries %lld\n", n, c.GetEntries());
  int run,ev,tr,valid; double polar,brho,vtxZ,vtxR;
  c.SetBranchAddress("run",&run); c.SetBranchAddress("event",&ev);
  c.SetBranchAddress("trackID",&tr); c.SetBranchAddress("valid",&valid);
  c.SetBranchAddress("polar",&polar); c.SetBranchAddress("brho",&brho);
  c.SetBranchAddress("vtxZ",&vtxZ); c.SetBranchAddress("vtxR",&vtxR);
  FILE*o=fopen("/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_pidpolar.csv","w");
  fprintf(o,"run,event,trackID,polar,brho,vtxZ,vtxR\n");
  long w=0;
  for(Long64_t i=0;i<c.GetEntries();i++){ c.GetEntry(i);
    if(!valid) continue;
    fprintf(o,"%d,%d,%d,%.6f,%.5f,%.2f,%.2f\n",run,ev,tr,polar,brho,vtxZ,vtxR); w++; }
  fclose(o); printf("written %ld valid tracks\n", w);
}
