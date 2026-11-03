// Flatten every per-run ATTPCROOT kin tree (which DOES carry phi, unlike the merged dd tree)
// into one CSV so the azimuthal test can be done against the Spyral parquets.
void tocsv(){
  TChain c("kin");
  int n=c.Add("/home/yassid/C15d_fit/run_*_kin_d.root");
  printf("files added: %d  entries: %lld\n", n, c.GetEntries());
  int run,ev,tr,ndf,fwd; double ke,th,phi,keX,thX,phiX,vx,vy,vz,c2;
  c.SetBranchAddress("run",&run); c.SetBranchAddress("event",&ev);
  c.SetBranchAddress("trackID",&tr); c.SetBranchAddress("ke",&ke);
  c.SetBranchAddress("theta",&th); c.SetBranchAddress("phi",&phi);
  c.SetBranchAddress("vx",&vx); c.SetBranchAddress("vy",&vy); c.SetBranchAddress("vz",&vz);
  c.SetBranchAddress("chi2ndf",&c2);
  FILE*o=fopen("/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_phi.csv","w");
  fprintf(o,"run,event,trackID,ke,theta,phi,vx,vy,vz,chi2ndf\n");
  for(Long64_t i=0;i<c.GetEntries();i++){ c.GetEntry(i);
    fprintf(o,"%d,%d,%d,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%.4f\n",run,ev,tr,ke,th,phi,vx,vy,vz,c2); }
  fclose(o); printf("written\n");
}
