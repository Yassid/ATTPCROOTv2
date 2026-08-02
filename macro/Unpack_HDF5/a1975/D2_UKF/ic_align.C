// Verify the IC file is entry-aligned with the reco/fit tree: the GET and FRIB timestamps
// must differ by the SAME constant at every entry of the run.
void ic_align(const char *runName = "run_0016")
{
   gSystem->Load("libAtReconstruction.so");
   TString run = runName;
   TString icf = "/mnt/f/a1975/ic_d2/" + run + "_IC.root";
   if (gSystem->AccessPathName(icf)) { printf("%s: no IC file yet\n", run.Data()); return; }
   TFile *fi = TFile::Open(icf);
   TTree *ti = (TTree *)fi->Get("ic");
   ULong64_t ts = 0; int ev = 0; float icm = 0;
   ti->SetBranchAddress("ts", &ts); ti->SetBranchAddress("evt", &ev); ti->SetBranchAddress("icmax", &icm);

   TFile *fr = TFile::Open("/mnt/f/a1975/reco_d2/" + run + "_reco.root");
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *ea = nullptr;
   tr->SetBranchAddress("AtEventCorrected", &ea);
   printf("%s: ic %lld entries, reco %lld entries\n", run.Data(), ti->GetEntries(), tr->GetEntries());

   Long64_t probes[6] = {0, 100, 5000, 20000, 30000, 40000};
   double off0 = 0; int nchk = 0, nbad = 0;
   for (int p = 0; p < 6; ++p) {
      Long64_t i = probes[p];
      if (i >= tr->GetEntries() || i >= ti->GetEntries()) continue;
      ti->GetEntry(i); tr->GetEntry(i);
      if (ea->GetEntries() == 0) continue;
      AtEvent *e = (AtEvent *)ea->At(0);
      // the GET clock wraps at 2^32 partway through a run, so compare modulo 2^32
      const double W = 4294967296.0;
      double d = (double)e->GetTimestamp(1) - (double)ts;
      d = fmod(fmod(d, W) + W, W);
      ++nchk;
      if (nchk == 1) off0 = d; else if (fabs(d - off0) > 2) ++nbad;
      printf("   entry %6lld  icEvt %6d  ic %7.0f  dTS %.0f\n", i, ev, icm, d);
   }
   printf("   -> %s\n", nbad ? "MISALIGNED" : "aligned (constant timestamp offset)");
}
