/// @file slim_cache_C14.C
/// @brief Clone ONLY AtPatternEvent (+EventHeader) from a big reco.root on F into a
///        tiny LOCAL file, so downstream PID/fit passes don't re-read 40 GB over drvfs.
///        Reads only the AtPatternEvent baskets from the source (skips AtEventH).
///
///   root -b -q 'pipeline/slim_cache_C14.C("run_0055","/mnt/f/a1954_C14_reco_hdb/","/home/yassid/a1954_C14_reco_hdb_slim/")'
void slim_cache_C14(TString run = "run_0055", TString inDir = "/mnt/f/a1954_C14_reco_hdb/",
                     TString outDir = "/home/yassid/a1954_C14_reco_hdb_slim/")
{
   TString rf = inDir + run + "_reco.root";
   if (gSystem->AccessPathName(rf)) {
      printf("MISSING %s\n", rf.Data());
      return;
   }
   gSystem->mkdir(outDir.Data(), kTRUE);
   TString of = outDir + run + "_slim.root";

   TFile *fin = TFile::Open(rf);
   TTree *t = (TTree *)fin->Get("cbmsim");
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchStatus("EventHeader*", 1);

   TFile *fout = new TFile(of, "RECREATE", "", 1); // low compression: files are tiny, favor speed
   TTree *nt = t->CloneTree(-1, "fast");
   nt->Write();
   printf("%s : %lld events -> %s (%.2f MB)\n", run.Data(), nt->GetEntries(), of.Data(), fout->GetSize() / 1e6);
   fout->Close();
   fin->Close();
}
