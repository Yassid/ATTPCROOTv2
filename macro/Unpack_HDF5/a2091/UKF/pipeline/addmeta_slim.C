/// @file addmeta_slim.C
/// @brief Inject FairRoot metadata (BranchList/FileHeader/TimeBasedBranchList/cbmout)
///        into a CloneTree-made slim file so FairFileSource (the fitters) can read it.
///        BranchList lists ONLY the branches actually present in the slim tree.
///
///   root -b -q 'addmeta_slim.C("/home/yassid/a2091_C15_reco_slim/run_0138_slim.root","/home/yassid/a2091_C15_reco/run_0138_reco.root")'
void addmeta_slim(TString slimFile, TString refReco)
{
   TFile *fref = TFile::Open(refReco, "READ");
   if (!fref || fref->IsZombie()) { printf("ERR ref %s\n", refReco.Data()); return; }
   TObject *fh = fref->Get("FileHeader");
   TList *tb = (TList *)fref->Get("TimeBasedBranchList");
   TObject *cb = fref->Get("cbmout");

   TFile *f = TFile::Open(slimFile, "UPDATE");
   if (!f || f->IsZombie()) { printf("ERR slim %s\n", slimFile.Data()); return; }
   TTree *t = (TTree *)f->Get("cbmsim");

   // fresh BranchList with only the branches in this tree
   TList bl;
   for (auto b : *t->GetListOfBranches()) bl.Add(new TObjString(b->GetName()));

   f->cd();
   bl.Write("BranchList", TObject::kSingleKey);
   if (tb) tb->Write("TimeBasedBranchList", TObject::kSingleKey); else { TList e; e.Write("TimeBasedBranchList", TObject::kSingleKey); }
   if (fh) fh->Write("FileHeader");
   if (cb) cb->Write("cbmout", TObject::kSingleKey);
   f->Write("", TObject::kOverwrite);
   printf("patched %s : branches", slimFile.Data());
   for (auto o : bl) printf(" %s", o->GetName());
   printf("\n");
   f->Close(); fref->Close();
}
