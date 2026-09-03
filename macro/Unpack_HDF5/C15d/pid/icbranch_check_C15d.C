/// @file icbranch_check_C15d.C
/// @brief Print HAS_ICTB if the IC summaries carry the 'ictb' branch, i.e. were made by the
///        CORRECTED icsum_C15d.C. ic_batch.sh skips runs whose output exists, so this is how a
///        driver tells "already done" from "done by the old, wrong version".
void icbranch_check_C15d()
{
   TString fn = gSystem->GetFromPipe("ls -1 /home/yassid/C15d_ic/*_ic.root 2>/dev/null | head -1");
   if (fn.Length() == 0) { printf("NO_FILES\n"); return; }
   TFile f(fn);
   auto *t = (TTree *)f.Get("ic");
   printf("%s\n", (t && t->GetBranch("ictb")) ? "HAS_ICTB" : "OLD_FORMAT");
}
