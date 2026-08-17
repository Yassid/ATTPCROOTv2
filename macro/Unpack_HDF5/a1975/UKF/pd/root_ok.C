/// @file root_ok.C
/// @brief Is a ROOT file a COMPLETE, readable cache, or wreckage left by a killed job?
///
/// Prints one line to stdout and nothing else:
///   VALID <entries>    the file opens, was closed cleanly, and holds a TTree
///   INVALID <reason>   anything else
///
/// The motivating failure: a system reboot killed cache_pd_catima.sh mid-write, leaving five
/// per-run caches that were 276 bytes of a 300-byte header. `[ -s "$out" ]` calls those files
/// "present", so a naive restart skips them and the merged cache silently loses five runs while
/// reporting success. Testing NON-EMPTY is not the same as testing USABLE.
///
/// kRecovered is rejected on purpose. ROOT can often reconstruct a tree from the baskets of a
/// file whose directory was never written, and it does so with only a warning -- that is exactly
/// the half-written state we are hunting, so a recovered file counts as failed, not as salvaged.
///
/// Zero entries is reported as VALID 0, not INVALID: a run that legitimately yields no tracks is
/// a physics statement, not a broken file, and failing it here would make every rerun rebuild it
/// forever. The caller decides whether 0 deserves a warning.
///
///   root -b -l -q 'pd/root_ok.C("/path/to/run_0115.root")'

void root_ok(const char *path)
{
   TFile *f = TFile::Open(path, "READ");
   if (!f || f->IsZombie()) {
      printf("INVALID unopenable\n");
      return;
   }
   if (f->TestBit(TFile::kRecovered)) {
      printf("INVALID recovered\n");
      f->Close();
      return;
   }

   TIter next(f->GetListOfKeys());
   while (TKey *k = (TKey *)next()) {
      TObject *o = k->ReadObj();
      if (o && o->InheritsFrom("TTree")) {
         printf("VALID %lld\n", ((TTree *)o)->GetEntries());
         f->Close();
         return;
      }
   }

   printf("INVALID no-tree\n");
   f->Close();
}
