/// @file dp_export_explorer.C
/// @brief Dump every (d,p) exres cache to raw float32 for the browser explorer.
///
///   root -b -q 'dp_export_explorer.C("/mnt/f/a1954_C14dp","/tmp/dpexp")'
///
/// Seven fields per event, interleaved, little-endian float32:
///   thReco keReco exReco thTrue cmTrue chi2ndf zReco
/// One .bin per sample plus an index.txt naming them, which the packer turns into one HTML page.
/// Truth is carried alongside reco on purpose: the point of the explorer is to zoom into a
/// feature of the kinematics and ask what it is, and that question needs the truth columns.
void dp_export_explorer(TString root = "/mnt/f/a1954_C14dp", TString outDir = "/tmp/dpexp")
{
   gSystem->mkdir(outDir, kTRUE);
   const char *cfgs[6] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
   const char *levs[3] = {"gs", "ex0740", "ex3103"};
   FILE *idx = fopen((outDir + "/index.txt").Data(), "w");
   long total = 0;
   for (auto cfg : cfgs) {
      for (auto lev : levs) {
         // the seed is baked into the tag, so glob for it rather than hard-coding
         TString pat = TString(root) + "/" + cfg + "/exres_" + lev + "_*_" + cfg + ".root";
         TString found = gSystem->GetFromPipe("ls " + pat + " 2>/dev/null | head -1");
         found = found.Strip(TString::kBoth);
         if (found.IsNull()) { printf("  MISSING %s %s\n", cfg, lev); continue; }
         TFile *f = TFile::Open(found);
         if (!f || f->IsZombie()) { printf("  BAD %s\n", found.Data()); continue; }
         TTree *t = (TTree *)f->Get("res");
         if (!t) { printf("  no res tree in %s\n", found.Data()); f->Close(); continue; }
         double exReco, exTrue, thTrue, thReco, keTrue, keReco, cmTrue, chi2ndf, zTrue, zReco;
         t->SetBranchAddress("exReco", &exReco);   t->SetBranchAddress("thTrue", &thTrue);
         t->SetBranchAddress("thReco", &thReco);   t->SetBranchAddress("keReco", &keReco);
         t->SetBranchAddress("cmTrue", &cmTrue);   t->SetBranchAddress("chi2ndf", &chi2ndf);
         t->SetBranchAddress("zReco", &zReco);
         TString base = TString(cfg) + "_" + lev;
         FILE *fo = fopen((outDir + "/" + base + ".bin").Data(), "wb");
         long n = 0;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (!std::isfinite(exReco) || !std::isfinite(thReco) || !std::isfinite(keReco)) continue;
            float v[7] = {(float)thReco, (float)keReco, (float)exReco,
                          (float)thTrue, (float)cmTrue, (float)chi2ndf, (float)zReco};
            fwrite(v, sizeof(float), 7, fo);
            ++n;
         }
         fclose(fo);
         fprintf(idx, "%s %s %s %ld\n", cfg, lev, base.Data(), n);
         printf("  %-22s %6ld events\n", base.Data(), n);
         total += n;
         f->Close();
      }
   }
   fclose(idx);
   printf("\n  total %ld events -> %s\n", total, outDir.Data());
}
