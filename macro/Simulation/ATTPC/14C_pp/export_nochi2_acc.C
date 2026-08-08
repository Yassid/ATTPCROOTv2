/// @file export_nochi2_acc.C
/// @brief Write the no-chi2-cut acceptance in the layout apply_acceptance_C14.C expects.
///
/// acceptance_split_C14.C already stores the counterfactual numerator hNoChi2_<level> (truth
/// match kept, chi2/ndf cut removed) alongside the denominator hGen_<level>. This sums the seeds
/// and writes <outDir>/acceptance_merged_<level>.root holding hGen_<level>_sum, hRec_<level>_sum
/// and hAcc_<level>_sum -- the same three names merge_acceptance.C produces -- so the correction
/// macro can consume it unchanged.
///
///   root -b -q 'export_nochi2_acc.C("diagnostics/split","gs","gf","/mnt/f/a1954_C14_acc_gf_nochi2")'

void export_nochi2_acc(TString dir, TString level, TString fitter, TString outDir, Int_t s0 = 1001, Int_t s1 = 1005)
{
   TH1D *gen = nullptr, *rec = nullptr;
   int nSeed = 0;
   for (int s = s0; s <= s1; ++s) {
      TFile *fi = TFile::Open(TString::Format("%s/acc_split_%s_%s_s%d.root", dir.Data(), level.Data(), fitter.Data(), s));
      if (!fi || fi->IsZombie()) {
         printf("\033[1;31mmissing %s %s seed %d\033[0m\n", level.Data(), fitter.Data(), s);
         continue;
      }
      auto add = [&](TH1D *&dst, const char *nm, const char *as) {
         auto *h = (TH1D *)fi->Get(nm);
         if (!h)
            return;
         if (!dst) {
            dst = (TH1D *)h->Clone(as);
            dst->SetDirectory(nullptr);
         } else
            dst->Add(h);
      };
      add(gen, "hGen_" + level, "gen");
      add(rec, "hNoChi2_" + level, "rec");
      ++nSeed;
      fi->Close();
   }
   if (!gen || !rec) {
      printf("\033[1;31mnothing to export\033[0m\n");
      return;
   }
   gSystem->mkdir(outDir, kTRUE);
   auto *acc = (TH1D *)rec->Clone("hAcc_" + level + "_sum");
   acc->Divide(rec, gen, 1, 1, "B");
   acc->SetTitle(TString::Format("14C(p,p') acceptance, %s, NO chi2 cut;#theta_{cm} [deg];acceptance", level.Data()));
   gen->SetName("hGen_" + level + "_sum");
   rec->SetName("hRec_" + level + "_sum");

   TFile fo(outDir + "/acceptance_merged_" + level + ".root", "RECREATE");
   gen->Write();
   rec->Write();
   acc->Write();
   fo.Close();
   printf("%s / %s: %d seeds, %.0f generated, %.0f reconstructed, overall %.3f -> %s\n", fitter.Data(), level.Data(),
          nSeed, gen->Integral(), rec->Integral(), gen->Integral() > 0 ? rec->Integral() / gen->Integral() : 0,
          outDir.Data());
}
