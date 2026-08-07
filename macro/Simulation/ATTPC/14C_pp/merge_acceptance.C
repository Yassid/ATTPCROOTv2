/// @file merge_acceptance.C
/// @brief Sum the per-seed acceptance runs of one level into a single curve.
///
/// Acceptance is a RATIO, so the seeds must be combined by summing the numerator and the
/// denominator and dividing ONCE at the end -- averaging the per-seed acceptances would weight
/// a thin backward bin the same as a well-populated forward one and would not give binomial
/// errors on the total. This adds hGen and hRec, then divides with option "B".
///
/// Each seed is an independent sample only because C14_pp_sim.C was given a distinct RNG seed
/// (it had NO seeding before 2026-08-07, so parallel jobs would otherwise have been byte
/// identical). The macro re-checks that the per-seed generated counts are not all equal, which
/// is what an unseeded run would look like.
///
///   root -b -q 'merge_acceptance.C("gs","1001,1002,1003,1004,1005")'
///   root -b -q 'merge_acceptance.C("ex1","1001,1002,1003,1004,1005",6.094)'

void merge_acceptance(TString level = "gs", TString seedsCSV = "1001,1002,1003,1004,1005", Double_t resEx = 0.0,
                      TString dir = "/mnt/f/a1954_C14_acc/")
{
   gStyle->SetOptStat(0);
   TObjArray *seeds = seedsCSV.Tokenize(",");
   TH1D *gen = nullptr, *rec = nullptr;
   std::vector<double> perSeedGen;
   int nUsed = 0;

   for (int i = 0; i < seeds->GetEntries(); ++i) {
      TString sd = ((TObjString *)seeds->At(i))->GetString().Strip(TString::kBoth);
      TString tag = level + "_s" + sd;
      TString fn = dir + "acceptance_" + tag + ".root";
      TFile *f = TFile::Open(fn);
      if (!f || f->IsZombie()) {
         printf("  skip %s (no %s)\n", tag.Data(), fn.Data());
         continue;
      }
      auto *g = (TH1D *)f->Get("hGen_" + tag);
      auto *r = (TH1D *)f->Get("hRec_" + tag);
      if (!g || !r) {
         printf("  skip %s (missing hGen/hRec)\n", tag.Data());
         f->Close();
         continue;
      }
      printf("  + %-14s gen %6.0f  rec %6.0f  acc %.3f\n", tag.Data(), g->Integral(), r->Integral(),
             g->Integral() > 0 ? r->Integral() / g->Integral() : 0);
      perSeedGen.push_back(g->Integral());
      if (!gen) {
         gen = (TH1D *)g->Clone("hGen_" + level + "_sum");
         rec = (TH1D *)r->Clone("hRec_" + level + "_sum");
         gen->SetDirectory(nullptr);
         rec->SetDirectory(nullptr);
      } else {
         gen->Add(g);
         rec->Add(r);
      }
      ++nUsed;
      f->Close();
   }
   if (!gen || nUsed == 0) {
      printf("\033[1;31mnothing to merge\033[0m\n");
      return;
   }

   // An unseeded parallel run would give byte-identical samples: same generated count every time.
   bool allEqual = true;
   for (size_t k = 1; k < perSeedGen.size(); ++k)
      if (perSeedGen[k] != perSeedGen[0])
         allEqual = false;
   if (allEqual && perSeedGen.size() > 1)
      printf("\033[1;31mWARNING: every seed generated exactly %.0f reactions -- suspicious, check "
             "the seeds really differed\033[0m\n",
             perSeedGen[0]);

   auto *acc = (TH1D *)rec->Clone("hAcc_" + level + "_sum");
   acc->Divide(rec, gen, 1, 1, "B");
   acc->SetTitle(TString::Format("14C(p,p') acceptance, Ex = %.2f MeV, %d seeds;#theta_{cm} [deg];acceptance",
                                 resEx, nUsed));
   acc->SetMinimum(0);
   acc->SetMaximum(1.05);

   printf("\n===== MERGED %s : %d seeds, %.0f reactions, %.0f reconstructed, overall %.4f =====\n", level.Data(),
          nUsed, gen->Integral(), rec->Integral(), gen->Integral() > 0 ? rec->Integral() / gen->Integral() : 0);
   printf("  theta_cm      gen     reco   acceptance\n");
   for (int b = 1; b <= gen->GetNbinsX(); ++b) {
      if (gen->GetBinContent(b) < 1)
         continue;
      printf("  %3.0f-%3.0f  %7.0f  %7.0f   %.4f +- %.4f\n", gen->GetBinLowEdge(b),
             gen->GetBinLowEdge(b) + gen->GetBinWidth(b), gen->GetBinContent(b), rec->GetBinContent(b),
             acc->GetBinContent(b), acc->GetBinError(b));
   }

   TFile fo(dir + "acceptance_merged_" + level + ".root", "RECREATE");
   gen->Write();
   rec->Write();
   acc->Write();
   fo.Close();

   TCanvas *c = new TCanvas("c", "acc", 1200, 480);
   c->Divide(2, 1);
   c->cd(1);
   acc->SetMarkerStyle(20);
   acc->SetLineColor(kAzure + 2);
   acc->SetMarkerColor(kAzure + 2);
   acc->Draw("E1");
   c->cd(2);
   gen->SetLineColor(kGray + 2);
   gen->Draw("hist");
   rec->SetLineColor(kOrange + 7);
   rec->Draw("hist same");
   auto *l = new TLegend(0.45, 0.75, 0.88, 0.88);
   l->AddEntry(gen, "generated (truth)", "l");
   l->AddEntry(rec, "reconstructed", "l");
   l->Draw();
   c->SaveAs(dir + "acceptance_merged_" + level + ".png");
   printf("\nwrote %sacceptance_merged_%s.{root,png}\n\n", dir.Data(), level.Data());
}
