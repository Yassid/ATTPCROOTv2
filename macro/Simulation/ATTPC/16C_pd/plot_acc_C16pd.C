/// @file plot_acc_C16pd.C
/// @brief The four 16C(p,d)15C acceptances on one figure, and their ratios to the ground state.
///
/// The reason to draw them together is that the INTEGRATED numbers are nearly identical -- 0.735,
/// 0.741, 0.734, 0.754 -- so a single figure per state, or a single number, hides the thing that
/// matters. The states differ at the two ends of the angular range and agree in the middle, and
/// the bottom panel is what shows it: the ratio to the ground state is flat to about 5 percent
/// between 20 and 140 degrees and departs badly outside that.
///
/// Acceptance here is truth-matched: a generated reaction counts in the numerator only if a
/// fitted track exists whose angle and energy match the true deuteron. Counting any converged fit
/// instead would credit the beam and the scattered ion, which in the a1954 analysis produced a
/// fake acceptance of about 0.9 on protons that had never made a track.
///
///   root -b -q 'plot_acc_C16pd.C()'

void plot_acc_C16pd(TString dir = "diagnostics/", TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const int NS = 4;
   const char *T[NS] = {"gs", "ex1", "ex2", "ex3"};
   const char *NM[NS] = {"g.s.", "0.740 MeV", "3.10 MeV", "4.66 MeV"};
   const int COL[NS] = {kBlack, kAzure + 2, kRed + 1, kGreen + 3};

   TH1D *acc[NS] = {nullptr, nullptr, nullptr, nullptr};
   for (int s = 0; s < NS; ++s) {
      TFile *f = TFile::Open(here + "/" + dir + "acceptance_" + T[s] + ".root");
      if (!f || f->IsZombie()) {
         printf("\033[1;31mmissing acceptance_%s.root\033[0m\n", T[s]);
         continue;
      }
      // the histogram name carries the tag the acceptance macro was run with
      TIter nx(f->GetListOfKeys());
      TKey *k;
      while ((k = (TKey *)nx())) {
         TString n = k->GetName();
         if (n.BeginsWith("hAcc")) {
            acc[s] = (TH1D *)f->Get(n);
            acc[s]->SetDirectory(nullptr);
            break;
         }
      }
      f->Close();
   }
   if (!acc[0]) {
      printf("\033[1;31mno ground-state acceptance -- nothing to normalise to\033[0m\n");
      return;
   }

   TCanvas *c = new TCanvas("cacc", "acceptance", 1250, 900);
   c->Divide(1, 2);

   c->cd(1);
   gPad->SetGridy();
   auto *fr = new TH1D("fr", "16C(p,d)^{15}C acceptance, per final state;#theta_{cm} [deg];acceptance", 1, 0, 180);
   fr->SetMinimum(0);
   fr->SetMaximum(1.15);
   fr->Draw();
   auto *lg = new TLegend(0.13, 0.14, 0.42, 0.36);
   for (int s = 0; s < NS; ++s) {
      if (!acc[s])
         continue;
      acc[s]->SetLineColor(COL[s]);
      acc[s]->SetMarkerColor(COL[s]);
      acc[s]->SetMarkerStyle(20 + s);
      acc[s]->SetMarkerSize(1.0);
      acc[s]->SetLineWidth(2);
      acc[s]->Draw("E1 same");
      lg->AddEntry(acc[s], NM[s], "lp");
   }
   lg->Draw();
   // the range the states agree over, and therefore the range a common curve could be used in
   for (double x : {20.0, 140.0}) {
      auto *l = new TLine(x, 0, x, 1.15);
      l->SetLineColor(kGray + 2);
      l->SetLineStyle(2);
      l->SetLineWidth(2);
      l->Draw();
   }

   c->cd(2);
   gPad->SetGridy();
   auto *fr2 = new TH1D("fr2", "ratio to the ground state;#theta_{cm} [deg];acceptance / acceptance(g.s.)", 1, 0, 180);
   fr2->SetMinimum(0);
   fr2->SetMaximum(2.0);
   fr2->Draw();
   printf("\n  ratio to the g.s., in 20 deg bins\n  theta_cm |");
   for (int s = 1; s < NS; ++s)
      printf("  %-9s", NM[s]);
   printf("\n");
   TH1D *rat[NS] = {nullptr, nullptr, nullptr, nullptr};
   for (int s = 1; s < NS; ++s) {
      if (!acc[s])
         continue;
      rat[s] = (TH1D *)acc[s]->Clone(TString::Format("r%d", s));
      rat[s]->SetDirectory(nullptr);
      rat[s]->Divide(acc[s], acc[0]);
      rat[s]->SetLineColor(COL[s]);
      rat[s]->SetMarkerColor(COL[s]);
      rat[s]->SetMarkerStyle(20 + s);
      rat[s]->SetLineWidth(2);
      rat[s]->Draw("E1 same");
   }
   auto *one = new TLine(0, 1, 180, 1);
   one->SetLineColor(kBlack);
   one->SetLineWidth(2);
   one->Draw();
   for (double x : {20.0, 140.0}) {
      auto *l = new TLine(x, 0, x, 2.0);
      l->SetLineColor(kGray + 2);
      l->SetLineStyle(2);
      l->SetLineWidth(2);
      l->Draw();
   }
   for (double lo = 0; lo < 180; lo += 20) {
      printf("   %3.0f-%3.0f |", lo, lo + 20);
      for (int s = 1; s < NS; ++s) {
         if (!rat[s]) {
            printf("     -    ");
            continue;
         }
         // average the ratio over the bins in this 20 deg window, skipping empty ones
         double sum = 0;
         int n = 0;
         for (int b = 1; b <= rat[s]->GetNbinsX(); ++b) {
            double x = rat[s]->GetBinCenter(b);
            if (x < lo || x >= lo + 20 || rat[s]->GetBinContent(b) <= 0)
               continue;
            sum += rat[s]->GetBinContent(b);
            ++n;
         }
         if (n)
            printf("   %.3f  ", sum / n);
         else
            printf("     -    ");
      }
      printf("\n");
   }

   TString png = here + "/plots/acceptance_C16pd" + tag + ".png";
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s\n\n", png.Data());
}
