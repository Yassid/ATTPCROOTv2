/// @file apply_acc_pp.C
/// @brief Correct the 16C(p,p) ground-state angular distribution with the simulated acceptance.
///
///     dsigma/dOmega (shape)  =  N(theta_cm) / acceptance(theta_cm) / dOmega(bin)
///
/// The yield comes from pp/angdist_gs_pp.C, i.e. the hand-drawn Ex-theta_cm cut, and the
/// acceptance from macro/Simulation/ATTPC/16C_pp_a1975/acceptance_C16pp.C, which is truth-matched:
/// a generated reaction counts as reconstructed only if a fitted track exists whose angle and
/// energy match the true proton. Counting any converged fit instead credits the beam and the
/// scattered ion and inflates the acceptance badly at forward angles.
///
/// RULES BAKED IN, from the equivalent a1954 analysis:
///   * bins where the acceptance is below accMin are DROPPED, not corrected. Dividing a yield by
///     a near-zero acceptance manufactures signal out of noise, and the two end regions here
///     (below 10 deg and above 160) have acceptance consistent with zero.
///   * the statistical error carries the acceptance's own binomial error, not just the counts.
///
/// WHAT IS STILL NOT CORRECTED. The cut fraction -- how much of the ground-state peak the drawn
/// polygon contains -- varies strongly with angle: 0.96-0.98 forward, 0.084 at theta_cm 55-60,
/// 0.67 at 85, 0.045 by 140. That is a selection efficiency on top of the detector acceptance and
/// this macro does NOT correct for it, because it depends on the peak shape rather than on the
/// apparatus. Where it is small the corrected point is a lower limit on the yield.
///
///   root -b -q 'apply_acc_pp.C()'

void apply_acc_pp(TString yieldFile = "plots/angdist_gs_pp.root",
                  TString accFile = "../../../../Simulation/ATTPC/16C_pp_a1975/diagnostics/acceptance_el.root",
                  TString accHist = "hAcc_el", Double_t accMin = 0.15, TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fy = TFile::Open(here + "/" + yieldFile);
   TH1D *Y = fy && !fy->IsZombie() ? (TH1D *)fy->Get("yield") : nullptr;
   if (!Y) {
      printf("\033[1;31mno 'yield' in %s -- run angdist_gs_pp.C first\033[0m\n", yieldFile.Data());
      return;
   }
   Y = (TH1D *)Y->Clone("Ycl");
   Y->SetDirectory(nullptr);
   fy->Close();

   TFile *fa = TFile::Open(here + "/" + accFile);
   TH1D *A = nullptr;
   if (fa && !fa->IsZombie()) {
      A = (TH1D *)fa->Get(accHist);
      if (!A) { // the name carries whatever tag the acceptance was run with
         TIter nx(fa->GetListOfKeys());
         TKey *k;
         while ((k = (TKey *)nx()))
            if (TString(k->GetName()).BeginsWith("hAcc")) {
               A = (TH1D *)fa->Get(k->GetName());
               break;
            }
      }
   }
   if (!A) {
      printf("\033[1;31mno acceptance histogram in %s\033[0m\n", accFile.Data());
      return;
   }
   A = (TH1D *)A->Clone("Acl");
   A->SetDirectory(nullptr);
   fa->Close();

   auto dOmega = [](double lo, double hi) {
      return 2 * TMath::Pi() * (std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad()));
   };

   auto *hD = (TH1D *)Y->Clone("hD");
   hD->Reset();
   hD->SetDirectory(nullptr);
   auto *hRaw = (TH1D *)Y->Clone("hRaw");
   hRaw->Reset();
   hRaw->SetDirectory(nullptr);

   printf("\n  theta_cm |  counts |  acc  +- err |  raw/sr    corrected/sr   ratio\n");
   int nDrop = 0;
   for (int b = 1; b <= Y->GetNbinsX(); ++b) {
      double lo = Y->GetBinLowEdge(b), hi = lo + Y->GetBinWidth(b), c = Y->GetBinCenter(b);
      double n = Y->GetBinContent(b), en = Y->GetBinError(b);
      if (n <= 0)
         continue;
      double dO = dOmega(lo, hi);
      double raw = n / dO, eraw = en / dO;
      hRaw->SetBinContent(b, raw);
      hRaw->SetBinError(b, eraw);
      int ab = A->FindBin(c);
      double a = A->GetBinContent(ab), ea = A->GetBinError(ab);
      if (a < accMin) {
         printf("  %3.0f-%3.0f  | %7.0f | %.3f       |  %9.4g      DROPPED (acc < %.2f)\n", lo, hi, n, a, raw,
                accMin);
         ++nDrop;
         continue;
      }
      double v = raw / a;
      // the acceptance error is part of the point's error, not a separate systematic
      double e = v * std::sqrt(std::pow(en / n, 2) + std::pow(ea / std::max(a, 1e-9), 2));
      hD->SetBinContent(b, v);
      hD->SetBinError(b, e);
      printf("  %3.0f-%3.0f  | %7.0f | %.3f %.3f |  %9.4g    %9.4g   %5.2f\n", lo, hi, n, a, ea, raw, v, v / raw);
   }
   if (nDrop)
      printf("\n  %d bins dropped for acceptance below %.2f (not corrected -- dividing by ~0 makes signal)\n", nDrop,
             accMin);

   TCanvas *c1 = new TCanvas("cacc", "acceptance-corrected", 1400, 600);
   c1->Divide(2, 1);
   c1->cd(1);
   gPad->SetGridy();
   A->SetTitle("simulated acceptance;#theta_{cm} [deg];acceptance");
   A->SetMinimum(0);
   A->SetMaximum(1.15);
   A->GetXaxis()->SetRangeUser(0, 160);
   A->SetMarkerStyle(20);
   A->SetLineWidth(2);
   A->Draw("E1");
   auto *l = new TLine(0, accMin, 160, accMin);
   l->SetLineColor(kRed + 1);
   l->SetLineStyle(2);
   l->SetLineWidth(2);
   l->Draw();

   c1->cd(2);
   gPad->SetLogy();
   gPad->SetGridy();
   hRaw->SetTitle("^{16}C(p,p) g.s.;#theta_{cm} [deg];counts / sr  [arb]");
   hRaw->SetMarkerStyle(24);
   hRaw->SetMarkerColor(kGray + 2);
   hRaw->SetLineColor(kGray + 2);
   hRaw->SetLineWidth(2);
   hRaw->SetMarkerSize(1.1);
   hRaw->Draw("E1");
   hD->SetMarkerStyle(20);
   hD->SetMarkerColor(kAzure + 2);
   hD->SetLineColor(kAzure + 2);
   hD->SetLineWidth(2);
   hD->SetMarkerSize(1.2);
   hD->Draw("E1 same");
   auto *lg = new TLegend(0.55, 0.74, 0.89, 0.88);
   lg->AddEntry(hRaw, "raw", "lp");
   lg->AddEntry(hD, "acceptance-corrected", "lp");
   lg->Draw();

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/angdist_gs_acc" + tag + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/angdist_gs_acc" + tag + ".root", "RECREATE");
   hD->Write("dsdo_corrected");
   hRaw->Write("dsdo_raw");
   A->Write("acceptance");
   fo.Close();
   printf("\n  wrote %s\n\n", png.Data());
}
