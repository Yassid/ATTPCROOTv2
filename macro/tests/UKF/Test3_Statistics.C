/// @file Test3_Statistics.C
/// @brief Test 3: Statistical consistency over many trials.
///
/// Smears the initial momentum by 5% (Gaussian) over N trials and checks:
///   - Convergence rate >95%
///   - Pull distribution ~ N(0,1)
///   - Momentum error distribution
///   - Chi2/ndf distribution
///
/// Run: root -b -q 'Test3_Statistics.C(200)'
///      root -b -q 'Test3_Statistics.C(50)'   // quick
#include "UKFTestHelpers.h"

void Test3_Statistics(int nTrials = 200)
{
   InitTest();
   std::cout << "\n===== Test 3: Statistical Consistency (" << nTrials << " trials) =====" << std::endl;

   double pTrue = fTrueMom.R();
   double smearSigma = 0.05 * pTrue;

   auto *hPull = new TH1F("hPull", "Momentum Pull;(p_{reco} - p_{true}) / #sigma_{p};Entries", 50, -5, 5);
   auto *hError = new TH1F("hError", "Momentum Error;(p_{reco} - p_{true}) / p_{true} [%];Entries", 50, -3, 3);
   auto *hChi2 = new TH1F("hChi2", "#chi^{2}/ndf;#chi^{2}/ndf;Entries", 50, 0, 5);

   auto *ukf = CreateUKF(1e-3, true);

   int nConverged = 0;
   int nFailed = 0;

   for (int t = 0; t < nTrials; ++t) {
      if (t % 50 == 0)
         std::cout << "  Trial " << t << "/" << nTrials << std::endl;

      double pSampled = gRandom->Gaus(pTrue, smearSigma);
      Polar3DVector momPolar(pSampled, fTrueMom.Theta(), fTrueMom.Phi());
      XYZVector momSmeared(momPolar);

      auto result = RunUKF(ukf, fTruePos, momSmeared);

      if (!result.converged) {
         nFailed++;
         continue;
      }
      nConverged++;

      double err = (result.pVertex - pTrue) / pTrue * 100;
      hError->Fill(err);
      hChi2->Fill(result.chi2ndf);

      auto &smoothedCov = ukf->GetSmoothedCovariances();
      if (!smoothedCov.empty()) {
         double sigmaP = std::sqrt(smoothedCov[0](3, 3));
         if (sigmaP > 0) {
            double pull = (result.pVertex - pTrue) / sigmaP;
            hPull->Fill(pull);
         }
      }
   }

   delete ukf;

   double convergenceRate = 100.0 * nConverged / nTrials;
   std::cout << "\n  Converged: " << nConverged << "/" << nTrials << " (" << convergenceRate << "%)" << std::endl;

   if (hPull->GetEntries() > 20) {
      hPull->Fit("gaus", "Q");
      auto *fit = hPull->GetFunction("gaus");
      if (fit)
         std::cout << "  Pull:  mean=" << fit->GetParameter(1) << "  sigma=" << fit->GetParameter(2)
                   << "  (expect ~0, ~1)" << std::endl;
   }
   if (hError->GetEntries() > 20) {
      hError->Fit("gaus", "Q");
      auto *fit = hError->GetFunction("gaus");
      if (fit)
         std::cout << "  Error: mean=" << fit->GetParameter(1) << "%  sigma=" << fit->GetParameter(2) << "%"
                   << std::endl;
   }
   std::cout << "  Chi2/ndf mean: " << hChi2->GetMean() << " (expect ~1)" << std::endl;

   bool pass = convergenceRate > 95.0;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires >95% convergence rate)" << std::endl;

   // Save plots
   TCanvas *c = new TCanvas("c3", "Statistical Validation", 1500, 500);
   c->Divide(3, 1);
   c->cd(1);
   hPull->Draw("hist");
   if (hPull->GetFunction("gaus"))
      hPull->GetFunction("gaus")->Draw("same");
   c->cd(2);
   hError->Draw("hist");
   if (hError->GetFunction("gaus"))
      hError->GetFunction("gaus")->Draw("same");
   c->cd(3);
   hChi2->Draw("hist");
   c->SaveAs("Test3_Statistics.png");
   std::cout << "  Saved: Test3_Statistics.png" << std::endl;
}
