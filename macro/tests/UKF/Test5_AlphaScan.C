/// @file Test5_AlphaScan.C
/// @brief Test 5: Alpha parameter sensitivity scan.
///
/// Runs the UKF with different alpha values and checks that reconstruction
/// accuracy and stability are acceptable across the range.
///
/// Run: root -b -q Test5_AlphaScan.C
#include "UKFTestHelpers.h"

void Test5_AlphaScan()
{
   InitTest();
   std::cout << "\n===== Test 5: Alpha Parameter Sensitivity =====" << std::endl;

   double alphas[] = {1e-3, 5e-3, 1e-2, 5e-2};
   int nAlphas = sizeof(alphas) / sizeof(alphas[0]);

   std::cout << std::setw(10) << "alpha" << std::setw(15) << "p_smooth" << std::setw(12) << "err(%)"
             << std::setw(12) << "chi2/ndf" << std::setw(10) << "nTouch" << std::setw(10) << "status" << std::endl;
   std::cout << std::string(69, '-') << std::endl;

   double pTrue = fTrueMom.R();
   bool allPass = true;

   for (int a = 0; a < nAlphas; ++a) {
      auto *ukf = CreateUKF(alphas[a], true);
      auto result = RunUKF(ukf, fTruePos, fTrueMom);
      delete ukf;

      std::string status;
      if (!result.converged) {
         status = "FAIL";
         allPass = false;
         std::cout << std::setw(10) << alphas[a] << std::setw(15) << "---" << std::setw(12) << "---"
                   << std::setw(12) << "---" << std::setw(10) << "---" << std::setw(10) << status << std::endl;
      } else {
         double err = (result.pVertex - pTrue) / pTrue * 100;
         bool ok = std::abs(err) < 2.0 && result.nTouch <= 2;
         if (!ok)
            allPass = false;
         status = ok ? "OK" : "WARN";

         std::cout << std::setw(10) << alphas[a] << std::setw(15) << result.pVertex << std::setw(12) << err
                   << std::setw(12) << result.chi2ndf << std::setw(10) << result.nTouch << std::setw(10) << status
                   << std::endl;
      }
   }

   std::cout << "\n  STATUS: " << (allPass ? "PASS" : "WARN") << std::endl;
}
