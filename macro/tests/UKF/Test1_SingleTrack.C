/// @file Test1_SingleTrack.C
/// @brief Test 1: Single track momentum reconstruction with known truth.
///
/// Verifies that the UKF reconstructs the vertex momentum to <1% error
/// with zero covariance regularizations.
///
/// Run: root -b -q Test1_SingleTrack.C
#include "UKFTestHelpers.h"

void Test1_SingleTrack()
{
   InitTest();
   std::cout << "\n===== Test 1: Single Track Momentum Reconstruction =====" << std::endl;

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, fTrueMom);
   delete ukf;

   if (!result.converged) {
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return;
   }

   double pTrue = fTrueMom.R();
   double errSmooth = (result.pVertex - pTrue) / pTrue * 100;
   double errFilt = (result.pVertexFilt - pTrue) / pTrue * 100;

   std::cout << "  True momentum:       " << pTrue << " MeV/c" << std::endl;
   std::cout << "  Smoothed momentum:   " << result.pVertex << " MeV/c (err: " << errSmooth << "%)" << std::endl;
   std::cout << "  Filtered momentum:   " << result.pVertexFilt << " MeV/c (err: " << errFilt << "%)" << std::endl;
   std::cout << "  Chi2/ndf:            " << result.chi2ndf << std::endl;
   std::cout << "  RMS resid (smooth):  " << result.rmsResidSmooth << " mm" << std::endl;
   std::cout << "  RMS resid (filter):  " << result.rmsResidFilt << " mm" << std::endl;
   std::cout << "  Regularizations:     " << result.nTouch << std::endl;

   bool pass = std::abs(errSmooth) < 1.0 && result.nTouch == 0;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires <1% error, 0 regularizations)" << std::endl;
}
