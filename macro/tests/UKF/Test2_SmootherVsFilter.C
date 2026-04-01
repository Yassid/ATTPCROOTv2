/// @file Test2_SmootherVsFilter.C
/// @brief Test 2: Smoother must improve vertex momentum over the filter.
///
/// Seeds the UKF with a 5% biased initial momentum, then checks that
/// the RTS smoother recovers the true momentum better than the forward
/// filter alone.
///
/// Run: root -b -q Test2_SmootherVsFilter.C
#include "UKFTestHelpers.h"

void Test2_SmootherVsFilter()
{
   InitTest();
   std::cout << "\n===== Test 2: Smoother vs Filter Quality =====" << std::endl;

   double pTrue = fTrueMom.R();
   double pSmeared = pTrue * 1.05; // 5% high
   Polar3DVector momPolar(pSmeared, fTrueMom.Theta(), fTrueMom.Phi());
   XYZVector momSmeared(momPolar);

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, momSmeared);
   delete ukf;

   if (!result.converged) {
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return;
   }

   double errSmooth = std::abs(result.pVertex - pTrue);
   double errFilt = std::abs(result.pVertexFilt - pTrue);

   std::cout << "  Initial momentum:    " << pSmeared << " MeV/c (5% high)" << std::endl;
   std::cout << "  True momentum:       " << pTrue << " MeV/c" << std::endl;
   std::cout << "  |p_smooth - p_true| = " << errSmooth << " MeV/c" << std::endl;
   std::cout << "  |p_filt   - p_true| = " << errFilt << " MeV/c" << std::endl;
   std::cout << "  RMS resid smooth:    " << result.rmsResidSmooth << " mm" << std::endl;
   std::cout << "  RMS resid filter:    " << result.rmsResidFilt << " mm" << std::endl;

   bool pass = errSmooth <= errFilt;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (smoother vertex error must be <= filter vertex error)" << std::endl;
}
