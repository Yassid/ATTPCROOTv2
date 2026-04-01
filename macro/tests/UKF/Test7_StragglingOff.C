/// @file Test7_StragglingOff.C
/// @brief Test 7: UKF works with energy straggling disabled.
///
/// Previously, disabling straggling set the augmented process noise to
/// exactly zero, making the augmented covariance singular and crashing
/// the Cholesky decomposition.
///
/// Run: root -b -q Test7_StragglingOff.C
#include "UKFTestHelpers.h"

void Test7_StragglingOff()
{
   InitTest();
   std::cout << "\n===== Test 7: Straggling OFF Stability =====" << std::endl;

   double pTrue = fTrueMom.R();

   // Test with true momentum
   auto *ukf1 = CreateUKF(1e-3, false);
   auto res1 = RunUKF(ukf1, fTruePos, fTrueMom);
   delete ukf1;

   double err1 = std::abs(res1.pVertex - pTrue) / pTrue * 100;
   std::cout << "  True seed:      p = " << res1.pVertex << " MeV/c, err = " << err1
             << "%, nTouch = " << res1.nTouch << ", converged = " << res1.converged << std::endl;

   // Test with 5% biased momentum
   double pBiased = pTrue * 1.05;
   Polar3DVector momP(pBiased, fTrueMom.Theta(), fTrueMom.Phi());
   XYZVector momBiased(momP);

   auto *ukf2 = CreateUKF(1e-3, false);
   auto res2 = RunUKF(ukf2, fTruePos, momBiased);
   delete ukf2;

   double err2 = std::abs(res2.pVertex - pTrue) / pTrue * 100;
   std::cout << "  Biased seed:    p = " << res2.pVertex << " MeV/c, err = " << err2
             << "%, nTouch = " << res2.nTouch << ", converged = " << res2.converged << std::endl;

   bool pass = res1.converged && res2.converged;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (both must converge with straggling OFF)" << std::endl;
}
