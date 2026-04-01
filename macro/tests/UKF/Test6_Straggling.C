/// @file Test6_Straggling.C
/// @brief Test 6: Energy straggling stability.
///
/// Tests that the UKF converges with energy straggling ON and a 5% biased
/// initial momentum. This was the original crash scenario before Bug 6
/// was fixed (stopped sigma points set momentum to zero).
///
/// Run: root -b -q Test6_Straggling.C
#include "UKFTestHelpers.h"

void Test6_Straggling()
{
   InitTest();
   std::cout << "\n===== Test 6: Energy Straggling Stability =====" << std::endl;

   double pTrue = fTrueMom.R();
   double pBiased = pTrue * 1.05;
   Polar3DVector momP(pBiased, fTrueMom.Theta(), fTrueMom.Phi());
   XYZVector momBiased(momP);

   auto *ukf = CreateUKF(1e-3, true);
   auto res = RunUKF(ukf, fTruePos, momBiased);
   delete ukf;

   double err = std::abs(res.pVertex - pTrue) / pTrue * 100;

   std::cout << "  Biased seed:    " << pBiased << " MeV/c (5% high)" << std::endl;
   std::cout << "  Straggling ON:  p = " << res.pVertex << " MeV/c, err = " << err
             << "%, chi2/ndf = " << res.chi2ndf << ", nTouch = " << res.nTouch
             << ", converged = " << res.converged << std::endl;

   bool pass = res.converged && err < 2.0;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (must converge with straggling ON + biased seed)" << std::endl;
}
