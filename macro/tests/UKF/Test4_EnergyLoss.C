/// @file Test4_EnergyLoss.C
/// @brief Test 4: Energy loss consistency.
///
/// Checks that the smoothed momentum decreases monotonically (proton
/// losing energy in H2 gas) and that total energy loss is physically
/// reasonable.
///
/// Run: root -b -q Test4_EnergyLoss.C
#include "UKFTestHelpers.h"

void Test4_EnergyLoss()
{
   InitTest();
   std::cout << "\n===== Test 4: Energy Loss Consistency =====" << std::endl;

   auto *ukf = CreateUKF(1e-3, true);
   auto result = RunUKF(ukf, fTruePos, fTrueMom);

   if (!result.converged) {
      delete ukf;
      std::cout << "  FAIL: UKF did not converge" << std::endl;
      return;
   }

   auto &smoothed = ukf->GetSmoothedStates();

   // Check monotonic momentum decrease
   int nInversions = 0;
   for (size_t i = 1; i < smoothed.size(); ++i) {
      if (smoothed[i][3] > smoothed[i - 1][3] * 1.01)
         nInversions++;
   }

   double pFirst = smoothed[0][3];
   double pLast = smoothed.back()[3];
   double totalEloss = AtTools::Kinematics::KE(pFirst, mass_p) - AtTools::Kinematics::KE(pLast, mass_p);

   std::cout << "  Vertex momentum:     " << pFirst << " MeV/c" << std::endl;
   std::cout << "  Final momentum:      " << pLast << " MeV/c" << std::endl;
   std::cout << "  Total energy loss:   " << totalEloss << " MeV" << std::endl;
   std::cout << "  Momentum inversions: " << nInversions << " (steps where p increased >1%)" << std::endl;
   std::cout << "  Track length:        " << smoothed.size() << " points" << std::endl;

   bool pass = totalEloss > 0 && totalEloss < 5.0 && nInversions < 3;
   std::cout << "\n  STATUS: " << (pass ? "PASS" : "FAIL");
   std::cout << " (requires dE>0, dE<5 MeV, <3 inversions)" << std::endl;

   delete ukf;
}
