/// @file dlc_prf_sigma.C
/// @brief Map the measured PUMA DLC sheet resistance to the resistive PRF sigma
///        used by AtPulse::SetChargeDispersion / the prfSigma sweep macros.
///
/// The resistive DLC anode spreads the avalanche charge (Dixit/telegraph) as
///     sigma = sqrt( 2 t / (R C) )
///   R = DLC sheet resistance      [Ohm/square]  -- PUMA measured 1.2-1.5e6
///   C = areal capacitance DLC<->pads [F/mm^2]
///   t = charge spreading time     [us]          -- ~ shaping/peaking time (0.5)
///
/// Use either:
///   (a) AtPulse::SetChargeDispersionFromDLC(R, C, t_us)  -- set it physically, or
///   (b) pass the sigma below as prfSigma to run_digi_ukf_genfit_test8 /
///       SetChargeDispersion(sigma_mm).
///
/// Run: root -b -q dlc_prf_sigma.C
double dlc_sigma_mm(double R_ohm_sq, double C_F_mm2, double t_us)
{
   double RC = R_ohm_sq * C_F_mm2; // s/mm^2
   return (RC > 0 && t_us > 0) ? std::sqrt(2.0 * (t_us * 1e-6) / RC) : 0.0;
}

void dlc_prf_sigma(double t_us = 0.5)
{
   printf("\nPUMA DLC -> PRF sigma   (t = %.2f us)\n", t_us);
   printf("Areal capacitance C set by the DLC insulator gap (CONFIRM from detector build).\n");
   printf("Columns: C from an ERAM-like insulator of thickness d (eps_r=3.4).\n\n");
   const double eps0 = 8.854e-12; // F/m
   auto C_of_d = [&](double d_um) { return eps0 * 3.4 / (d_um * 1e-6) / 1e6; }; // F/mm^2
   printf("%-14s %-14s %-14s %-14s\n", "R [MOhm/sq]", "sigma(d=25um)", "sigma(d=50um)", "sigma(d=75um)");
   for (double Rm : {1.20, 1.35, 1.50}) {
      printf("%-14.2f %-14.3f %-14.3f %-14.3f\n", Rm, dlc_sigma_mm(Rm * 1e6, C_of_d(25), t_us),
             dlc_sigma_mm(Rm * 1e6, C_of_d(50), t_us), dlc_sigma_mm(Rm * 1e6, C_of_d(75), t_us));
   }
   printf("\nT2K-empirical cross-check (RC=100 ns/mm^2 at R=0.4 MOhm, RC ~ R):\n");
   for (double Rm : {1.20, 1.35, 1.50}) {
      double RC_ns = 100.0 * (Rm / 0.4); // ns/mm^2
      printf("  R=%.2f MOhm/sq -> sigma = %.3f mm\n", Rm, std::sqrt(2.0 * t_us * 1e3 / RC_ns));
   }
   printf("\n=> sigma_DLC ~ 1.0-1.8 mm (central ~1.1 mm) -> charge over ~1-2 pads (2 mm pitch).\n");
   printf("   Example: atPulse->SetChargeDispersionFromDLC(1.35e6, %.3e, %.2f); // sigma=%.3f mm\n",
          C_of_d(50), t_us, dlc_sigma_mm(1.35e6, C_of_d(50), t_us));
}
