/// @file run_ukf_display.C
/// @brief Launch the standalone UKF fitter display for the pi-in-TPC sandbox.
///
/// Same `AtUKFDisplay` GUI as 16C_pp — works on this digi/UKF output without
/// changes (square pads, P10 gas, B=0.5 T are read from the parameter
/// containers, geometry is loaded from data/geofile_full.root).
///
/// Prerequisites:
///   1. pi_TPC_sim.C    -> data/attpcsim.root
///   2. run_digi_attpc.C -> data/output_digi.root
///   3. run_ukf_only.C   -> data/output_ukf_only.root  (optional)
///
/// Run: root -l run_ukf_display.C
///      root -l 'run_ukf_display.C(3)'   # start at event 3
///
/// Tip: AtUKFDisplay supports live re-fitting from the GUI. Particle/charge
/// dropdown defaults to proton; switch to "pi+" / "pi-" before clicking
/// "Fit" to use the right mass/charge for these events.

void run_ukf_display(int startEvent = 1)
{
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   auto *display = AtUKFDisplay::GetInstance();
   display->LoadFiles("data/output_digi.root", "data/output_ukf_only.root", "data/attpcsim.root");
   // Per-setup defaults so the "Fit" button uses sensible values out of the box.
   display->SetParticle(3);          // 3 = Pi+
   display->SetBField(2.0);          // T  (matches pi_TPC_sim.C default)
   display->SetGasDensity(1.654e-3); // g/cm^3  (P10 at 1 bar, 273 K)
   display->GotoEvent(startEvent);
}
