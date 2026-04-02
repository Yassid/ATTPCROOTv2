/// @file run_ukf_display.C
/// @brief Launch the standalone UKF fitter display.
///
/// Prerequisites:
///   1. Run run_digi_attpc.C → output_digi.root
///   2. Optionally run run_ukf_only.C → output_ukf_only.root
///
/// Run: root -l run_ukf_display.C
///      root -l 'run_ukf_display.C(3)' // start at event 3
///
/// For Eve/OpenGL mode (requires working OpenGL):
///   root -l -e 'TEveManager::Create()' run_ukf_display.C

void run_ukf_display(int startEvent = 1)
{
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   auto *display = AtUKFDisplay::GetInstance();
   display->LoadFiles("data/output_digi.root", "data/output_ukf_only.root", "data/attpcsim.root");
   display->GotoEvent(startEvent);
}
