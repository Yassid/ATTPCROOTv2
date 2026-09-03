/// @file run_gate_dt.C
/// @brief Launcher for draw_gate_dt.C that KEEPS THE WINDOW OPEN when ROOT is started
///        detached (no terminal).
///
/// Why this exists: `root -l macro.C` detached with stdin from /dev/null runs the macro,
/// hits EOF on the interactive prompt, and exits instantly -- the canvas flashes up and
/// disappears. gApplication->Run() enters the GUI event loop instead, so the window stays
/// until it is closed or the process is killed. Closing the canvas terminates the app.
void run_gate_dt(TString cache = "pid/pid_plane_dt_dv1104.root",
                 TString outJson = "pid/triton_d2_dv1104.json",
                 TString gateName = "triton_d2_dv1104", int Z = 1, int A = 3, TString overlay = "")
{
   gROOT->ProcessLine(".L pid/draw_gate_dt.C");
   gROOT->ProcessLine(Form("draw_gate_dt(\"%s\",\"%s\",\"%s\",%d,%d,\"%s\")",
                           cache.Data(), outJson.Data(), gateName.Data(), Z, A, overlay.Data()));
   if (gPad)
      ((TCanvas *)gPad->GetCanvas())->Connect("Closed()", "TApplication", gApplication, "Terminate()");
   printf("\n  window is open; it stays until you close it.\n");
   gApplication->Run();
}
