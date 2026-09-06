/// @file run_gate_dp.C
/// @brief Launcher for the a1975 16C(d,p)17C PROTON gate drawer.
///
/// IT MUST NOT CALL gApplication->Run(). It used to, and that caused the drawer to RESPAWN every
/// time the window was closed: `root -l pid/run_gate_dp.C` makes TRint process this file, and a
/// gApplication->Run() from inside it is a NESTED event loop. Quit()'s Terminate(0) then only
/// unwound that nested loop; TRint resumed its file list, re-executed this macro and popped a new
/// window -- measured at ~2 windows/second, unstoppable from inside the session. Both halves are
/// now fixed (Quit() is gSystem->Exit(0)), but the nesting is gone too, belt and braces.
///
/// gate_draw_dp() is NON-BLOCKING -- it builds a TGMainFrame and returns -- so ROOT simply drops
/// to its prompt with the GUI live and processing events. The only requirement is that stdin stays
/// open, or ROOT hits EOF and exits instantly with the canvas flashing past:
///
///   cd .../macro/Unpack_HDF5/a1975/D2_UKF
///   root -l 'pid/run_gate_dp.C'                                        # interactive terminal
///   setsid bash -c 'tail -f /dev/null | root -l pid/run_gate_dp.C' &   # detached
///
/// Quit, or the window's [X], exits the process. Nothing is rebuilt: the cache is the
/// species-blind plane the (d,t) gate was drawn on.
///
/// Cross-check on the ADOPTED (d,t) channel, i.e. that the locus machinery reproduces a band that
/// is already trusted:
///   root -l 'pid/run_gate_dp.C("pid/pid_plane_dt_dv1104.root","/tmp/scratch_triton.json",1,3)'
void run_gate_dp(TString cache = "pid/pid_plane_dt_dv1104.root",
                 TString outJson = "pid/proton_d2_dv1104.json", int Z = 1, int A = 1,
                 TString refP = "pid/proton_band_d2_v2.json", TString refD = "",
                 bool showLocus = true, double eBeam = 184.25, bool flipPolar = true,
                 double xMax = 55.0, double yMax = 1.6, double icLo = 900, double icHi = 1300)
{
   gROOT->ProcessLine(".L pid/gate_draw_dp.C");
   gROOT->ProcessLine(Form("gate_draw_dp(\"%s\",\"%s\",\"%s\",\"%s\",%g,%g,%g,%g,%s,%g,%s,%d,%d)",
                           outJson.Data(), cache.Data(), refP.Data(), refD.Data(), xMax, yMax, icLo,
                           icHi, showLocus ? "true" : "false", eBeam, flipPolar ? "true" : "false", Z, A));
   printf("\n  window is open; it stays until you close it (Quit, or the [X]).\n");
}
