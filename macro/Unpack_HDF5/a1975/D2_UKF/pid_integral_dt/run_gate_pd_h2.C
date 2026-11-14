/// @file run_gate_pd_h2.C
/// @brief The buttoned gate drawer (pid/gate_draw_dp.C) on the H2 float-tau planes, for the 16C(p,d)15C deuteron
///        (2026-11-11): proton target, 15C residual, Ebeam 185 (the (p,d) working point), IC 950-1350.
///        One process per plane (the drawer's Quit exits its process):
///   cd .../D2_UKF
///   setsid bash -c "tail -f /dev/null | root -l 'pid_integral_dt/run_gate_pd_h2.C(\"peak\")'" &
///   setsid bash -c "tail -f /dev/null | root -l 'pid_integral_dt/run_gate_pd_h2.C(\"integral\")'" &
/// Caches from make_h2_pid_caches.C. Gates are saved next to the cache (gate_pd_h2_<plane>.json), never over an
/// adopted gate. No reference gate is overlaid: the adopted (p,d) gate lives on the TC/PSAMax plane, another scale.
void run_gate_pd_h2(TString plane = "peak", double xMax = 45.0, double yMax = 2.0)
{
   TString dir = "/mnt/f/a1975/reco_h2_tau/";
   gROOT->ProcessLine(".L /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/gate_draw_dp.C");
   gROOT->ProcessLine("gDpTargetAmu = 1.00727646688; gDpM4Override = 15.0105993;");
   gROOT->ProcessLine(Form("gDpTitleTag = \"   ***  H2 (p,d), %s  ***\";",
                           plane == "peak" ? "PEAK AMPLITUDE" : "TRACE INTEGRAL"));
   gROOT->ProcessLine(Form("gate_draw_dp(\"%sgate_pd_h2_%s.json\",\"%spid_plane_h2_%s.root\",\"\",\"\",%g,%g,950,1350,true,"
                           "185.0,true,1,2)",
                           dir.Data(), plane.Data(), dir.Data(), plane.Data(), xMax, yMax));
   printf("\n  H2 %s plane open; it stays until you close it (Quit, or the [X]).\n", plane.Data());
}
