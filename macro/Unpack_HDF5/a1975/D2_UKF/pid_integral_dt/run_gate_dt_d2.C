/// @file run_gate_dt_d2.C
/// @brief The buttoned gate drawer (pid/gate_draw_dp.C) on the D2 float-tau planes, for the 16C(d,t)15C TRITON (2026-11-17):
///        deuterium target (drawer default), 15C residual (from Z=1 A=3), Ebeam 184.0 (the (d,t) working point), IC 900-1400.
///        Twin of run_gate_pd_h2.C. One process per plane (the drawer's Quit exits its process):
///   cd .../D2_UKF
///   setsid bash -c "tail -f /dev/null | root -l 'pid_integral_dt/run_gate_dt_d2.C(\"integral\")'" &
///   setsid bash -c "tail -f /dev/null | root -l 'pid_integral_dt/run_gate_dt_d2.C(\"peak\")'" &
/// Caches from make_d2_pid_caches.C. Gates are saved next to the cache (gate_dt_d2_<plane>.json), never over the adopted
/// pid/triton_d2_dv1104.json. That gate is overlaid on the PEAK plane only (it was drawn on the fixed-tau peak plane; the integral
/// plane is another scale).
void run_gate_dt_d2(TString plane = "integral", double xMax = 55.0, double yMax = 2.6)
{
   TString dir = "/mnt/f/a1975/reco_d2_hdbtau/";
   gROOT->ProcessLine(".L /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/gate_draw_dp.C");
   gROOT->ProcessLine(Form("gDpTitleTag = \"   ***  D2 (d,t) TRITON, %s  ***\";", plane == "peak" ? "PEAK AMPLITUDE" : "TRACE INTEGRAL"));
   TString ref = plane == "peak" ? "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/pid/triton_d2_dv1104.json" : "";
   gROOT->ProcessLine(Form("gate_draw_dp(\"%sgate_dt_d2_%s.json\",\"%spid_plane_d2_%s.root\",\"%s\",\"\",%g,%g,900,1400,true,"
                           "184.0,true,1,3)",
                           dir.Data(), plane.Data(), dir.Data(), plane.Data(), ref.Data(), xMax, yMax));
   printf("\n  D2 (d,t) %s plane open; it stays until you close it (Quit, or the [X]).\n", plane.Data());
}
