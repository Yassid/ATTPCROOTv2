/// @file make_proton_gate_a1975.C
/// @brief Build a proton AtParticleID gate (sqrt_dEdx vs brho) and write it as
/// spyral_utils-format JSON. Also round-trips it back and verifies a real Spyral
/// gate file loads, demonstrating the C++ PID port is Spyral-compatible.
///
/// Run: root -b -q make_proton_gate_a1975.C

void make_proton_gate_a1975()
{
   gSystem->Load("libAtTools.so");

   // Polygon enclosing the proton locus seen in run_0116_pid.png
   // (sqrt_dEdx on x, brho [T m] on y).
   std::vector<std::pair<double, double>> v = {
      {0.4, 0.25}, {2.0, 0.18}, {5.0, 0.12}, {9.0, 0.09}, {13.0, 0.06},
      {13.0, 0.18}, {9.0, 0.25}, {5.0, 0.35}, {2.0, 0.50}, {0.4, 0.60},
   };
   AtTools::AtCut2D cut("proton_cut", v, "sqrt_dEdx", "brho");
   AtTools::AtParticleID pid(cut, 1, 1); // Z=1, A=1 -> proton
   const char *out = "proton_gate_a1975.json";
   pid.WriteJSON(out);
   printf("Wrote %s\n", out);

   // Round-trip: load it back and sanity-check a point inside/outside.
   auto reloaded = AtTools::AtParticleID::LoadJSON(out);
   printf("Reloaded: name='%s' Z=%d A=%d mass=%.3f MeV, %zu vertices\n", reloaded.GetName().c_str(), reloaded.GetZ(),
          reloaded.GetA(), reloaded.GetMassMeV(), reloaded.GetCut().GetVertices().size());
   printf("  point (3.0, 0.30) inside? %d   (point 3.0, 2.0) inside? %d\n", reloaded.IsInside(3.0, 0.30),
          reloaded.IsInside(3.0, 2.0));

   // Verify a REAL Spyral gate file loads (compatibility check).
   const char *spyralGate = "/home/yassid/C16_dp/C16_dp/programs/Spyral-1.0/PID/gate_pd.json";
   if (gSystem->AccessPathName(spyralGate) == 0) {
      auto sp = AtTools::AtParticleID::LoadJSON(spyralGate);
      printf("Spyral gate '%s': name='%s' axes=%s/%s, %zu vertices  (Z=%d A=%d)\n", spyralGate, sp.GetName().c_str(),
             sp.GetXAxis().c_str(), sp.GetYAxis().c_str(), sp.GetCut().GetVertices().size(), sp.GetZ(), sp.GetA());
   }
}
