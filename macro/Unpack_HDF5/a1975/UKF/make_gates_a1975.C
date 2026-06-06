/// @file make_gates_a1975.C
/// @brief Write proton + deuteron AtParticleID gates (dedx vs brho) as Spyral JSON.
/// Polygons follow the bands seen in combined_pid_gates.png. Tune and re-run.
void make_gates_a1975()
{
   gSystem->Load("libAtTools.so");

   // Proton: tight band around the bright hyperbola (brho ~0.08-0.22)
   std::vector<std::pair<double, double>> p = {{20, 0.255}, {80, 0.19},  {200, 0.15},  {400, 0.118}, {700, 0.102},
                                               {1200, 0.094}, {2000, 0.090}, {2000, 0.072}, {1200, 0.078}, {700, 0.085},
                                               {400, 0.095}, {200, 0.115}, {80, 0.145}, {20, 0.185}};
   AtTools::AtParticleID(AtTools::AtCut2D("proton_cut", p, "dedx", "brho"), 1, 1).WriteJSON("proton_pid.json");

   // Deuteron: the band above the proton (brho ~0.20-0.36, dedx ~150-800)
   std::vector<std::pair<double, double>> d = {{120, 0.36}, {260, 0.31}, {450, 0.275}, {650, 0.24}, {820, 0.215},
                                               {820, 0.165}, {600, 0.185}, {420, 0.215}, {260, 0.25}, {120, 0.29}};
   AtTools::AtParticleID(AtTools::AtCut2D("deuteron_cut", d, "dedx", "brho"), 1, 2).WriteJSON("deuteron_pid.json");

   printf("wrote proton_pid.json, deuteron_pid.json\n");
}
