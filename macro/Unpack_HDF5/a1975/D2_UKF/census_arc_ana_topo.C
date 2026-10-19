/// @file census_arc_ana_topo.C -- census_arc_ana.C with the VERTEX-ARC topology required (dirB, d0<10, rzmin<30, zmax>=850): reads census_arc_dp.C output (data_*.root + sim_gs_s4001.root) for the bound window.
/// Per theta_cm bin, vz range [z0,z1): pattern-level candidates from the HIT helix estimate (C), how many AtSpyralPID
/// accepts (V), how many it also puts in the bound window itself (VB), and the failCode breakdown of the rest.
/// Then the pattern-level and Spyral-level data/sim ratios, each normalised to the 35-40 bin.
void census_arc_ana_topo(double z0 = 200, double z1 = 400, double exHi = 0.729)
{
   TChain d("c"), s("c");
   d.Add("/mnt/f/a1975/census_arc/data_*.root");
   s.Add("/mnt/f/a1975/census_arc/sim_gs_s4001.root");
   const int NT = 7;
   double TB[NT + 1] = {10, 15, 20, 25, 30, 35, 40, 45};
   TString base = Form("dirB&&d0<10&&rzmin<30&&zmax>=850&&thl>89&&thl<161&&ke>1.2&&ex>-1.5&&ex<%g&&vz>=%g&&vz<%g", exHi, z0, z1);
   TString spyB = Form("sv&&sth>89&&sth<161&&ske>1.2&&sex>-1.5&&sex<0.729&&svz>=%g&&svz<%g", z0, z1);
   double C[2][NT], V[2][NT], VB[2][NT], SB[2][NT];
   printf("vz %g-%g, my Ex < %g;  C = hit-helix bound candidates, V = Spyral valid, VB = Spyral also bound\n", z0, z1, exHi);
   printf("SB = all Spyral-bound tracks in the bin (Spyral's own kinematics, what effchain counts as N0)\n");
   for (int k = 0; k < 2; ++k) {
      TChain &t = k ? s : d;
      printf("\n=== %s  (%lld tracks)\n thcm   |     C     V/C   VB/C |  fail: 1     2     3     8     9  other |   SB\n", k ? "SIM g.s." : "DATA", t.GetEntries());
      for (int b = 0; b < NT; ++b) {
         TString c = base + Form("&&tcm>=%g&&tcm<%g", TB[b], TB[b + 1]);
         TString cs = Form("&&stcm>=%g&&stcm<%g", TB[b], TB[b + 1]);
         C[k][b] = t.GetEntries(c);
         V[k][b] = t.GetEntries(c + "&&sv");
         VB[k][b] = t.GetEntries(c + "&&" + spyB + cs);
         SB[k][b] = t.GetEntries(spyB + cs);
         double fc[6];
         int codes[5] = {1, 2, 3, 8, 9};
         double sumf = 0;
         for (int j = 0; j < 5; ++j) { fc[j] = t.GetEntries(c + Form("&&!sv&&sf==%d", codes[j])); sumf += fc[j]; }
         fc[5] = C[k][b] - V[k][b] - sumf;
         auto fr = [&](double x) { return C[k][b] > 0 ? x / C[k][b] : 0; };
         printf(" %2.0f-%2.0f | %5.0f  %.3f  %.3f | %.3f %.3f %.3f %.3f %.3f %.3f | %5.0f\n", TB[b], TB[b + 1], C[k][b],
                fr(V[k][b]), fr(VB[k][b]), fr(fc[0]), fr(fc[1]), fr(fc[2]), fr(fc[3]), fr(fc[4]), fr(fc[5]), SB[k][b]);
      }
   }
   int r = 5; // 35-40 reference
   printf("\n data/sim, normalised to 35-40:\n thcm   | pattern C   Spyral SB   (SB ratio)/(C ratio)\n");
   for (int b = 0; b < NT; ++b) {
      double rc = (C[1][b] > 0 && C[0][r] > 0) ? (C[0][b] / C[1][b]) / (C[0][r] / C[1][r]) : 0;
      double rs = (SB[1][b] > 0 && SB[0][r] > 0) ? (SB[0][b] / SB[1][b]) / (SB[0][r] / SB[1][r]) : 0;
      printf(" %2.0f-%2.0f |   %.3f       %.3f        %.3f\n", TB[b], TB[b + 1], rc, rs, rc > 0 ? rs / rc : 0);
   }
}
