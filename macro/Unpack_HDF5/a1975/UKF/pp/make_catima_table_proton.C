/// @file make_catima_table_proton.C
/// @brief Proton-in-H CATIMA dE/dx table in genfit SRIM-table format, the proton
/// analogue of pd/make_catima_table.C (deuteron). Lets genfit run material effects
/// for protons with a tabulated eloss instead of its internal Bethe-Bloch (which
/// returns degenerate fits for PDG 2212). Writes resources/energy_loss/proton_H2_catima.txt.
///
///   root -b -q 'pp/make_catima_table_proton.C'
///
/// Density matches genfit's gasMediumDensity for the a1975 H gas (same as the deuteron table).

void make_catima_table_proton(double densGcm3 = 8.3147e-5,
                              TString outFile = "../../../../resources/energy_loss/proton_H2_catima.txt")
{
   gSystem->Load("libAtTools.so");
   auto eloss = new AtTools::AtELossCATIMA(densGcm3);
   eloss->SetProjectile(1, 1, 1.00727646688); // A, Z, mass(amu) -- proton
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen
   eloss->SetMaterial(mat);

   std::vector<double> egrid;
   for (double e = 0.01; e < 0.1; e += 0.01) egrid.push_back(e);
   for (double e = 0.1; e < 1.0; e += 0.05) egrid.push_back(e);
   for (double e = 1.0; e < 10.0; e += 0.25) egrid.push_back(e);
   for (double e = 10.0; e <= 250.0; e += 1.0) egrid.push_back(e);

   FILE *f = fopen(outFile.Data(), "w");
   fprintf(f, "   Ion        dE/dx      dE/dx     Projected  Longitudinal   Lateral\n");
   fprintf(f, "  Energy      Elec.      Nuclear     Range     Straggling   Straggling\n");
   fprintf(f, "-----------  ---------- ---------- ----------  ----------  ----------\n");
   for (double E : egrid) {
      double dedx_mm = eloss->GetdEdx(E);        // MeV/mm at densGcm3
      double table = dedx_mm / (densGcm3 * 100); // MeV/(mg/cm2) mass stopping
      fprintf(f, "  %8.3f MeV  %.4E  %.4E  %8.3f mm   %6.2f um   %6.2f um  \n", E, table, 0.0, 1.0, 0.0, 0.0);
   }
   fclose(f);
   printf("wrote %s (%d points)\n", outFile.Data(), (int)egrid.size());
}
