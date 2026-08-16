/// @file make_eloss_table.C
/// @brief Write a genfit-readable dE/dx table generated from CATIMA.
///
/// WHY THIS EXISTS. genfit's MaterialEffects::dEdx sends every step with beta*gamma < 0.05 to
/// dEdxParam(), and dEdxParam() returns a hard 0.0 when no curve is loaded for the species. With
/// no energy-loss file that is our configuration, so BELOW beta*gamma = 0.05 genfit applies NO
/// energy loss at all -- not an approximate one. For a triton beta*gamma = 0.05 is KE = 3.5 MeV,
/// and the 16C(d,t)15C low branch is 0.8-6 MeV, so most of the physics sits in the dead region.
///
/// Loading a table fixes exactly that, and ONLY that: the curve is consulted only below
/// beta*gamma = 0.05, with Bethe-Bloch above, which is the hybrid we want. It does NOT fix
/// multiple scattering (noiseCoulomb is Highland and never looks at dEdx_), and it fixes energy
/// straggling only partly (the Urban terms scale with dEdx_, the Bohr term does not).
///
/// FORMAT, as parsed by MaterialEffects::setEnergyLossFile: three header lines, then
///   ener unit dEdx_elec dEdx_nucl range unit lonStra unit latStra unit
/// Only columns 1, 3 and 4 are used; the curve stores dEdx_elec + dEdx_nucl and genfit later
/// multiplies by gasMediumDensity_ (mg/cm3), so the table must be in MeV/(mg/cm2).
/// AtELossCATIMA::GetdEdx returns MeV/mm with the density already folded in, hence
///   table = GetdEdx * 10 / density_mg_per_cm3.
///
/// matA IS THE ISOTOPE MASS AND MUST BE THE REAL ONE. Stopping power per g/cm2 goes as Z/A of
/// the target, so deuterium's 2.014 against a round 2 is a flat +0.70% on every point of the
/// table -- verified, the ratio to CATIMA truth was 1.0070 at every energy from 0.5 to 40 MeV.
/// It has to be a double, and it cannot go through AtELossCATIMA::SetMaterial's
/// vector<tuple<int,int,int>> overload, which would silently truncate it back to 2. Build the
/// catima::Material directly and use the Material overload instead.
///
///   root -b -q 'make_eloss_table.C("triton_D2_300torr.txt", 3, 1, 3.01550072, 6.61e-5, 2.014)'
void make_eloss_table(TString out = "triton_D2_300torr.txt", int projA = 3, int projZ = 1,
                      double projMassAmu = 3.01550072, double densityGCm3 = 6.61e-5, double matA = 2.014,
                      double keMin = 0.01, double keMax = 60.0, int nPts = 600)
{
   gSystem->Load("libAtTools.so");
   AtTools::AtELossCATIMA el(densityGCm3);
   el.SetProjectile(projA, projZ, projMassAmu);
   catima::Material mat;
   mat.add_element(matA, 1, 1); // (A, Z, stoichiometry); Z=1 for both H and D
   el.SetMaterial(mat);

   const double densMg = densityGCm3 * 1000.0; // mg/cm3
   std::ofstream f(out.Data());
   f << "CATIMA dE/dx table for genfit  (A=" << projA << " Z=" << projZ << " in matA=" << matA << ")\n";
   f << "density " << densityGCm3 << " g/cm3 = " << densMg << " mg/cm3\n";
   f << "Energy Unit  dEdx_elec  dEdx_nucl  Range Unit  lonStra Unit  latStra Unit\n";

   // log spacing: the interesting region is the sub-MeV end where genfit currently has nothing
   for (int i = 0; i < nPts; ++i) {
      double ke = keMin * std::pow(keMax / keMin, double(i) / (nPts - 1));
      double dedxMeVmm = el.GetdEdx(ke);              // MeV/mm, density folded in
      double table = dedxMeVmm * 10.0 / densMg;       // -> MeV/(mg/cm2)
      double rangeMM = el.GetRange(ke);               // mm
      double lon = el.GetRangeStraggling(ke);         // mm
      if (!std::isfinite(table) || table <= 0)
         continue;
      f << Form("%.6g MeV %.6g 0 %.6g mm %.6g mm %.6g mm\n", ke, table, rangeMM,
                std::isfinite(lon) ? lon : 0.0, 0.0);
   }
   f.close();
   printf("wrote %s\n", out.Data());
   // sanity: genfit will reconstruct dedx = table * densMg (MeV/cm)
   for (double ke : {0.5, 1.0, 2.0, 3.5, 10.0}) {
      double d = el.GetdEdx(ke);
      printf("  KE %5.2f MeV : CATIMA %8.4f MeV/mm   (range %7.1f mm)\n", ke, d, el.GetRange(ke));
   }
}
