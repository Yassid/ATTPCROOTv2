/// @file make_catima_table.C
/// @brief Generate a CATIMA dE/dx table in GenFit's SRIM-table format, so genfit
/// can be run with CATIMA energy loss (to isolate fitter vs eloss-model effects).
/// Also overlays CATIMA vs the SRIM deuteron table as a quick diagnostic.
///
/// genfit uses table value as MASS stopping [MeV/(mg/cm2)] and multiplies by
/// gasMediumDensity_[mg/cm3]. AtELossCATIMA::GetdEdx returns linear MeV/mm at the
/// CATIMA density. So table_value = GetdEdx[MeV/mm] / (density_gcm3 * 100); then
/// table_value * gasMediumDensity_[mg/cm3=density_gcm3*1000] = GetdEdx*10 = MeV/cm.
///
///   root -b -q 'pd/make_catima_table.C'

// A,Z,mass are arguments now: the (p,t) channel needs a TRITON curve in the same H2 gas, and the
// species was hard-coded. Defaults reproduce the deuteron table exactly, so existing calls are
// unaffected.
void make_catima_table(double densGcm3 = 8.3147e-5, // = 0.083147 mg/cm3 (genfit gasMediumDensity)
                       TString outFile = "../../../../resources/energy_loss/deuteron_H2_catima.txt",
                       int projA = 2, int projZ = 1, double projMass = 2.0135532)
{
   gSystem->Load("libAtTools.so");
   // deuteron in H2 gas
   auto eloss = new AtTools::AtELossCATIMA(densGcm3);
   eloss->SetProjectile(projA, projZ, projMass); // A,Z,mass(amu)
   std::vector<std::tuple<int, int, int>> mat;
   mat.push_back(std::make_tuple(1, 1, 1)); // hydrogen
   eloss->SetMaterial(mat);

   // energy grid (MeV): fine at low E, coarser high; up to 250 MeV (> data range)
   std::vector<double> egrid;
   for (double e = 0.01; e < 0.1; e += 0.01) egrid.push_back(e);
   for (double e = 0.1; e < 1.0; e += 0.05) egrid.push_back(e);
   for (double e = 1.0; e < 10.0; e += 0.25) egrid.push_back(e);
   for (double e = 10.0; e <= 250.0; e += 1.0) egrid.push_back(e);

   FILE *f = fopen(outFile.Data(), "w");
   fprintf(f, "   Ion        dE/dx      dE/dx     Projected  Longitudinal   Lateral\n");
   fprintf(f, "  Energy      Elec.      Nuclear     Range     Straggling   Straggling\n");
   fprintf(f, "-----------  ---------- ---------- ----------  ----------  ----------\n");

   TGraph *gC = new TGraph();
   for (double E : egrid) {
      double dedx_mm = eloss->GetdEdx(E);        // MeV/mm (linear, at densGcm3)
      double table = dedx_mm / (densGcm3 * 100); // MeV/(mg/cm2) mass stopping
      // genfit reads: ener enerUnit dEdx_elec dEdx_nucl range rangeUnit lonStra u latStra u
      fprintf(f, "  %8.3f MeV  %.4E  %.4E  %8.3f mm   %6.2f um   %6.2f um  \n", E, table, 0.0, 1.0, 0.0, 0.0);
      gC->SetPoint(gC->GetN(), E, table);
   }
   fclose(f);
   printf("wrote %s (%d points)\n", outFile.Data(), (int)egrid.size());

   // --- diagnostic overlay vs SRIM deuteron_D2_1bar.txt (mass stopping = elec+nucl) ---
   TString srim = "../../../../resources/energy_loss/deuteron_D2_1bar.txt";
   TGraph *gS = new TGraph();
   std::ifstream in(srim.Data());
   std::string line;
   for (int i = 0; i < 3 && std::getline(in, line); ++i) {}
   double en, de, dn, rg, l1, l2;
   std::string eu, ru, u1, u2;
   while (std::getline(in, line)) {
      std::istringstream ss(line);
      if (!(ss >> en >> eu >> de >> dn >> rg >> ru >> l1 >> u1 >> l2 >> u2)) continue;
      if (eu.find("keV") != std::string::npos) en /= 1000.0;
      gS->SetPoint(gS->GetN(), en, de + dn);
   }
   gStyle->SetOptStat(0);
   TCanvas *c = new TCanvas("c", "eloss", 800, 600);
   c->SetLogx(); c->SetLogy();
   gC->SetLineColor(kRed + 1); gC->SetLineWidth(3);
   gC->SetTitle("Deuteron dE/dx: CATIMA (H2) vs SRIM (D2);E [MeV];mass stopping [MeV/(mg/cm^{2})]");
   gC->Draw("AL");
   gS->SetLineColor(kAzure + 2); gS->SetLineWidth(3); gS->SetLineStyle(2); gS->Draw("L same");
   TLegend *lg = new TLegend(0.55, 0.72, 0.88, 0.88);
   lg->AddEntry(gC, "CATIMA (deuteron, H2)", "l");
   lg->AddEntry(gS, "SRIM (deuteron, D2 1bar)", "l");
   lg->Draw();
   c->SaveAs("pd/plots/eloss_catima_vs_srim.png");
   // print ratio at a few energies relevant to recoil deuterons (1-30 MeV)
   printf("\n  E[MeV]   CATIMA    SRIM(D2)   ratio\n");
   for (double E : {1.0, 2.0, 5.0, 10.0, 20.0, 30.0})
      printf("  %5.1f   %.4f   %.4f   %.3f\n", E, gC->Eval(E), gS->Eval(E), gC->Eval(E) / gS->Eval(E));
}
