/// @file far_cmp_dp.C -- far-end test (2026-10-12). Data (census_arc data_*.root, proton gate) vs the maxELoss 11.5 g.s.
/// sim (a1975_C16_dp_sim_far/far_s41*_census.root), both on SPYRAL kinematics, bound window, theta_lab 89-161, KE>1.2.
/// The sim generates uniform in beam energy loss (~uniform in vz) and isotropic in cm over 2-50, so sim counts per
/// (vz slice, thcm bin) are efficiency x solid angle. Prints, per cell, data, sim and data/sim normalised so the
/// 35-40 bin over vz 200-900 is 1 -- a flat 1 everywhere would mean data = sim shape in vz AND angle.
/// Then the far-end (vz 900-1000) forward (thcm 10-20) events: Ex centroid/width and fine vz, data vs sim.
void far_cmp_dp()
{
   TChain d("c"), s("c");
   d.Add("/mnt/f/a1975/census_arc/data_*.root");
   s.Add("/mnt/f/a1975_C16_dp_sim_far/far_s41*_census.root");
   TGraph G;
   {
      std::ifstream in("pid/proton_d2_dv1104.json");
      std::string str((std::istreambuf_iterator<char>(in)), {});
      auto p = str.find("vertices"); p = str.find('[', p) + 1;
      while (true) {
         auto a = str.find('[', p); if (a == std::string::npos) break;
         auto b = str.find(']', a); double x, y; sscanf(str.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &x, &y);
         G.AddPoint(x, y); p = b + 1; if (str[str.find_first_not_of(" \n\t,", p)] == ']') break;
      }
   }
   auto cg = new TCutG("gp", G.GetN(), G.GetX(), G.GetY()); cg->SetVarX("ssq"); cg->SetVarY("sbr");
   TString B = "sv&&sth>89&&sth<161&&ske>1.2&&sex>-1.5&&sex<0.729";
   TString BD = B + "&&gp";
   printf("data rows %lld, sim rows %lld\n", d.GetEntries(), s.GetEntries());
   const int NT = 7, NZ = 11;
   double TB[NT + 1] = {10, 15, 20, 25, 30, 35, 40, 45};
   double ZB[NZ + 1] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1000, 1050};
   double D[NZ][NT], S[NZ][NT];
   for (int z = 0; z < NZ; ++z)
      for (int t = 0; t < NT; ++t) {
         TString c = Form("&&svz>=%g&&svz<%g&&stcm>=%g&&stcm<%g", ZB[z], ZB[z + 1], TB[t], TB[t + 1]);
         D[z][t] = d.GetEntries(BD + c);
         S[z][t] = s.GetEntries(B + c);
      }
   double dr = 0, sr = 0;
   for (int z = 1; z < 8; ++z) { dr += D[z][5]; sr += S[z][5]; }
   double norm = sr > 0 ? dr / sr : 0;
   printf("normalisation: data/sim in thcm 35-40, vz 200-900 = %.4f\n", norm);
   for (int k = 0; k < 3; ++k) {
      printf("\n%s\n  vz        |", k == 0 ? "DATA counts" : (k == 1 ? "SIM counts" : "data / (sim x norm)"));
      for (int t = 0; t < NT; ++t) printf("  %2.0f-%2.0f ", TB[t], TB[t + 1]);
      printf("\n");
      for (int z = 0; z < NZ; ++z) {
         printf("  %4.0f-%4.0f |", ZB[z], ZB[z + 1]);
         for (int t = 0; t < NT; ++t) {
            if (k == 0) printf("  %6.0f", D[z][t]);
            else if (k == 1) printf("  %6.0f", S[z][t]);
            else if (S[z][t] > 0) printf("  %6.2f", D[z][t] / (S[z][t] * norm));
            else printf("     -  ");
         }
         printf("\n");
      }
   }
   // far-end forward: are they normal g.s. events?
   auto stat = [](TChain &t, TString c, const char *v) {
      t.Draw(v, c, "goff"); int n = t.GetSelectedRows(); if (n < 3) { printf("    %-6s n=%d\n", v, n); return; }
      std::vector<double> x(t.GetV1(), t.GetV1() + n); std::sort(x.begin(), x.end());
      printf("    %-6s n=%4d  p16 %8.3f  med %8.3f  p84 %8.3f\n", v, n, x[n * 16 / 100], x[n / 2], x[n * 84 / 100]);
   };
   TString WB = "sv&&sth>89&&sth<161&&ske>1.2&&sex>-3&&sex<3";
   for (auto zr : {std::pair<double, double>{200, 400}, {900, 1000}}) {
      TString c = Form("&&svz>=%g&&svz<%g&&stcm>=10&&stcm<20", zr.first, zr.second);
      printf("\nthcm 10-20, vz %.0f-%.0f, Ex -3..3 (wide, to see the peak):\n  DATA\n", zr.first, zr.second);
      for (auto v : {"sex", "svz", "ske", "sth"}) stat(d, WB + "&&gp" + c, v);
      printf("  SIM\n");
      for (auto v : {"sex", "svz", "ske", "sth"}) stat(s, WB + c, v);
   }
   printf("\nfine vz, thcm 10-20, bound window:  vz   data   sim\n");
   for (double z = 850; z < 1050; z += 10)
      printf("   %4.0f-%4.0f   %4lld   %4lld\n", z, z + 10, d.GetEntries(BD + Form("&&svz>=%g&&svz<%g&&stcm>=10&&stcm<20", z, z + 10)),
             s.GetEntries(B + Form("&&svz>=%g&&svz<%g&&stcm>=10&&stcm<20", z, z + 10)));
}
