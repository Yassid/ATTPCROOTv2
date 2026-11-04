/// @file embed_vz_ana.C
/// Summary of embed_vz_dp.C: survival of sim protons overlaid on real data events per (vz, theta_cm), closure, how the lost ones fail, and survival against the observed data/sim band.
///   root -l -b -q embed_vz_ana.C
namespace eva { const int NT = 5, NZ = 7; const double TB[NT + 1] = {10, 20, 25, 30, 35, 40}, ZB[NZ + 1] = {200, 300, 400, 500, 600, 700, 800, 900};
/// observed, genfit frame B, data/sim normalised on 30-40 (genfit_vz_map_dp.C)
const double OBS[NZ][NT] = {{.38, 1.13, 2.10, 2.10, 1.29}, {.38, .86, 1.57, 1.41, 1.07}, {.20, .38, .49, .77, 1.07}, {.04, .17, .14, .46, 1.27}, {.08, .06, .04, .45, .64}, {.25, .10, .08, .65, .93}, {1.32, .11, .41, .96, .98}};
TString cell(int z, int t) { return Form("vz0>=%g&&vz0<%g&&tc0>=%g&&tc0<%g", ZB[z], ZB[z + 1], TB[t], TB[t + 1]); } }
void embed_vz_ana(TString pat = "/mnt/f/a1975/sim_rebuild/dp_embed/e_r*.root") { using namespace eva;
   TChain e("e"); e.Add(pat); printf("%lld rows\n", e.GetEntries());
   double S[NZ][NT], N[NZ][NT];
   printf("\nclosure, mode 1, proton hits alone:\n"); for (int z = 0; z < NZ; ++z) { printf(" %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int t = 0; t < NT; ++t) { double n = e.GetEntries("mode==1&&" + cell(z, t)); printf(" %6.3f", n > 0 ? e.GetEntries("mode==1&&found&&" + cell(z, t)) / n : -1.); } printf("\n"); }
   printf("\nSURVIVAL, mode 2, overlaid on data; columns 10-20 20-25 25-30 30-35 35-40 (35-40 sim tracks end at the wall)\n"); for (int z = 0; z < NZ; ++z) { printf(" %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int t = 0; t < NT; ++t) { N[z][t] = e.GetEntries("mode==2&&" + cell(z, t)); S[z][t] = N[z][t] > 0 ? e.GetEntries("mode==2&&found&&" + cell(z, t)) / N[z][t] : 0; printf(" %6.3f", S[z][t]); } printf("\n"); }
   printf("\noverlays per cell:\n"); for (int z = 0; z < NZ; ++z) { printf(" %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int t = 0; t < NT; ++t) printf(" %6.0f", N[z][t]); printf("\n"); }
   printf("\nobserved data/sim relative to the front rows (vz 200-400) of the same column, then survival on the same scale:\n"); for (int z = 0; z < NZ; ++z) { printf(" %3.0f-%3.0f", ZB[z], ZB[z + 1]); for (int t = 0; t < NT; ++t) printf("  %5.2f/%5.2f", OBS[z][t] / (0.5 * (OBS[0][t] + OBS[1][t])), S[z][t] / (0.5 * (S[0][t] + S[1][t]))); printf("\n"); }
   const char *CL[6] = {"cleaner took >50%% of proton hits", "no track holds 50%% of proton hits", "merged: track >30%% foreign hits", "Spyral invalid", "valid, kinematics fail or not matched", "all lost"};
   TString LOST = "mode==2&&!found"; TString A = LOST + "&&nPc<0.5*nhP", B = LOST + "&&nPc>=0.5*nhP&&bestP<0.5*nPc", C = LOST + "&&nPc>=0.5*nhP&&bestP>=0.5*nPc&&bestN>1.3*bestP", D = LOST + "&&nPc>=0.5*nhP&&bestP>=0.5*nPc&&bestN<=1.3*bestP&&!bestV", E = LOST + "&&nPc>=0.5*nhP&&bestP>=0.5*nPc&&bestN<=1.3*bestP&&bestV"; TString CC[5] = {A, B, C, D, E};
   printf("\nhow the lost ones fail, fraction of all overlays in the group:\n%-40s", "");
   const double GL[3] = {200, 400, 800}; const char *GN[3] = {"vz 200-400", "vz 400-800", "vz 800-900"}; const double GH[3] = {400, 800, 900}; for (auto n : GN) printf(" %14s", n); printf("\n");
   for (int t = 0; t < NT; ++t) { printf("theta_cm %.0f-%.0f\n", TB[t], TB[t + 1]); for (int c = 0; c < 5; ++c) { printf("  %-38s", Form(CL[c])); for (int g = 0; g < 3; ++g) { TString w = Form("&&vz0>=%g&&vz0<%g&&tc0>=%g&&tc0<%g", GL[g], GH[g], TB[t], TB[t + 1]); double n = e.GetEntries("mode==2" + w); printf(" %14.3f", n > 0 ? e.GetEntries(CC[c] + w) / n : -1.); } printf("\n"); } }
   printf("\nSpyral failCode of invalid case (all lost, vz 400-800, theta_cm 20-35):\n"); for (int f : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}) { double n = e.GetEntries(D + Form("&&bestF==%d&&vz0>=400&&vz0<800&&tc0>=20&&tc0<35", f)); if (n > 0) printf("  code %2d: %.0f\n", f, n); }
   const double NB[7] = {0, 50, 150, 300, 600, 1200, 1e6}; printf("\nsurvival vs data-event size, share; vz 200-400 then vz 400-800, theta_cm 20-35:\n"); for (int j = 0; j < 6; ++j) { TString s = Form("mode==2&&nhD>=%g&&nhD<%g&&tc0>=20&&tc0<35", NB[j], NB[j + 1]); printf("  nhD %5.0f-%-7.0f", NB[j], NB[j + 1]); for (auto zw : {"&&vz0>=200&&vz0<400", "&&vz0>=400&&vz0<800"}) { double n = e.GetEntries(s + zw), a = e.GetEntries(TString("mode==2&&tc0>=20&&tc0<35") + zw); printf("   share %.3f found %.3f", n / a, n > 0 ? e.GetEntries(s + zw + "&&found") / n : 0.); } printf("\n"); }
}
