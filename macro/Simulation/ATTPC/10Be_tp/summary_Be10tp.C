/// @file summary_Be10tp.C
/// @brief The 10Be(t,p)12Be field x pad-pitch campaign, in the binning this channel needs.
///
///   root -b -q 'summary_Be10tp.C("/mnt/f/Be10_tp")'
///
/// Adapted from summary_C14dp.C so the two channels are read the same way. Two things are
/// different and both matter:
///
/// 1. FOUR levels, not two, and they are close together: 0+ g.s., 2+ 2.109, 0+_2 2.251, 1- 2.715.
///    The pair that decides whether this experiment is worth doing is 2.109/2.251, only 142 keV
///    apart. That is far below anything the (d,p) study measured, so this campaign is asking a
///    much harder question of the same detector.
/// 2. The (t,p) proton is FAST. Where the (d,p) proton came out backward at 3-8 MeV, this one
///    runs 40 MeV at forward angles down to ~5 MeV past theta_lab 120. Where in theta_cm the
///    resolution is good is therefore not inherited from (d,p) and has to be measured.
///
/// Both reconstructions are reported, as in the (d,p) campaign: one constant beam energy for every
/// vertex (what the adopted analysis does) and the beam energy at the RECONSTRUCTED vertex, whose
/// z profile is solved out of MC truth here rather than assumed.

#include <algorithm>
#include <map>
#include <vector>

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const int NL = 4;
static const char *LVL[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
static const double LEX[NL] = {0.0, 2.109, 2.251, 2.715};
static const char *LJP[NL] = {"0+", "2+", "0+_2", "1-"};

// 10Be(t,p)12Be
static const double U = 931.49401;
static const double M1 = 10.0135341 * U;  // beam 10Be
static const double M2 = 3.0160493 * U;   // target t
static const double M3 = 1.007825 * U;    // ejectile p
static const double M4 = 12.0269221 * U;  // residual 12Be, ground state

static double tp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double tp_ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + M1, Et3 = Ke + M3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double a = (std::cos(th) * tp_om2(s, M1 * M1, M2 * M2) * tp_om2(uu, M2 * M2, M3 * M3) -
               (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                 (2 * M2 * M2) +
              s + uu - M2 * M2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
static double tp_q(std::vector<double> v, double p)
{
   if (v.size() < 20) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static TString tp_find(const TString &dir, const char *cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg, pre, lvl, cfg));
   return f.Strip(TString::kBoth);
}

struct Cell {
   bool ok{false};
   long n{0};
   double sConst{NAN}, mConst{NAN}, sZ{NAN}, mZ{NAN}, floor{NAN}, sKE{NAN}, sTh{NAN}, keMed{NAN}, thLabMed{NAN};
};

void summary_Be10tp(TString root = "/mnt/f/Be10_tp", Double_t Ebeam = 112.20)
{
   const int NS = 4;
   const double cmLo[NS] = {2, 20, 45, 90};
   const double cmHi[NS] = {20, 45, 90, 180};
   std::vector<std::map<TString, std::map<TString, Cell>>> S(NS); // [slice][cfg][level]

   for (int c = 0; c < NC; ++c)
      for (int l = 0; l < NL; ++l) {
         TString f = tp_find(root, CFG[c], LVL[l], "exres");
         if (f.IsNull()) continue;
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie()) continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t) { fr->Close(); continue; }
         double exReco, thTrue, thReco, keTrue, keReco, zTrue, zReco, cmTrue;
         t->SetBranchAddress("exReco", &exReco);
         t->SetBranchAddress("thTrue", &thTrue);
         t->SetBranchAddress("thReco", &thReco);
         t->SetBranchAddress("keTrue", &keTrue);
         t->SetBranchAddress("keReco", &keReco);
         t->SetBranchAddress("zTrue", &zTrue);
         t->SetBranchAddress("zReco", &zReco);
         t->SetBranchAddress("cmTrue", &cmTrue);
         const double ex0 = LEX[l], mres = M4 + ex0;

         // The beam-energy-vs-z profile this sample implies, solved from TRUTH. Fitted over all
         // angles, because the profile is a property of the gas and the vertex distribution, not
         // of an angular slice -- taking it per slice would just make it noisier.
         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = 80., hi = 130.;
            double flo = tp_ex(mres, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = tp_ex(mres, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double m = 0.5 * (lo + hi), fm = tp_ex(mres, m, thTrue * TMath::DegToRad(), keTrue);
               if (std::isnan(fm)) break;
               if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
            }
            double e = 0.5 * (lo + hi);
            if (e > 85 && e < 125) { eb.push_back(e); zz.push_back(zTrue); }
         }
         if (eb.size() < 100) { fr->Close(); continue; }
         TGraph g((int)eb.size(), zz.data(), eb.data());
         TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
         fEb.SetParameters(Ebeam, -0.006, 0.);
         g.Fit(&fEb, "QN");

         for (int s = 0; s < NS; ++s) {
            std::vector<double> dC, dZ, dF, dK, dT, kev, thv;
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               if (cmTrue < cmLo[s] || cmTrue >= cmHi[s]) continue;
               dC.push_back(exReco - ex0);
               dK.push_back(keReco - keTrue);
               dT.push_back(thReco - thTrue);
               kev.push_back(keTrue);
               thv.push_back(thTrue);
               double ebR = fEb.Eval(zReco);
               double a = tp_ex(M4, ebR, thReco * TMath::DegToRad(), keReco);
               double b = tp_ex(M4, ebR, thTrue * TMath::DegToRad(), keTrue);
               if (!std::isnan(a)) dZ.push_back(a - ex0);
               if (!std::isnan(b)) dF.push_back(b - ex0);
            }
            Cell cc;
            cc.n = dC.size();
            cc.ok = cc.n >= 20;
            if (!cc.ok) { S[s][CFG[c]][LVL[l]] = cc; continue; }
            cc.sConst = (tp_q(dC, .75) - tp_q(dC, .25)) / 1.349;
            cc.mConst = tp_q(dC, .5);
            cc.sZ = (tp_q(dZ, .75) - tp_q(dZ, .25)) / 1.349;
            cc.mZ = tp_q(dZ, .5);
            cc.floor = (tp_q(dF, .75) - tp_q(dF, .25)) / 1.349;
            cc.sKE = (tp_q(dK, .75) - tp_q(dK, .25)) / 1.349;
            cc.sTh = (tp_q(dT, .75) - tp_q(dT, .25)) / 1.349;
            cc.keMed = tp_q(kev, .5);
            cc.thLabMed = tp_q(thv, .5);
            S[s][CFG[c]][LVL[l]] = cc;
         }
         fr->Close();
      }

   // acceptance, per configuration and level
   std::map<TString, std::map<TString, double>> acc;
   for (int c = 0; c < NC; ++c)
      for (int l = 0; l < NL; ++l) {
         acc[CFG[c]][LVL[l]] = NAN;
         TString fa = tp_find(root, CFG[c], LVL[l], "acceptance");
         if (fa.IsNull()) continue;
         TFile *f = TFile::Open(fa);
         if (!f || f->IsZombie()) continue;
         TH1D *g = nullptr, *r = nullptr;
         TIter nx(f->GetListOfKeys());
         while (auto *k = (TKey *)nx()) {
            TString n = k->GetName();
            if (n.BeginsWith("hGen_")) g = (TH1D *)f->Get(n);
            if (n.BeginsWith("hRec_")) r = (TH1D *)f->Get(n);
         }
         if (g && r && g->Integral() > 0) acc[CFG[c]][LVL[l]] = r->Integral() / g->Integral();
         f->Close();
      }

   printf("\n=========== 10Be(t,p)12Be field x pitch campaign : %s ===========\n", root.Data());
   printf("chi2/ndf < 5, truth-matched, Ebeam(const) = %.2f MeV at the vertex.\n", Ebeam);
   printf("Levels: 0+ 0.000 | 2+ 2.109 | 0+_2 2.251 | 1- 2.715   (S_n = 3.171, all bound)\n");

   for (int s = 0; s < NS; ++s) {
      double kem = NAN, thm = NAN;
      for (int c = 0; c < NC && std::isnan(kem); ++c)
         if (S[s].count(CFG[c]) && S[s][CFG[c]].count("gs") && S[s][CFG[c]]["gs"].ok) {
            kem = S[s][CFG[c]]["gs"].keMed;
            thm = S[s][CFG[c]]["gs"].thLabMed;
         }
      printf("\n---- theta_cm %.0f-%.0f deg   (proton: theta_lab ~ %.0f deg, KE ~ %.1f MeV) ----\n", cmLo[s], cmHi[s],
             thm, kem);
      printf("  %-12s %6s | %9s %9s | %9s %9s | %8s %8s\n", "config", "n", "sig const", "sig E(z)", "med const",
             "med E(z)", "sig(KE)", "sig(th)");
      for (int c = 0; c < NC; ++c) {
         auto &cc = S[s][CFG[c]]["gs"];
         if (!cc.ok) { printf("  %-12s %6s | %9s %9s | %9s %9s | %8s %8s\n", CFG[c], "-", "-", "-", "-", "-", "-", "-"); continue; }
         printf("  %-12s %6ld | %9.3f %9.3f | %+9.3f %+9.3f | %8.3f %8.3f\n", CFG[c], cc.n, cc.sConst, cc.sZ,
                cc.mConst, cc.mZ, cc.sKE, cc.sTh);
      }
   }

   // THE question this campaign exists to answer.
   printf("\n---- SEPARATION of the 2.109/2.251 doublet (142 keV):  dE/(sigma_2109 + sigma_2251) ----\n");
   printf("   a value of 1 means the two peaks are one sigma-sum apart, i.e. barely distinguishable;\n");
   printf("   the (d,p) campaign's 15C doublet (740 keV) reached 5-7 by this measure.\n");
   printf("  %-12s", "config");
   for (int s = 0; s < NS; ++s) printf("   %3.0f-%-3.0f const  E(z)", cmLo[s], cmHi[s]);
   printf("\n");
   for (int c = 0; c < NC; ++c) {
      printf("  %-12s", CFG[c]);
      for (int s = 0; s < NS; ++s) {
         auto &a = S[s][CFG[c]]["ex2109"];
         auto &b = S[s][CFG[c]]["ex2251"];
         if (!a.ok || !b.ok) { printf(" %19s", "-"); continue; }
         printf(" %10.2f %8.2f", 0.142 / (a.sConst + b.sConst), 0.142 / (a.sZ + b.sZ));
      }
      printf("\n");
   }

   printf("\n---- separation of g.s. from the 2.109 (2109 keV), the easy case ----\n");
   printf("  %-12s", "config");
   for (int s = 0; s < NS; ++s) printf("   %3.0f-%-3.0f const  E(z)", cmLo[s], cmHi[s]);
   printf("\n");
   for (int c = 0; c < NC; ++c) {
      printf("  %-12s", CFG[c]);
      for (int s = 0; s < NS; ++s) {
         auto &a = S[s][CFG[c]]["gs"];
         auto &b = S[s][CFG[c]]["ex2109"];
         if (!a.ok || !b.ok) { printf(" %19s", "-"); continue; }
         printf(" %10.2f %8.2f", 2.109 / (a.sConst + b.sConst), 2.109 / (a.sZ + b.sZ));
      }
      printf("\n");
   }

   printf("\n---- sigma(Ex) per level, all theta_cm, constant Ebeam / vertex Ebeam ----\n");
   printf("  %-12s", "config");
   for (int l = 0; l < NL; ++l) printf(" %16s", Form("%s %.3f", LJP[l], LEX[l]));
   printf("\n");
   for (int c = 0; c < NC; ++c) {
      printf("  %-12s", CFG[c]);
      for (int l = 0; l < NL; ++l) {
         // recombine the four slices for an all-angle figure
         double num = 0, den = 0, numZ = 0;
         for (int s = 0; s < NS; ++s) {
            auto &cc = S[s][CFG[c]][LVL[l]];
            if (!cc.ok) continue;
            num += cc.sConst * cc.n; numZ += cc.sZ * cc.n; den += cc.n;
         }
         if (den <= 0) printf(" %16s", "-");
         else printf(" %8.3f %7.3f", num / den, numZ / den);
      }
      printf("\n");
   }

   printf("\n---- overall acceptance ----\n  %-12s", "config");
   for (int l = 0; l < NL; ++l) printf(" %10s", LVL[l]);
   printf("\n");
   for (int c = 0; c < NC; ++c) {
      printf("  %-12s", CFG[c]);
      for (int l = 0; l < NL; ++l) {
         double v = acc[CFG[c]][LVL[l]];
         if (std::isnan(v)) printf(" %10s", "-"); else printf(" %10.3f", v);
      }
      printf("\n");
   }
   printf("\nsummary done\n\n");
}
