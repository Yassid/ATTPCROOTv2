/// @file summary_C14dp.C
/// @brief The 14C(d,p)15C field x pad-pitch campaign, in the binning the channel actually needs.
///
///   root -b -q 'summary_C14dp.C("/mnt/f/a1954_C14dp_hf")'
///
/// Unlike (p,p'), this channel's interesting region is BACKWARD in the laboratory: at theta_cm 20
/// deg the proton comes out at theta_lab 125 deg with 3.1 MeV. So everything here is sliced in
/// theta_cm rather than theta_lab, and the first slice -- where a transfer angular distribution
/// has its yield -- is the one to read.
///
/// Both reconstructions are reported, for the same reason as in the (p,p') campaign: one constant
/// beam energy for every vertex (what the adopted analysis does) and the beam energy at the
/// reconstructed vertex. Here the two are expected to differ a great deal, because dEx/dE_beam at
/// theta_cm 20 deg is ~0.047 against 0.004 for (p,p').

#include <algorithm>
#include <map>
#include <vector>

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const int NL = 2;
static const char *LVL[NL] = {"gs", "ex0740"};
static const double LEX[NL] = {0.0, 0.740};

// 14C(d,p)15C
static const double U = 931.49401;
static const double M1 = 14.003242 * U;   // beam
static const double M2 = 2.0141018 * U;   // target d
static const double M3 = 1.007825 * U;    // ejectile p
static const double M4 = 15.0105993 * U;  // residual 15C, ground state

static double dp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double dp_ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + M1, Et3 = Ke + M3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double a = (std::cos(th) * dp_om2(s, M1 * M1, M2 * M2) * dp_om2(uu, M2 * M2, M3 * M3) -
               (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                 (2 * M2 * M2) +
              s + uu - M2 * M2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
static double dp_q(std::vector<double> v, double p)
{
   if (v.size() < 20) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static TString dp_find(const TString &dir, const char *cfg, const char *lvl, const char *pre)
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

void summary_C14dp(TString root = "/mnt/f/a1954_C14dp_hf", Double_t Ebeam = 159.75)
{
   const int NS = 4;
   const double cmLo[NS] = {8, 30, 60, 100};
   const double cmHi[NS] = {30, 60, 100, 180};
   // [slice][cfg][level]
   std::vector<std::map<TString, std::map<TString, Cell>>> S(NS);

   for (int c = 0; c < NC; ++c)
      for (int l = 0; l < NL; ++l) {
         TString f = dp_find(root, CFG[c], LVL[l], "exres");
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

         // the beam-energy profile this sample requires, solved from truth (all angles, so the fit
         // has the statistics; the profile is a property of the gas, not of the slice)
         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = 100., hi = 200.;
            double flo = dp_ex(mres, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = dp_ex(mres, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double m = 0.5 * (lo + hi), fm = dp_ex(mres, m, thTrue * TMath::DegToRad(), keTrue);
               if (std::isnan(fm)) break;
               if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
            }
            double e = 0.5 * (lo + hi);
            if (e > 105 && e < 195) { eb.push_back(e); zz.push_back(zTrue); }
         }
         if (eb.size() < 100) { fr->Close(); continue; }
         TGraph g((int)eb.size(), zz.data(), eb.data());
         TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
         fEb.SetParameters(Ebeam, -0.010, 0.);
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
               double a = dp_ex(M4, ebR, thReco * TMath::DegToRad(), keReco);
               double b = dp_ex(M4, ebR, thTrue * TMath::DegToRad(), keTrue);
               if (!std::isnan(a)) dZ.push_back(a - ex0);
               if (!std::isnan(b)) dF.push_back(b - ex0);
            }
            Cell cc;
            cc.n = dC.size();
            cc.ok = cc.n >= 20;
            if (!cc.ok) { S[s][CFG[c]][LVL[l]] = cc; continue; }
            cc.sConst = (dp_q(dC, .75) - dp_q(dC, .25)) / 1.349;
            cc.mConst = dp_q(dC, .5);
            cc.sZ = (dp_q(dZ, .75) - dp_q(dZ, .25)) / 1.349;
            cc.mZ = dp_q(dZ, .5);
            cc.floor = (dp_q(dF, .75) - dp_q(dF, .25)) / 1.349;
            cc.sKE = (dp_q(dK, .75) - dp_q(dK, .25)) / 1.349;
            cc.sTh = (dp_q(dT, .75) - dp_q(dT, .25)) / 1.349;
            cc.keMed = dp_q(kev, .5);
            cc.thLabMed = dp_q(thv, .5);
            S[s][CFG[c]][LVL[l]] = cc;
         }
         fr->Close();
      }

   // acceptance, per configuration and level
   std::map<TString, std::map<TString, double>> acc;
   for (int c = 0; c < NC; ++c)
      for (int l = 0; l < NL; ++l) {
         acc[CFG[c]][LVL[l]] = NAN;
         TString fa = dp_find(root, CFG[c], LVL[l], "acceptance");
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

   printf("\n=========== 14C(d,p)15C field x pitch campaign : %s ===========\n", root.Data());
   printf("chi2/ndf < 5, truth-matched. Sliced in theta_cm (analysis convention): small theta_cm is\n"
          "the BACKWARD-lab proton, which is where a transfer angular distribution has its yield.\n");

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

   printf("\n---- separation of the 15C doublet, 0.740 MeV apart:  dE/(sigma_gs + sigma_0740) ----\n");
   printf("  %-12s", "config");
   for (int s = 0; s < NS; ++s) printf("   %3.0f-%-3.0f const  E(z)", cmLo[s], cmHi[s]);
   printf("\n");
   for (int c = 0; c < NC; ++c) {
      printf("  %-12s", CFG[c]);
      for (int s = 0; s < NS; ++s) {
         auto &a = S[s][CFG[c]]["gs"];
         auto &b = S[s][CFG[c]]["ex0740"];
         if (!a.ok || !b.ok) { printf(" %19s", "-"); continue; }
         printf(" %10.2f %8.2f", 0.740 / (a.sConst + b.sConst), 0.740 / (a.sZ + b.sZ));
      }
      printf("\n");
   }

   printf("\n---- overall acceptance ----\n  %-12s %10s %10s\n", "config", "gs", "0.740");
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
