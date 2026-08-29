/// @file dp_budget_C14.C
/// @brief Why is the 14C(d,p) excitation-energy resolution what it is? Decompose it, event by
/// event, into the three things Ex is built from.
///
///   root -b -q 'dp_budget_C14.C()'
///
/// Ex is a function of the ejectile energy, the ejectile angle, and the beam energy at the vertex.
/// Each contributes  (dEx/dx) * (x_reco - x_true),  with the derivative evaluated per event on the
/// true kinematics. Summing the three should reproduce the actual Ex residual; that it does is the
/// check that the decomposition is complete rather than a story.
///
/// The beam term uses the difference between the CONSTANT beam energy the analysis assumes and the
/// true energy at that event's vertex, taken from the energy-loss profile -- so it is the cost of
/// the assumption, not of any measurement.

#include <algorithm>
#include <vector>

static const double U = 931.49401;
static const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;

static double bd_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double bd_ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + M1, Et3 = Ke + M3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double a = (std::cos(th) * bd_om2(s, M1 * M1, M2 * M2) * bd_om2(uu, M2 * M2, M3 * M3) -
               (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                 (2 * M2 * M2) +
              s + uu - M2 * M2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
static double bd_q(std::vector<double> v, double p)
{
   if (v.size() < 15) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static double bd_sig(std::vector<double> v) { return (bd_q(v, .75) - bd_q(v, .25)) / 1.349; }

void dp_budget_C14(TString root = "/mnt/f/a1954_C14dp_hf", Double_t Ebeam = 159.75, Double_t cmLo = 8,
                   Double_t cmHi = 30)
{
   const int NC = 6;
   const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};

   printf("\n  14C(d,p)15C, theta_cm %.0f-%.0f deg (the transfer peak). Where the E_x width comes from.\n", cmLo,
          cmHi);
   printf("  All entries are IQR/1.349 in MeV unless marked.\n\n");
   printf("  %-12s %5s | %8s %8s | %8s %8s %8s | %8s %8s | %8s\n", "config", "n", "dEx/dKE", "sig(KE)", "KE term",
          "th term", "beam term", "quad sum", "actual", "bias");
   for (int c = 0; c < NC; ++c) {
      std::vector<double> dKE, dTH, dEB, tot, act, derK;
      for (const char *lv : {"gs", "ex0740"}) {
         double ex0 = (TString(lv) == "gs") ? 0.0 : 0.740;
         TString f = gSystem->GetFromPipe(TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1",
                                                          root.Data(), CFG[c], lv, CFG[c]));
         f = f.Strip(TString::kBoth);
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
         const double mres = M4 + ex0;
         // the true beam energy at the vertex, from this sample's own truth
         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = 100., hi = 200.;
            double flo = bd_ex(mres, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = bd_ex(mres, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double m = 0.5 * (lo + hi), fm = bd_ex(mres, m, thTrue * TMath::DegToRad(), keTrue);
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

         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (cmTrue < cmLo || cmTrue >= cmHi) continue;
            double thR = thTrue * TMath::DegToRad();
            // derivatives on the TRUE kinematics of this event
            double dK = (bd_ex(M4, Ebeam, thR, keTrue + 0.05) - bd_ex(M4, Ebeam, thR, keTrue - 0.05)) / 0.1;
            double dT = (bd_ex(M4, Ebeam, thR + 0.1 * TMath::DegToRad(), keTrue) -
                         bd_ex(M4, Ebeam, thR - 0.1 * TMath::DegToRad(), keTrue)) / 0.2;
            double dB = (bd_ex(M4, Ebeam + 0.5, thR, keTrue) - bd_ex(M4, Ebeam - 0.5, thR, keTrue));
            if (std::isnan(dK) || std::isnan(dT) || std::isnan(dB)) continue;
            double tK = dK * (keReco - keTrue);
            double tT = dT * (thReco - thTrue);
            double tB = dB * (Ebeam - fEb.Eval(zTrue)); // cost of assuming a constant beam energy
            dKE.push_back(tK); dTH.push_back(tT); dEB.push_back(tB);
            tot.push_back(tK + tT + tB);
            act.push_back(exReco - ex0);
            derK.push_back(dK);
         }
         fr->Close();
      }
      if (dKE.size() < 15) { printf("  %-12s %5zu | (too few)\n", CFG[c], dKE.size()); continue; }
      std::vector<double> sK = dKE, sT = dTH, sB = dEB;
      double a = bd_sig(sK), b = bd_sig(sT), d = bd_sig(sB);
      printf("  %-12s %5zu | %8.2f %8.3f | %8.3f %8.3f %8.3f | %8.3f %8.3f | %+8.3f\n", CFG[c], dKE.size(),
             bd_q(derK, .5), bd_sig(*(new std::vector<double>(dKE))) / std::fabs(bd_q(derK, .5)), a, b, d,
             std::sqrt(a * a + b * b + d * d), bd_sig(act), bd_q(act, .5));
   }
   printf("\n  KE term = dEx/dKE x (KE_reco - KE_true);  th term likewise;  beam term = dEx/dE_beam x\n"
          "  (159.75 - E_true(z_vertex)), i.e. the cost of the constant-beam-energy assumption.\n"
          "  'quad sum' against 'actual' is the check that nothing else is contributing.\n\n");
}
