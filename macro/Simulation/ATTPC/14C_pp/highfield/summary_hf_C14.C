/// @file summary_hf_C14.C
/// @brief The whole 14C(p,p') field x pitch campaign in one set of tables.
///
///   root -b -q 'summary_hf_C14.C("/mnt/f/a1954_C14_hf")'
///
/// Everything is quoted TWICE, under the two ways of reconstructing the excitation energy:
///
///   AS ANALYSED  -- one constant beam energy for every vertex, which is what the adopted a1954
///                   analysis does.
///   E_beam(z)    -- the beam energy evaluated at the RECONSTRUCTED vertex, using the energy-loss
///                   profile that exres_ebeamz_C14.C extracts from the simulation's own truth.
///
/// The two differ by a lot and the difference is the point: under the first, no configuration is
/// detector-limited and the matrix says almost nothing; under the second the detector is what is
/// left and the matrix says everything. Reporting only one of them would be misleading either
/// way -- the first alone hides the hardware, the second alone implies an analysis that does not
/// exist yet.
///
/// Median and IQR, never a walk-outward FWHM (46Ar: that estimator returned 0.765 vs 1.800 MeV
/// for two histograms differing by a constant shift of every entry).

#include <algorithm>
#include <map>
#include <vector>

static const std::vector<TString> SM_CFG = {"b285_attpc", "b285_2mm", "b400_attpc",
                                            "b400_2mm",   "b700_attpc", "b700_2mm"};
static const std::vector<TString> SM_LVL = {"gs", "ex6094", "ex6728", "ex7012", "ex8317"};
static double sm_ex0(const TString &l)
{
   if (l == "gs") return 0.0;
   if (l == "ex6094") return 6.094;
   if (l == "ex6728") return 6.728;
   if (l == "ex7012") return 7.012;
   if (l == "ex8317") return 8.317;
   return -1;
}

static double sm_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double sm_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(th) * sm_om2(s, m1 * m1, m2 * m2) * sm_om2(uu, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                   (2 * m2 * m2) +
                s + uu - m2 * m2;
   return arg > 0 ? std::sqrt(arg) - m4 : NAN;
}
static double sm_q(std::vector<double> v, double p)
{
   if (v.empty()) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static double sm_sig(std::vector<double> v)
{
   return v.size() < 20 ? NAN : (sm_q(v, .75) - sm_q(v, .25)) / 1.349;
}
static double sm_med(std::vector<double> v)
{
   return v.size() < 20 ? NAN : sm_q(v, .50);
}
static TString sm_find(const TString &dir, const TString &cfg, const TString &lvl, const TString &pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg.Data(), pre.Data(),
                      lvl.Data(), cfg.Data()));
   return f.Strip(TString::kBoth);
}

struct SmCell {
   bool ok{false};
   long n{0};
   double sigConst{NAN}, medConst{NAN}, sigZ{NAN}, medZ{NAN}, floorZ{NAN};
   double sigKE{NAN}, sigTh{NAN}, acc{NAN}, eff{NAN};
};

void summary_hf_C14(TString rootDir = "/mnt/f/a1954_C14_hf", Double_t thMin = 20., Double_t thMax = 90.,
                    Double_t Ebeam0 = 159.75)
{
   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;
   std::map<TString, std::map<TString, SmCell>> S;

   for (auto &cfg : SM_CFG)
      for (auto &lvl : SM_LVL) {
         TString f = sm_find(rootDir, cfg, lvl, "exres");
         if (f.IsNull()) continue;
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie()) continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t) { fr->Close(); continue; }
         double exReco, exTrue, thTrue, thReco, keTrue, keReco, zTrue, zReco;
         t->SetBranchAddress("exReco", &exReco);
         t->SetBranchAddress("exTrue", &exTrue);
         t->SetBranchAddress("thTrue", &thTrue);
         t->SetBranchAddress("thReco", &thReco);
         t->SetBranchAddress("keTrue", &keTrue);
         t->SetBranchAddress("keReco", &keReco);
         t->SetBranchAddress("zTrue", &zTrue);
         t->SetBranchAddress("zReco", &zReco);
         const double ex0 = sm_ex0(lvl), m_res = m_C14 + ex0;

         // pass 1: the beam energy each truth event requires, fitted against the true vertex
         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (thTrue < thMin || thTrue >= thMax) continue;
            double lo = 100., hi = 200.;
            double flo = sm_ex(m_C14, m_p, m_p, m_res, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = sm_ex(m_C14, m_p, m_p, m_res, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double mid = 0.5 * (lo + hi);
               double fm = sm_ex(m_C14, m_p, m_p, m_res, mid, thTrue * TMath::DegToRad(), keTrue);
               if (std::isnan(fm)) break;
               if (fm * flo <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
            }
            double e = 0.5 * (lo + hi);
            if (e > 105 && e < 195) { eb.push_back(e); zz.push_back(zTrue); }
         }
         if (eb.size() < 100) { fr->Close(); continue; }
         TGraph g((int)eb.size(), zz.data(), eb.data());
         TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
         fEb.SetParameters(Ebeam0, -0.012, 0.);
         g.Fit(&fEb, "QN");

         // pass 2: both reconstructions, plus the tracking metrics
         std::vector<double> dC, dZ, dF, dK, dT;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (thTrue < thMin || thTrue >= thMax) continue;
            dC.push_back(exReco - ex0);
            dK.push_back(keReco - keTrue);
            dT.push_back(thReco - thTrue);
            double ebR = fEb.Eval(zReco);
            double a = sm_ex(m_C14, m_p, m_p, m_C14, ebR, thReco * TMath::DegToRad(), keReco);
            double b = sm_ex(m_C14, m_p, m_p, m_C14, ebR, thTrue * TMath::DegToRad(), keTrue);
            if (!std::isnan(a)) dZ.push_back(a - ex0);
            if (!std::isnan(b)) dF.push_back(b - ex0);
         }
         SmCell c;
         c.ok = dC.size() >= 20;
         c.n = dC.size();
         c.sigConst = sm_sig(dC); c.medConst = sm_med(dC);
         c.sigZ = sm_sig(dZ);     c.medZ = sm_med(dZ);
         c.floorZ = sm_sig(dF);
         c.sigKE = sm_sig(dK);    c.sigTh = sm_sig(dT);

         TString fa = sm_find(rootDir, cfg, lvl, "acceptance");
         if (!fa.IsNull()) {
            TFile *fo = TFile::Open(fa);
            if (fo && !fo->IsZombie()) {
               TH1D *hg = nullptr, *hr = nullptr;
               TIter nx(fo->GetListOfKeys());
               while (auto *k = (TKey *)nx()) {
                  TString n = k->GetName();
                  if (n.BeginsWith("hGen_")) hg = (TH1D *)fo->Get(n);
                  if (n.BeginsWith("hRec_")) hr = (TH1D *)fo->Get(n);
               }
               if (hg && hr && hg->Integral() > 0) c.acc = hr->Integral() / hg->Integral();
               fo->Close();
            }
         }
         S[cfg][lvl] = c;
         fr->Close();
      }

   auto table = [&](const char *title, const char *unit, std::function<double(const SmCell &)> get) {
      printf("\n---- %s [%s] ----\n", title, unit);
      printf("  %-12s", "config");
      for (auto &l : SM_LVL) printf(" %9s", l.Data());
      printf("\n");
      for (auto &cfg : SM_CFG) {
         printf("  %-12s", cfg.Data());
         for (auto &l : SM_LVL) {
            double v = (S.count(cfg) && S[cfg].count(l) && S[cfg][l].ok) ? get(S[cfg][l]) : NAN;
            if (std::isnan(v)) printf(" %9s", "-");
            else printf(" %9.3f", v);
         }
         printf("\n");
      }
   };

   printf("\n============ 14C(p,p') field x pitch campaign : %s ============\n", rootDir.Data());
   printf("theta_lab %.0f-%.0f deg, chi2/ndf < 5, truth-matched fits\n", thMin, thMax);

   table("sigma(Ex), AS ANALYSED (one constant beam energy)", "MeV", [](const SmCell &c) { return c.sigConst; });
   table("centroid bias, AS ANALYSED", "MeV", [](const SmCell &c) { return c.medConst; });
   table("sigma(Ex), beam energy at the RECONSTRUCTED vertex", "MeV", [](const SmCell &c) { return c.sigZ; });
   table("centroid bias, E_beam(z_reco)", "MeV", [](const SmCell &c) { return c.medZ; });
   table("METHOD FLOOR: E_beam(z_reco) with perfect tracking", "MeV", [](const SmCell &c) { return c.floorZ; });
   table("proton energy resolution sigma(KE_reco - KE_true)", "MeV", [](const SmCell &c) { return c.sigKE; });
   table("proton angle resolution sigma(theta_reco - theta_true)", "deg", [](const SmCell &c) { return c.sigTh; });
   table("acceptance (fraction of truth reactions with a good fit)", "-", [](const SmCell &c) { return c.acc; });

   // separation, on the corrected reconstruction -- the constant-beam one is method-limited and
   // its separation numbers would be a statement about the analysis, not the detector
   struct P { const char *a, *b; };
   std::vector<P> pairs = {{"ex6094", "ex6728"}, {"ex6728", "ex7012"}, {"ex7012", "ex8317"}};
   for (int mode = 0; mode < 2; ++mode) {
      printf("\n---- separation dE/(sigma_a+sigma_b), %s  (>2 resolvable, <1 one bump) ----\n",
             mode ? "E_beam(z_reco)" : "AS ANALYSED");
      printf("  %-12s", "config");
      for (auto &p : pairs) printf("  %5.3f-%-5.3f", sm_ex0(p.a), sm_ex0(p.b));
      printf("\n");
      for (auto &cfg : SM_CFG) {
         printf("  %-12s", cfg.Data());
         for (auto &p : pairs) {
            bool have = S.count(cfg) && S[cfg].count(p.a) && S[cfg].count(p.b) && S[cfg][p.a].ok && S[cfg][p.b].ok;
            if (!have) { printf(" %13s", "-"); continue; }
            double sa = mode ? S[cfg][p.a].sigZ : S[cfg][p.a].sigConst;
            double sb = mode ? S[cfg][p.b].sigZ : S[cfg][p.b].sigConst;
            printf(" %13.2f", (sm_ex0(p.b) - sm_ex0(p.a)) / (sa + sb));
         }
         printf("\n");
      }
   }
   printf("\nsummary done\n\n");
}
