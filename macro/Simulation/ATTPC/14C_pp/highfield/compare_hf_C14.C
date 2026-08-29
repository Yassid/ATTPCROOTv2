/// @file compare_hf_C14.C
/// @brief Read the whole 14C(p,p') field x pad-pitch campaign and answer the two questions it
/// was run for: how much does each configuration narrow Ex, and does any of them separate the
/// 6.728 / 7.012 MeV pair.
///
///   root -b -q 'compare_hf_C14.C("/mnt/f/a1954_C14_hf")'
///
/// Reads whatever exists -- a configuration or a level that has not finished is skipped and said
/// to be missing, never silently treated as zero.
///
/// THE SEPARATION METRIC. For two levels dE apart, quote  dE / (sigma_1 + sigma_2)  with
/// sigma = IQR/1.349 of the reconstructed Ex of each. Above ~2 the pair is resolvable in a fit,
/// below ~1 it is one bump. This is a resolution statement only: it says nothing about whether
/// the yields are large enough, which is what the acceptance columns are for.

#include <algorithm>
#include <map>
#include <vector>

static std::vector<TString> hf_configs()
{
   return {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
}
static std::vector<TString> hf_levels()
{
   return {"gs", "ex6094", "ex6728", "ex7012", "ex8317"};
}
static double hf_levelEx(const TString &l)
{
   if (l == "gs") return 0.0;
   if (l == "ex6094") return 6.094;
   if (l == "ex6728") return 6.728;
   if (l == "ex7012") return 7.012;
   if (l == "ex8317") return 8.317;
   return -1;
}

/// the one file of this (level, config), found by pattern -- the seed is part of the name and
/// differs between fields by construction
static TString hf_find(const TString &dir, const TString &cfg, const TString &lvl, const TString &prefix)
{
   TString cmd = TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg.Data(),
                                 prefix.Data(), lvl.Data(), cfg.Data());
   TString f = gSystem->GetFromPipe(cmd);
   f = f.Strip(TString::kBoth);
   return f;
}

struct HfStat {
   long n{0};
   double med{0}, iqr{0}, sigma{0}, frac05{0}, acc{0};
   bool ok{false};
};

static double hf_quant(std::vector<double> &v, double q)
{
   if (v.empty()) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, q * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void compare_hf_C14(TString rootDir = "/mnt/f/a1954_C14_hf", Double_t thMin = 20., Double_t thMax = 90.,
                    TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   auto cfgs = hf_configs();
   auto lvls = hf_levels();
   std::map<TString, std::map<TString, HfStat>> S; // [cfg][level]
   std::map<TString, std::map<TString, TH1D *>> H; // Ex spectra, kept for the figures

   printf("\n================ 14C(p,p') field x pitch campaign : %s ================\n", rootDir.Data());
   printf("theta_lab window %.0f-%.0f deg\n", thMin, thMax);

   for (auto &cfg : cfgs) {
      for (auto &lvl : lvls) {
         TString f = hf_find(rootDir, cfg, lvl, "exres");
         if (f.IsNull())
            continue;
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie())
            continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t)
            continue;
         double exReco, thTrue;
         t->SetBranchAddress("exReco", &exReco);
         t->SetBranchAddress("thTrue", &thTrue);
         std::vector<double> d;
         TString hn = "hx_" + cfg + "_" + lvl;
         auto *h = new TH1D(hn, "", 300, -3, 13);
         h->SetDirectory(nullptr);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (thTrue < thMin || thTrue >= thMax)
               continue;
            d.push_back(exReco - hf_levelEx(lvl));
            h->Fill(exReco);
         }
         HfStat s;
         s.n = d.size();
         if (s.n >= 20) {
            double q25 = hf_quant(d, .25), q50 = hf_quant(d, .50), q75 = hf_quant(d, .75);
            s.med = q50;
            s.iqr = q75 - q25;
            s.sigma = s.iqr / 1.349;
            s.frac05 = 100.0 * std::count_if(d.begin(), d.end(), [](double x) { return std::fabs(x) < 0.5; }) / s.n;
            s.ok = true;
         }
         // overall acceptance, straight out of the acceptance file of the same sample
         TString fa = hf_find(rootDir, cfg, lvl, "acceptance");
         if (!fa.IsNull()) {
            TFile *faf = TFile::Open(fa);
            if (faf && !faf->IsZombie()) {
               // the histogram names carry the full sample tag (seed included), so find them by
               // prefix rather than reconstructing the name
               TIter next(faf->GetListOfKeys());
               TH1D *g = nullptr, *r = nullptr;
               while (auto *k = (TKey *)next()) {
                  TString n = k->GetName();
                  if (n.BeginsWith("hGen_")) g = (TH1D *)faf->Get(n);
                  if (n.BeginsWith("hRec_")) r = (TH1D *)faf->Get(n);
               }
               if (g && r && g->Integral() > 0)
                  s.acc = r->Integral() / g->Integral();
               faf->Close();
            }
         }
         S[cfg][lvl] = s;
         H[cfg][lvl] = h;
         fr->Close();
      }
   }

   // ---- table 1 : width and centroid per level -------------------------------------------------
   for (auto &lvl : lvls) {
      printf("\n---- level %s (Ex = %.3f MeV) ----\n", lvl.Data(), hf_levelEx(lvl));
      printf("  %-12s %8s %10s %9s %10s %9s %9s\n", "config", "n", "median", "IQR", "sigma", "|d|<0.5", "acc");
      for (auto &cfg : cfgs) {
         if (!S.count(cfg) || !S[cfg].count(lvl) || !S[cfg][lvl].ok) {
            printf("  %-12s %8s %10s %9s %10s %9s %9s\n", cfg.Data(), "-", "-", "-", "-", "-", "-");
            continue;
         }
         auto &s = S[cfg][lvl];
         printf("  %-12s %8ld %+10.3f %9.3f %10.3f %8.1f%% %9.3f\n", cfg.Data(), s.n, s.med, s.iqr, s.sigma, s.frac05,
                s.acc);
      }
   }

   // ---- table 2 : can the multiplet be separated ------------------------------------------------
   struct Pair {
      const char *a, *b;
   };
   std::vector<Pair> pairs = {{"ex6094", "ex6728"}, {"ex6728", "ex7012"}, {"ex7012", "ex8317"}};
   printf("\n---- separation  dE / (sigma_a + sigma_b)  (>2 resolvable, <1 one bump) ----\n");
   printf("  %-12s", "config");
   for (auto &p : pairs)
      printf(" %6.3f-%-6.3f", hf_levelEx(p.a), hf_levelEx(p.b));
   printf("\n");
   for (auto &cfg : cfgs) {
      printf("  %-12s", cfg.Data());
      for (auto &p : pairs) {
         if (!S[cfg].count(p.a) || !S[cfg].count(p.b) || !S[cfg][p.a].ok || !S[cfg][p.b].ok) {
            printf(" %13s", "-");
            continue;
         }
         double dE = hf_levelEx(p.b) - hf_levelEx(p.a);
         printf(" %13.2f", dE / (S[cfg][p.a].sigma + S[cfg][p.b].sigma));
      }
      printf("\n");
   }

   // ---- figure : sigma vs theta_lab, all configurations on one pad per level --------------------
   {
      auto *cs = new TCanvas("hfsig", "hfsig", 1500, 900);
      cs->Divide(3, 2);
      const int col[6] = {kBlack, kGray + 2, kBlue + 1, kAzure + 7, kRed + 1, kOrange + 7};
      int ip = 0;
      for (auto &lvl : lvls) {
         cs->cd(++ip);
         bool any = false;
         auto *leg = new TLegend(0.50, 0.62, 0.89, 0.89);
         leg->SetBorderSize(0);
         leg->SetFillStyle(0);
         for (size_t ic = 0; ic < cfgs.size(); ++ic) {
            TString f = hf_find(rootDir, cfgs[ic], lvl, "exres");
            if (f.IsNull())
               continue;
            TFile *fr = TFile::Open(f);
            if (!fr || fr->IsZombie())
               continue;
            TH1D *h = nullptr;
            TIter next(fr->GetListOfKeys());
            while (auto *k = (TKey *)next())
               if (TString(k->GetName()).BeginsWith("hIQR_"))
                  h = (TH1D *)fr->Get(k->GetName());
            if (!h)
               continue;
            h = (TH1D *)h->Clone(TString("sig_") + cfgs[ic] + "_" + lvl);
            h->SetDirectory(nullptr);
            h->Scale(1.0 / 1.349); // IQR -> gaussian-equivalent sigma
            h->SetLineColor(col[ic]);
            h->SetMarkerColor(col[ic]);
            h->SetMarkerStyle(20 + (int)ic % 4);
            h->SetLineWidth(2);
            h->SetTitle(lvl + " (E_{x} = " + TString::Format("%.3f", hf_levelEx(lvl)) +
                        " MeV);#theta_{lab} [deg];#sigma(E_{x}) = IQR/1.349 [MeV]");
            h->GetYaxis()->SetRangeUser(0, 1.2);
            h->Draw(any ? "e1 same" : "e1");
            leg->AddEntry(h, cfgs[ic], "lp");
            any = true;
            fr->Close();
         }
         if (any)
            leg->Draw();
      }
      cs->SaveAs(outDir + "hf_sigma_vs_theta.png");
      printf("wrote %shf_sigma_vs_theta.png\n", outDir.Data());
   }

   // ---- figure : the multiplet as each configuration would see it -------------------------------
   auto *c = new TCanvas("hfcmp", "hfcmp", 1500, 900);
   c->Divide(3, 2);
   int ipad = 0;
   for (auto &cfg : cfgs) {
      c->cd(++ipad);
      // normalise every level to unit area first, so the widths are comparable and the sum below
      // is not dominated by whichever sample happened to have more statistics
      for (auto &lvl : lvls)
         if (H[cfg].count(lvl) && H[cfg][lvl]->Integral() > 0)
            H[cfg][lvl]->Scale(1.0 / H[cfg][lvl]->Integral());

      // THE SUM, at equal weights. The individual curves say how wide each level is; only the
      // sum says whether an experiment looking at this spectrum would see one bump or two.
      // Equal weights are deliberate -- the real relative cross sections are an input this study
      // does not have, and weighting by them would fold a physics assumption into what is meant
      // to be a detector statement.
      TH1D *hSum = nullptr;
      for (auto &lvl : {TString("ex6094"), TString("ex6728"), TString("ex7012")}) {
         if (!H[cfg].count(lvl))
            continue;
         if (!hSum) {
            hSum = (TH1D *)H[cfg][lvl]->Clone("sum_" + cfg);
            hSum->SetDirectory(nullptr);
         } else
            hSum->Add(H[cfg][lvl]);
      }

      double ymax = 0;
      for (auto &lvl : lvls)
         if (H[cfg].count(lvl))
            ymax = std::max(ymax, H[cfg][lvl]->GetMaximum());
      if (hSum)
         ymax = std::max(ymax, hSum->GetMaximum());

      bool any = false;
      int ic = 0;
      const int col[5] = {kBlack, kBlue + 1, kGreen + 2, kRed + 1, kMagenta + 1};
      for (auto &lvl : lvls) {
         if (!H[cfg].count(lvl)) {
            ++ic;
            continue;
         }
         auto *h = H[cfg][lvl];
         h->SetLineColor(col[ic % 5]);
         h->SetLineWidth(2);
         h->SetTitle(cfg + ";E_{x} reconstructed [MeV];normalised");
         h->GetYaxis()->SetRangeUser(0, 1.15 * ymax);
         h->Draw(any ? "hist same" : "hist");
         any = true;
         ++ic;
      }
      if (hSum && any) {
         hSum->SetLineColor(kGray + 2);
         hSum->SetLineWidth(3);
         hSum->Draw("hist same");
      }
      if (!any) {
         auto *tx = new TLatex(0.2, 0.5, cfg + ": no data");
         tx->SetNDC();
         tx->Draw();
      }
   }
   c->SaveAs(outDir + "hf_multiplet.png");
   printf("\nwrote %shf_multiplet.png\n", outDir.Data());
   printf("\ncompare done\n\n");
}
