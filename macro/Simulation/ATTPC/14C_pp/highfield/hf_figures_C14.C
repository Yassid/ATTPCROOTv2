/// @file hf_figures_C14.C
/// @brief The two figures the campaign is for: the 14C multiplet as each configuration would
/// record it, under BOTH ways of reconstructing the excitation energy.
///
///   root -b -q 'hf_figures_C14.C("/mnt/f/a1954_C14_hf")'
///
/// Left column of each pad pair: one constant beam energy for every vertex, which is what the
/// adopted a1954 analysis does. Right: the beam energy at the reconstructed vertex, using the
/// energy-loss profile extracted from the simulation's own truth (see exres_ebeamz_C14.C).
///
/// The three levels are added with EQUAL WEIGHTS. The real relative cross sections are an input
/// this study does not have, and weighting by them would fold a physics assumption into what is
/// meant to be a statement about the detector.

#include <algorithm>
#include <map>
#include <vector>

static const std::vector<TString> FG_CFG = {"b285_attpc", "b285_2mm", "b400_attpc",
                                            "b400_2mm",   "b700_attpc", "b700_2mm"};
static const std::vector<TString> FG_LVL = {"ex6094", "ex6728", "ex7012"};
static double fg_ex0(const TString &l)
{
   if (l == "ex6094") return 6.094;
   if (l == "ex6728") return 6.728;
   if (l == "ex7012") return 7.012;
   if (l == "gs") return 0.0;
   if (l == "ex8317") return 8.317;
   return -1;
}
static double fg_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double fg_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(th) * fg_om2(s, m1 * m1, m2 * m2) * fg_om2(uu, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                   (2 * m2 * m2) +
                s + uu - m2 * m2;
   return arg > 0 ? std::sqrt(arg) - m4 : NAN;
}
static TString fg_find(const TString &dir, const TString &cfg, const TString &lvl)
{
   TString f = gSystem->GetFromPipe(TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(),
                                                    cfg.Data(), lvl.Data(), cfg.Data()));
   return f.Strip(TString::kBoth);
}

void hf_figures_C14(TString rootDir = "/mnt/f/a1954_C14_hf", Double_t thMin = 20., Double_t thMax = 90.,
                    Double_t Ebeam0 = 159.75, TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;

   // [cfg][method][level] ; method 0 = constant E_beam, 1 = E_beam(z_reco)
   std::map<TString, std::array<std::map<TString, TH1D *>, 2>> H;

   for (auto &cfg : FG_CFG)
      for (auto &lvl : FG_LVL) {
         TString f = fg_find(rootDir, cfg, lvl);
         if (f.IsNull()) continue;
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie()) continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t) { fr->Close(); continue; }
         double exReco, thTrue, thReco, keTrue, keReco, zTrue, zReco;
         t->SetBranchAddress("exReco", &exReco);
         t->SetBranchAddress("thTrue", &thTrue);
         t->SetBranchAddress("thReco", &thReco);
         t->SetBranchAddress("keTrue", &keTrue);
         t->SetBranchAddress("keReco", &keReco);
         t->SetBranchAddress("zTrue", &zTrue);
         t->SetBranchAddress("zReco", &zReco);
         const double ex0 = fg_ex0(lvl), m_res = m_C14 + ex0;

         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (thTrue < thMin || thTrue >= thMax) continue;
            double lo = 100., hi = 200.;
            double flo = fg_ex(m_C14, m_p, m_p, m_res, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = fg_ex(m_C14, m_p, m_p, m_res, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double mid = 0.5 * (lo + hi);
               double fm = fg_ex(m_C14, m_p, m_p, m_res, mid, thTrue * TMath::DegToRad(), keTrue);
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

         auto *h0 = new TH1D("h0_" + cfg + "_" + lvl, "", 260, 5.2, 7.8);
         auto *h1 = new TH1D("h1_" + cfg + "_" + lvl, "", 260, 5.2, 7.8);
         h0->SetDirectory(nullptr);
         h1->SetDirectory(nullptr);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (thTrue < thMin || thTrue >= thMax) continue;
            h0->Fill(exReco);
            double v = fg_ex(m_C14, m_p, m_p, m_C14, fEb.Eval(zReco), thReco * TMath::DegToRad(), keReco);
            if (!std::isnan(v)) h1->Fill(v);
         }
         if (h0->Integral() > 0) h0->Scale(1. / h0->Integral());
         if (h1->Integral() > 0) h1->Scale(1. / h1->Integral());
         H[cfg][0][lvl] = h0;
         H[cfg][1][lvl] = h1;
         fr->Close();
      }

   const int col[3] = {kBlue + 1, kGreen + 2, kRed + 1};
   const char *mname[2] = {"constant E_{beam} (as analysed)", "E_{beam}(z_{reco})"};
   for (int m = 0; m < 2; ++m) {
      auto *c = new TCanvas(TString::Format("mult%d", m), "mult", 1600, 950);
      c->Divide(3, 2);
      int ip = 0;
      for (auto &cfg : FG_CFG) {
         c->cd(++ip);
         if (!H.count(cfg) || H[cfg][m].empty()) continue;
         TH1D *hs = nullptr;
         double ymax = 0;
         for (auto &lvl : FG_LVL)
            if (H[cfg][m].count(lvl)) {
               if (!hs) { hs = (TH1D *)H[cfg][m][lvl]->Clone("s_" + cfg + TString::Format("_%d", m));
                          hs->SetDirectory(nullptr); }
               else hs->Add(H[cfg][m][lvl]);
            }
         for (auto &lvl : FG_LVL)
            if (H[cfg][m].count(lvl)) ymax = std::max(ymax, H[cfg][m][lvl]->GetMaximum());
         if (hs) ymax = std::max(ymax, hs->GetMaximum());
         bool any = false;
         int ic = 0;
         for (auto &lvl : FG_LVL) {
            if (!H[cfg][m].count(lvl)) { ++ic; continue; }
            auto *h = H[cfg][m][lvl];
            h->SetLineColor(col[ic % 3]);
            h->SetLineWidth(2);
            h->SetTitle(cfg + ", " + mname[m] + ";E_{x} [MeV];normalised");
            h->GetYaxis()->SetRangeUser(0, 1.2 * ymax);
            h->Draw(any ? "hist same" : "hist");
            any = true;
            ++ic;
         }
         if (hs && any) {
            hs->SetLineColor(kBlack);
            hs->SetLineWidth(3);
            hs->Draw("hist same");
         }
         // the three true level positions, so a bump can be read against where it should be
         for (auto &lvl : FG_LVL) {
            auto *l = new TLine(fg_ex0(lvl), 0, fg_ex0(lvl), 1.2 * ymax);
            l->SetLineStyle(2);
            l->SetLineColor(kGray + 1);
            l->Draw();
         }
      }
      c->SaveAs(outDir + (m ? "hf_multiplet_ebeamz.png" : "hf_multiplet_asanalysed.png"));
   }
   printf("\nwrote %shf_multiplet_asanalysed.png and %shf_multiplet_ebeamz.png\nfigures done\n\n", outDir.Data(),
          outDir.Data());
}
