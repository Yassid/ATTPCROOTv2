/// @file dp_plots_C17.C
/// @brief The 17C(d,p)18C proposal figure: acceptance, the excitation-energy spectrum, and the
///        resolution, for the four bound states of 18C.
///
///   root -b -q 'dp_plots_C17.C'
///   root -b -q 'dp_plots_C17.C("/mnt/f/C17dp/b285_attpc/","./plots/")'
///
/// PANELS
///   A  acceptance vs theta_cm, per level. Read from the acceptance_*.root the campaign already
///      wrote. The per-level split is not cosmetic: an excited level puts the same theta_cm at a
///      different lab angle and energy, so its acceptance is a different curve.
///   B  reconstructed Ex spectrum, the four levels OVERLAID and each normalised to unit area.
///      Overlaid and equal-area rather than summed, because with a flat generator the relative
///      heights would be an artefact of how many events each level was given, not physics. When
///      the DWBA calculations arrive they set the weights; the shapes here do not change.
///   C  the same spectrum with the vertex beam-energy correction applied, which is the free
///      software gain (see dp_summary_C17.C).
///   D  sigma(Ex) vs theta_lab over the FULL range, with the method floor beside it.
///
/// EVERY PANEL IS FRAMED ON ITS OWN DATA. A fixed axis chosen in advance hides the shape of
/// whatever actually came out, and has cost a wrong physics claim before.
#include <algorithm>
#include <vector>

static double p17_omega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}
static double p17_ex(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * p17_omega2(s, m1 * m1, m2 * m2) * p17_omega2(u, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                         (2 * m2 * m2) +
                      s + u - m2 * m2;
   return arg <= 0 ? -1e9 : std::sqrt(arg) - m4;
}
static double p17_iqr(std::vector<double> v)
{
   if (v.size() < 20)
      return -1;
   std::sort(v.begin(), v.end());
   return (v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())]) / 1.349;
}

void dp_plots_C17(TString dir = "/mnt/f/C17dp/b285_attpc/", TString outDir = "./plots/",
                  Double_t EbeamConst = 135.0, Double_t chi2Cut = 5.0, Double_t mBeamAmu = 17.0225787,
                  Double_t mTgtAmu = 2.0141018, Double_t mEjAmu = 1.007825, Double_t mResAmu = 18.0267519)
{
   gStyle->SetOptStat(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gSystem->mkdir(outDir, kTRUE);

   const double u = 931.49401;
   const double m1 = mBeamAmu * u, m2 = mTgtAmu * u, m3 = mEjAmu * u, m4 = mResAmu * u;

   const int nL = 4;
   const double ExGen[nL] = {0.0, 1.588, 2.515, 3.972};
   const char *stTag[nL] = {"gs", "ex1588", "ex2515", "ex3972"};
   const char *stName[nL] = {"0^{+} g.s.", "2^{+} 1.588", "(2^{+}) 2.515", "(2,3)^{+} 3.972"};
   const int col[nL] = {kBlack, kRed + 1, kAzure + 2, kGreen + 3};

   // E_beam(z), measured from truth on the ground state (same solve as dp_summary_C17.C)
   double ebz_a = EbeamConst, ebz_b = 0;
   bool haveEbz = false;

   const int nSl = 12;
   const double slLo[nSl] = {0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165};
   const double slHi[nSl] = {15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180};

   TH1D *hExC[nL] = {nullptr}, *hExV[nL] = {nullptr};
   TGraph *gAcc[nL] = {nullptr};
   TGraph *gSig[nL] = {nullptr}, *gFlo[nL] = {nullptr};

   for (int l = 0; l < nL; ++l) {
      TString found = gSystem->GetFromPipe(TString("ls -1 ") + dir + "exres_" + stTag[l] + "_s*_b285_attpc.root" +
                                           " 2>/dev/null | head -1");
      found = found.Strip(TString::kBoth);
      if (found.IsNull()) {
         printf("\033[1;31m  %-8s exres MISSING\033[0m\n", stTag[l]);
         continue;
      }
      TFile *f = TFile::Open(found);
      TTree *t = f ? (TTree *)f->Get("res") : nullptr;
      if (!t) {
         printf("\033[1;31m  %-8s no res tree\033[0m\n", stTag[l]);
         continue;
      }
      double exR, exT, thT, thR, keT, keR, cmT, c2n, zT, zR;
      t->SetBranchAddress("exReco", &exR);
      t->SetBranchAddress("exTrue", &exT);
      t->SetBranchAddress("thTrue", &thT);
      t->SetBranchAddress("thReco", &thR);
      t->SetBranchAddress("keTrue", &keT);
      t->SetBranchAddress("keReco", &keR);
      t->SetBranchAddress("cmTrue", &cmT);
      t->SetBranchAddress("chi2ndf", &c2n);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);

      if (!haveEbz && l == 0) { // solve E_beam(z) against the GENERATED Ex -- see dp_summary_C17.C
         TGraph g;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = EbeamConst - 40, hi = EbeamConst + 40;
            if ((p17_ex(m1, m2, m3, m4, lo, thT, keT) - ExGen[l]) *
                   (p17_ex(m1, m2, m3, m4, hi, thT, keT) - ExGen[l]) > 0)
               continue;
            for (int it = 0; it < 60; ++it) {
               const double mid = 0.5 * (lo + hi);
               if ((p17_ex(m1, m2, m3, m4, lo, thT, keT) - ExGen[l]) *
                      (p17_ex(m1, m2, m3, m4, mid, thT, keT) - ExGen[l]) <= 0)
                  hi = mid;
               else
                  lo = mid;
            }
            g.SetPoint(g.GetN(), zT, 0.5 * (lo + hi));
         }
         if (g.GetN() > 100) {
            TF1 lin("lin", "pol1");
            g.Fit(&lin, "QN");
            ebz_a = lin.GetParameter(0);
            ebz_b = lin.GetParameter(1);
            haveEbz = true;
         }
      }

      hExC[l] = new TH1D(TString("hExC") + stTag[l], "", 220, -2.0, 7.0);
      hExV[l] = new TH1D(TString("hExV") + stTag[l], "", 220, -2.0, 7.0);
      std::vector<double> rsl[nSl], fsl[nSl];
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (c2n >= chi2Cut)
            continue;
         hExC[l]->Fill(exR);
         const double exV = haveEbz ? p17_ex(m1, m2, m3, m4, ebz_a + ebz_b * zR, thR, keR) : exR;
         hExV[l]->Fill(exV);
         for (int s = 0; s < nSl; ++s)
            if (thT >= slLo[s] && thT < slHi[s]) {
               rsl[s].push_back(exR - ExGen[l]);
               fsl[s].push_back(exT); // exTrue IS the method-floor residual, see dp_summary_C17.C
            }
      }
      gSig[l] = new TGraph();
      gFlo[l] = new TGraph();
      for (int s = 0; s < nSl; ++s) {
         const double sm = p17_iqr(rsl[s]), sf = p17_iqr(fsl[s]);
         if (sm > 0)
            gSig[l]->SetPoint(gSig[l]->GetN(), 0.5 * (slLo[s] + slHi[s]), sm);
         if (sf > 0)
            gFlo[l]->SetPoint(gFlo[l]->GetN(), 0.5 * (slLo[s] + slHi[s]), sf);
      }

      // acceptance, from what the campaign already wrote
      TString af = gSystem->GetFromPipe(TString("ls -1 ") + dir + "acceptance_" + stTag[l] + "_s*_b285_attpc.root" +
                                        " 2>/dev/null | head -1");
      af = af.Strip(TString::kBoth);
      if (!af.IsNull()) {
         TFile *fa = TFile::Open(af);
         if (fa && !fa->IsZombie()) {
            // acceptance_C14.C:207 writes it as hAcc_<tag>, alongside hGen_/hRec_/hLab_/hLabR_.
            // Match the prefix explicitly rather than by a substring: a loose match would silently
            // pick up a different histogram if that file ever gains one.
            TH1 *ha = nullptr;
            TIter nx(fa->GetListOfKeys());
            while (auto *k = (TKey *)nx()) {
               if (!TString(k->GetName()).BeginsWith("hAcc_"))
                  continue;
               TObject *o = k->ReadObj();
               if (o->InheritsFrom("TH1")) {
                  ha = (TH1 *)o;
                  break;
               }
            }
            if (!ha)
               printf("\033[1;31m  %-8s no hAcc_ histogram in %s\033[0m\n", stTag[l], af.Data());
            if (ha) {
               gAcc[l] = new TGraph();
               for (int b = 1; b <= ha->GetNbinsX(); ++b)
                  if (ha->GetBinContent(b) > 0)
                     gAcc[l]->SetPoint(gAcc[l]->GetN(), ha->GetBinCenter(b), ha->GetBinContent(b));
            }
         }
      }
   }

   TCanvas *c = new TCanvas("c17dp", "17C(d,p)18C", 1400, 1000);
   c->Divide(2, 2);

   // --- A: acceptance
   c->cd(1);
   gPad->SetGridy();
   bool first = true;
   for (int l = 0; l < nL; ++l) {
      if (!gAcc[l] || gAcc[l]->GetN() == 0)
         continue;
      gAcc[l]->SetLineColor(col[l]);
      gAcc[l]->SetMarkerColor(col[l]);
      gAcc[l]->SetMarkerStyle(20 + l);
      gAcc[l]->SetLineWidth(2);
      if (first) {
         gAcc[l]->SetTitle("acceptance;#theta_{cm} [deg];acceptance");
         gAcc[l]->GetYaxis()->SetRangeUser(0, 1.05);
         gAcc[l]->Draw("ALP");
         first = false;
      } else
         gAcc[l]->Draw("LP SAME");
   }
   {
      TLegend *lg = new TLegend(0.55, 0.15, 0.88, 0.42);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      for (int l = 0; l < nL; ++l)
         if (gAcc[l] && gAcc[l]->GetN())
            lg->AddEntry(gAcc[l], stName[l], "lp");
      lg->Draw();
   }

   // --- B/C: Ex spectra, equal area
   auto drawSpec = [&](TH1D *h[nL], const char *title) {
      double mx = 0;
      for (int l = 0; l < nL; ++l)
         if (h[l] && h[l]->Integral() > 0) {
            h[l]->Scale(1.0 / h[l]->Integral());
            mx = std::max(mx, h[l]->GetMaximum());
         }
      bool f1 = true;
      for (int l = 0; l < nL; ++l) {
         if (!h[l] || h[l]->Integral() <= 0)
            continue;
         h[l]->SetLineColor(col[l]);
         h[l]->SetLineWidth(2);
         if (f1) {
            h[l]->SetTitle(TString(title) + ";E_{x}(^{18}C) reconstructed [MeV];normalised counts");
            h[l]->GetYaxis()->SetRangeUser(0, 1.25 * mx);
            h[l]->Draw("HIST");
            f1 = false;
         } else
            h[l]->Draw("HIST SAME");
      }
      for (int l = 0; l < nL; ++l) { // where each level actually is
         TLine *ln = new TLine(ExGen[l], 0, ExGen[l], 1.25 * mx);
         ln->SetLineColor(col[l]);
         ln->SetLineStyle(2);
         ln->Draw();
      }
   };
   c->cd(2);
   drawSpec(hExC, "E_{x} spectrum, constant E_{beam}");
   c->cd(3);
   drawSpec(hExV, "E_{x} spectrum, E_{beam} at the reconstructed vertex");

   // --- D: sigma(Ex) vs theta_lab, with the floor
   c->cd(4);
   gPad->SetGridy();
   double smx = 0;
   for (int l = 0; l < nL; ++l)
      if (gSig[l])
         for (int i = 0; i < gSig[l]->GetN(); ++i)
            smx = std::max(smx, gSig[l]->GetY()[i]);
   first = true;
   for (int l = 0; l < nL; ++l) {
      if (!gSig[l] || gSig[l]->GetN() == 0)
         continue;
      gSig[l]->SetLineColor(col[l]);
      gSig[l]->SetMarkerColor(col[l]);
      gSig[l]->SetMarkerStyle(20 + l);
      gSig[l]->SetLineWidth(2);
      if (first) {
         gSig[l]->SetTitle("#sigma(E_{x}) vs #theta_{lab}  (dashed = method floor)"
                           ";#theta_{lab} true [deg];#sigma(E_{x}) [MeV]  (IQR/1.349)");
         gSig[l]->GetYaxis()->SetRangeUser(0, 1.15 * smx);
         gSig[l]->GetXaxis()->SetLimits(0, 180);
         gSig[l]->Draw("ALP");
         first = false;
      } else
         gSig[l]->Draw("LP SAME");
      if (gFlo[l] && gFlo[l]->GetN()) {
         gFlo[l]->SetLineColor(col[l]);
         gFlo[l]->SetLineStyle(2);
         gFlo[l]->SetLineWidth(2);
         gFlo[l]->Draw("L SAME");
      }
   }

   c->SaveAs(outDir + "C17dp_summary.png");
   printf("\nwrote %sC17dp_summary.png\n", outDir.Data());
   if (haveEbz)
      printf("E_beam(z) = %.4f %+.6f * z[mm], measured from truth\n\n", ebz_a, ebz_b);
}
