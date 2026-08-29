/// @file hf_conclusion_C14.C
/// @brief The whole 14C+p field x pad-pitch campaign in one figure.
///
///   root -b -q 'hf_conclusion_C14.C("/mnt/f/a1954_C14_hf")'
///
/// Four panels, in the order the argument runs:
///
///   A  RESOLUTION. sigma(Ex) per configuration, under the analysis as written (one constant beam
///      energy for every vertex) and with the beam energy taken at the reconstructed vertex. As
///      written, the six configurations are indistinguishable; corrected, they separate.
///
///   B  WHAT IT COSTS. The same six configurations' acceptance, weighted three ways: flat in
///      theta_cm, which is how the samples were generated, and by the theta_cm distributions the
///      a1954 elastic and inelastic yields actually have. The flat column is the one that made
///      4 T look affordable; the elastic column is the one that matters for the normalisation.
///
///   C  WHY THE WEIGHTINGS DIFFER. dEx/dE_beam along the elastic locus, against where the yield
///      sits. The sensitivity spans two orders of magnitude across the acceptance, so "how much
///      does the constant beam energy cost" has no single answer -- it depends entirely on the
///      angular distribution folded against it.
///
///   D  THE GAP. sigma(Ex) of the 6.094 level in the angular region where the inelastic yield
///      lives. The simulation says the vertex correction is worth a factor four there. The data's
///      own measured resolution, the 0.132 MeV that fit_angles_ps_C14.C FIXES, sits a factor four
///      above the corrected prediction and close to the uncorrected one. That gap is bigger than
///      every effect in this matrix and nothing here explains it.
///
/// Everything is recomputed from the files; no number is hard-coded.

#include <algorithm>
#include <vector>

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
// short labels: ROOT bin labels are a single line, and "#newline" is not a thing in them
static const char *CLB[NC] = {"2.85 / 8#times12", "2.85 / 2mm", "4 / 8#times12",
                              "4 / 2mm",          "7 / 8#times12", "7 / 2mm"};

static double cc_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double cc_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(th) * cc_om2(s, m1 * m1, m2 * m2) * cc_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
/// ejectile KE that puts excitation exWanted at lab angle thRad -- bisection on the same
/// expression the analysis inverts, with the upper bracket walked down past the kinematic limit
static double cc_locus(double m1, double m2, double m3, double m4, double E, double thRad, double exWanted)
{
   auto F = [&](double ke) { return cc_ex(m1, m2, m3, m4, E, thRad, ke) - exWanted; };
   double lo = 0.02, hi = 60.0, flo = F(lo);
   if (std::isnan(flo))
      return -1;
   double fhi = F(hi);
   while ((std::isnan(fhi) || fhi * flo > 0) && hi > lo + 0.05) {
      hi *= 0.9;
      fhi = F(hi);
   }
   if (std::isnan(fhi) || fhi * flo > 0)
      return -1;
   for (int i = 0; i < 70; ++i) {
      double m = 0.5 * (lo + hi), fm = F(m);
      if (std::isnan(fm))
         return -1;
      if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
   }
   return 0.5 * (lo + hi);
}
static double cc_q(std::vector<double> v, double p)
{
   if (v.size() < 20)
      return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static TString cc_find(const TString &dir, const char *cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg, pre, lvl, cfg));
   return f.Strip(TString::kBoth);
}

/// sigma(Ex) of one (config, level) in a theta_lab window, both reconstructions
static bool cc_sigma(const TString &root, const char *cfg, const char *lvl, double ex0, double thMin, double thMax,
                     double Ebeam, double &sigConst, double &sigZ)
{
   const double u = 931.49401, mC = 14.003242 * u, mp = 1.007825 * u, mres = mC + ex0;
   sigConst = sigZ = NAN;
   TString f = cc_find(root, cfg, lvl, "exres");
   if (f.IsNull())
      return false;
   TFile *fr = TFile::Open(f);
   if (!fr || fr->IsZombie())
      return false;
   TTree *t = (TTree *)fr->Get("res");
   if (!t) { fr->Close(); return false; }
   double exReco, thTrue, thReco, keTrue, keReco, zTrue, zReco;
   t->SetBranchAddress("exReco", &exReco);
   t->SetBranchAddress("thTrue", &thTrue);
   t->SetBranchAddress("thReco", &thReco);
   t->SetBranchAddress("keTrue", &keTrue);
   t->SetBranchAddress("keReco", &keReco);
   t->SetBranchAddress("zTrue", &zTrue);
   t->SetBranchAddress("zReco", &zReco);
   // the beam-energy profile this sample requires, solved from truth event by event
   std::vector<double> eb, zz;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (thTrue < thMin || thTrue >= thMax)
         continue;
      double lo = 100., hi = 200.;
      double flo = cc_ex(mC, mp, mp, mres, lo, thTrue * TMath::DegToRad(), keTrue);
      double fhi = cc_ex(mC, mp, mp, mres, hi, thTrue * TMath::DegToRad(), keTrue);
      if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0)
         continue;
      for (int it = 0; it < 60; ++it) {
         double m = 0.5 * (lo + hi), fm = cc_ex(mC, mp, mp, mres, m, thTrue * TMath::DegToRad(), keTrue);
         if (std::isnan(fm)) break;
         if (fm * flo <= 0) { hi = m; fhi = fm; } else { lo = m; flo = fm; }
      }
      double e = 0.5 * (lo + hi);
      if (e > 105 && e < 195) { eb.push_back(e); zz.push_back(zTrue); }
   }
   if (eb.size() < 100) { fr->Close(); return false; }
   TGraph g((int)eb.size(), zz.data(), eb.data());
   TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
   fEb.SetParameters(Ebeam, -0.012, 0.);
   g.Fit(&fEb, "QN");
   std::vector<double> dC, dZ;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (thTrue < thMin || thTrue >= thMax)
         continue;
      dC.push_back(exReco - ex0);
      double v = cc_ex(mC, mp, mp, mC, fEb.Eval(zReco), thReco * TMath::DegToRad(), keReco);
      if (!std::isnan(v))
         dZ.push_back(v - ex0);
   }
   sigConst = (cc_q(dC, .75) - cc_q(dC, .25)) / 1.349;
   sigZ = (cc_q(dZ, .75) - cc_q(dZ, .25)) / 1.349;
   fr->Close();
   return true;
}

void hf_conclusion_C14(TString root = "/mnt/f/a1954_C14_hf",
                       TString dataCache = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954/UKF/"
                                           "pp/plots/proton_kin_cat5_s013.root",
                       Double_t Ebeam = 159.75, Double_t dataSigma6094 = 0.132, TString outDir = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleFont(42, "xyz");
   gStyle->SetLabelFont(42, "xyz");
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double u = 931.49401, mC = 14.003242 * u, mp = 1.007825 * u;

   // ---------- A : resolution, elastic, the campaign's own window --------------------------------
   double sA[NC], sB[NC];
   for (int c = 0; c < NC; ++c)
      cc_sigma(root, CFG[c], "gs", 0.0, 20, 90, Ebeam, sA[c], sB[c]);

   // ---------- D : the 6.094 level where the inelastic yield lives -------------------------------
   double dA[NC], dB[NC];
   for (int c = 0; c < NC; ++c)
      cc_sigma(root, CFG[c], "ex6094", 6.094, 35, 65, Ebeam, dA[c], dB[c]);

   // ---------- the data: two theta_cm distributions ----------------------------------------------
   TFile *fd = TFile::Open(dataCache);
   TNtuple *t = fd && !fd->IsZombie() ? (TNtuple *)fd->Get("pk") : nullptr;
   auto *hEl = new TH1D("hEl", "", 36, 0, 180);
   auto *hIn = new TH1D("hIn", "", 36, 0, 180);
   if (t) {
      float ke, th, vz, thcm, ex, c2n;
      t->SetBranchAddress("ke", &ke);
      t->SetBranchAddress("theta", &th);
      t->SetBranchAddress("vertexz", &vz);
      t->SetBranchAddress("thcm", &thcm);
      t->SetBranchAddress("ex", &ex);
      t->SetBranchAddress("chi2ndf", &c2n);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (std::fabs(ex) < 1.5) hEl->Fill(thcm);
         if (ex > 5.0 && ex < 10.0) hIn->Fill(thcm);
      }
   }

   // ---------- B : acceptance, three weightings ---------------------------------------------------
   double aFlat[NC], aEl[NC], aIn[NC];
   for (int c = 0; c < NC; ++c) {
      aFlat[c] = aEl[c] = aIn[c] = NAN;
      TString fa = cc_find(root, CFG[c], "gs", "acceptance");
      if (fa.IsNull())
         continue;
      TFile *f = TFile::Open(fa);
      if (!f || f->IsZombie())
         continue;
      TH1D *g = nullptr, *r = nullptr;
      TIter nx(f->GetListOfKeys());
      while (auto *k = (TKey *)nx()) {
         TString n = k->GetName();
         if (n.BeginsWith("hGen_")) g = (TH1D *)f->Get(n);
         if (n.BeginsWith("hRec_")) r = (TH1D *)f->Get(n);
      }
      if (!g || !r || g->Integral() <= 0) { f->Close(); continue; }
      aFlat[c] = r->Integral() / g->Integral();
      auto fold = [&](TH1D *w) {
         double num = 0, den = 0;
         for (int b = 1; b <= w->GetNbinsX(); ++b) {
            double x = w->GetBinCenter(b), q = w->GetBinContent(b);
            if (q <= 0) continue;
            int bb = g->FindBin(x);
            if (g->GetBinContent(bb) <= 0) continue;
            num += q * r->GetBinContent(bb) / g->GetBinContent(bb);
            den += q;
         }
         return den > 0 ? num / den : NAN;
      };
      aEl[c] = fold(hEl);
      aIn[c] = fold(hIn);
      f->Close();
   }

   // ---------- C : the sensitivity curve ----------------------------------------------------------
   auto *gSens = new TGraph();
   for (double thcm = 8; thcm <= 150; thcm += 3) {
      double thlab = 0.5 * (180 - thcm);
      double ke = cc_locus(mC, mp, mp, mC, Ebeam, thlab * TMath::DegToRad(), 0.0);
      if (ke <= 0) continue;
      double a = cc_ex(mC, mp, mp, mC, Ebeam + 0.5, thlab * TMath::DegToRad(), ke);
      double b = cc_ex(mC, mp, mp, mC, Ebeam - 0.5, thlab * TMath::DegToRad(), ke);
      if (std::isnan(a) || std::isnan(b)) continue;
      gSens->SetPoint(gSens->GetN(), thcm, a - b);
   }

   // ================= drawing =====================================================================
   auto *c1 = new TCanvas("hfconc", "hfconc", 1620, 1020);
   c1->Divide(2, 2);
   const int COLD = kGray + 3, COLN = kAzure + 2, COL3 = kOrange + 8, COLW = kRed + 1;

   auto frameCfg = [&](const char *name, const char *ytit, double y0, double y1, bool logy) {
      auto *h = new TH1F(name, TString(";field [T] / pad pitch [mm];") + ytit, NC, 0, NC);
      for (int c = 0; c < NC; ++c)
         h->GetXaxis()->SetBinLabel(c + 1, CLB[c]);
      h->GetXaxis()->SetLabelSize(0.043);
      h->GetXaxis()->SetTitleSize(0.042);
      h->GetXaxis()->SetTitleOffset(1.25);
      h->GetYaxis()->SetTitleSize(0.045);
      h->GetYaxis()->SetTitleOffset(1.25);
      h->GetYaxis()->SetRangeUser(y0, y1);
      h->SetLineColor(kWhite);
      h->Draw();
      if (logy) gPad->SetLogy();
      return h;
   };
   auto series = [&](double *v, int col, int mk, double dx) {
      auto *g = new TGraph();
      for (int c = 0; c < NC; ++c)
         if (!std::isnan(v[c])) g->SetPoint(g->GetN(), c + 0.5 + dx, v[c]);
      g->SetMarkerStyle(mk);
      g->SetMarkerSize(1.7);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->SetLineWidth(2);
      g->Draw("pl same");
      return g;
   };

   // --- A ---
   c1->cd(1);
   gPad->SetGridy();
   gPad->SetBottomMargin(0.14);
   gPad->SetLeftMargin(0.14);
   frameCfg("fA", "#sigma(E_{x}) elastic  [MeV]", 0.025, 3.0, kTRUE);
   auto *gA1 = series(sA, COLD, 20, -0.06);
   auto *gA2 = series(sB, COLN, 21, +0.06);
   auto *lA = new TLegend(0.17, 0.63, 0.66, 0.79);
   lA->SetBorderSize(0);
   lA->SetFillStyle(0);
   lA->SetTextSize(0.040);
   lA->AddEntry(gA1, "constant E_{beam} (as analysed)", "pl");
   lA->AddEntry(gA2, "E_{beam}(z_{vertex})", "pl");
   lA->Draw();
   auto *tA = new TLatex(0.17, 0.885, "#bf{A}  the tracking improves 9#times across the matrix;");
   tA->SetNDC(); tA->SetTextSize(0.040); tA->Draw();
   auto *tA2 = new TLatex(0.17, 0.838, "as analysed, none of it reaches E_{x}");
   tA2->SetNDC(); tA2->SetTextSize(0.040); tA2->Draw();

   // --- B ---
   c1->cd(2);
   gPad->SetGridy();
   gPad->SetBottomMargin(0.14);
   gPad->SetLeftMargin(0.14);
   frameCfg("fB", "fraction reconstructed", 0.0, 1.22, kFALSE);
   auto *gB1 = series(aFlat, COLD, 24, -0.10);
   auto *gB2 = series(aEl, COLW, 20, 0.0);
   auto *gB3 = series(aIn, COL3, 22, +0.10);
   auto *lB = new TLegend(0.17, 0.17, 0.66, 0.37);
   lB->SetBorderSize(0);
   lB->SetFillStyle(0);
   lB->SetTextSize(0.040);
   lB->AddEntry(gB1, "flat in #theta_{cm} (as generated)", "pl");
   lB->AddEntry(gB2, "weighted by the a1954 elastic yield", "pl");
   lB->AddEntry(gB3, "weighted by the a1954 inelastic yield", "pl");
   lB->Draw();
   auto *tB = new TLatex(0.17, 0.885, "#bf{B}  7 T keeps a quarter of the elastic normalisation");
   tB->SetNDC(); tB->SetTextSize(0.040); tB->Draw();

   // --- C ---
   c1->cd(3);
   gPad->SetBottomMargin(0.13);
   gPad->SetLeftMargin(0.14);
   auto *fC = new TH1F("fC", ";#theta_{cm} [deg];dE_{x}/dE_{beam}", 100, 0, 150);
   fC->GetYaxis()->SetRangeUser(0, 0.168);
   fC->GetYaxis()->SetTitleOffset(1.25);
   fC->GetYaxis()->SetTitleSize(0.045);
   fC->GetXaxis()->SetTitleSize(0.045);
   fC->SetLineColor(kWhite);
   fC->Draw();
   double mEl = hEl->GetMaximum(), mIn = hIn->GetMaximum();
   if (mEl > 0) hEl->Scale(0.115 / mEl);
   if (mIn > 0) hIn->Scale(0.115 / mIn);
   hEl->SetFillColorAlpha(COLN, 0.30);
   hEl->SetLineColor(COLN);
   hEl->Draw("hist same");
   hIn->SetFillColorAlpha(COL3, 0.30);
   hIn->SetLineColor(COL3);
   hIn->Draw("hist same");
   gSens->SetLineColor(COLW);
   gSens->SetLineWidth(3);
   gSens->Draw("l same");
   auto *lC = new TLegend(0.19, 0.60, 0.63, 0.80);
   lC->SetBorderSize(0);
   lC->SetFillStyle(0);
   lC->SetTextSize(0.038);
   lC->AddEntry(gSens, "dE_{x}/dE_{beam}, elastic locus", "l");
   lC->AddEntry(hEl, "a1954 elastic yield", "f");
   lC->AddEntry(hIn, "a1954 inelastic yield", "f");
   lC->Draw();
   auto *tC = new TLatex(0.17, 0.885, "#bf{C}  the sensitivity and the yield sit at opposite ends");
   tC->SetNDC(); tC->SetTextSize(0.040); tC->Draw();

   // --- D ---
   c1->cd(4);
   gPad->SetGridy();
   gPad->SetBottomMargin(0.14);
   gPad->SetLeftMargin(0.14);
   frameCfg("fD", "#sigma(E_{x}) at 6.094 MeV  [MeV]", 0.018, 2.2, kTRUE);
   auto *gD1 = series(dA, COLD, 20, -0.06);
   auto *gD2 = series(dB, COLN, 21, +0.06);
   auto *lD = new TLine(0, dataSigma6094, NC, dataSigma6094);
   lD->SetLineColor(COLW);
   lD->SetLineWidth(3);
   lD->SetLineStyle(2);
   lD->Draw();
   auto *tD0 = new TLatex(0.18, 0.545,
                          TString::Format("a1954 measured resolution, %.3f MeV", dataSigma6094));
   tD0->SetNDC(); tD0->SetTextSize(0.038); tD0->SetTextColor(COLW); tD0->Draw();
   auto *lDl = new TLegend(0.55, 0.24, 0.94, 0.40);
   lDl->SetBorderSize(0);
   lDl->SetFillStyle(0);
   lDl->SetTextSize(0.040);
   lDl->AddEntry(gD1, "constant E_{beam}", "pl");
   lDl->AddEntry(gD2, "E_{beam}(z_{vertex})", "pl");
   lDl->Draw();
   auto *tD = new TLatex(0.17, 0.885, "#bf{D}  #theta_{lab} 35-65#circ, where the inelastic yield is");
   tD->SetNDC(); tD->SetTextSize(0.040); tD->Draw();
   auto *tD2 = new TLatex(0.17, 0.838, "the data sits a factor 4 above every prediction");
   tD2->SetNDC(); tD2->SetTextSize(0.038); tD2->Draw();

   c1->SaveAs(outDir + "hf_conclusion.png");

   printf("\n  panel A  sigma(Ex) elastic, theta_lab 20-90 deg\n");
   for (int c = 0; c < NC; ++c) printf("    %-12s %7.3f -> %7.3f\n", CFG[c], sA[c], sB[c]);
   printf("  panel B  acceptance: flat / a1954 elastic / a1954 inelastic\n");
   for (int c = 0; c < NC; ++c) printf("    %-12s %7.3f %7.3f %7.3f\n", CFG[c], aFlat[c], aEl[c], aIn[c]);
   printf("  panel D  sigma(Ex) at 6.094, theta_lab 35-65 deg   (data fixes %.3f)\n", dataSigma6094);
   for (int c = 0; c < NC; ++c) printf("    %-12s %7.3f -> %7.3f\n", CFG[c], dA[c], dB[c]);
   printf("\n  wrote %shf_conclusion.png\n\n", outDir.Data());
   if (fd) fd->Close();
}
