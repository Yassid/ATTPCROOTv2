/// @file dp_conclusion_C14.C
/// @brief The 14C(d,p)15C campaign in one figure, laid against the (p,p') result.
///
///   root -b -q 'dp_conclusion_C14.C()'
///
/// The question this answers: (p,p') gained almost nothing from a higher field because its
/// accepted protons are slow, forward, and already well measured. A (d,p) proton at the angles
/// where transfer peaks goes BACKWARD at 125 deg with 3 MeV. Does the matrix pay off there?
///
/// It does, and not the way the kinematic estimate predicted. Extrapolating the (p,p') resolution
/// to a 3 MeV backward proton gave sigma(KE) ~ 0.02 MeV; the real value is 0.43, twenty times
/// worse, because a backward low-energy track at 2.85 T is nearly straight over its short length
/// and its curvature is barely measured. Raising the field is what fixes that -- the pad pitch
/// cannot, because at 2.85 T transverse diffusion (3.7 mm at full drift) is already larger than a
/// 2 mm pad.

#include <algorithm>
#include <vector>

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const char *CLB[NC] = {"2.85 / 8#times12", "2.85 / 2mm", "4 / 8#times12",
                              "4 / 2mm",          "7 / 8#times12", "7 / 2mm"};
static const int NS = 4;
static const double cmLo[NS] = {8, 30, 60, 100};
static const double cmHi[NS] = {30, 60, 100, 180};
static const char *SLB[NS] = {"#theta_{cm} 8-30#circ  (p backward, 3 MeV)", "#theta_{cm} 30-60#circ  (85#circ, 10 MeV)",
                              "#theta_{cm} 60-100#circ  (60#circ, 22 MeV)", "#theta_{cm} 100-180#circ  (34#circ, 40 MeV)"};

static const double U = 931.49401;
static const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;

static double c_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double c_ex(double m4, double K, double th, double Ke)
{
   double Et1 = K + M1, Et3 = Ke + M3;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * Et1;
   double uu = M2 * M2 + M3 * M3 - 2 * M2 * Et3;
   double a = (std::cos(th) * c_om2(s, M1 * M1, M2 * M2) * c_om2(uu, M2 * M2, M3 * M3) -
               (s - M1 * M1 - M2 * M2) * (M2 * M2 + M3 * M3 - uu)) /
                 (2 * M2 * M2) +
              s + uu - M2 * M2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
static double c_q(std::vector<double> v, double p)
{
   if (v.size() < 20) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_conclusion_C14(TString root = "/mnt/f/a1954_C14dp_hf", Double_t Ebeam = 159.75, TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   // [slice][cfg] for each level
   double sConst[2][NS][NC], sZ[2][NS][NC], sKE[NS][NC];
   for (int l = 0; l < 2; ++l)
      for (int s = 0; s < NS; ++s)
         for (int c = 0; c < NC; ++c) { sConst[l][s][c] = sZ[l][s][c] = NAN; if (l == 0) sKE[s][c] = NAN; }

   const char *LVL[2] = {"gs", "ex0740"};
   const double LEX[2] = {0.0, 0.740};
   for (int l = 0; l < 2; ++l)
      for (int c = 0; c < NC; ++c) {
         TString f = gSystem->GetFromPipe(TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1",
                                                          root.Data(), CFG[c], LVL[l], CFG[c]));
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
         const double mres = M4 + LEX[l];
         std::vector<double> eb, zz;
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double lo = 100., hi = 200.;
            double flo = c_ex(mres, lo, thTrue * TMath::DegToRad(), keTrue);
            double fhi = c_ex(mres, hi, thTrue * TMath::DegToRad(), keTrue);
            if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
            for (int it = 0; it < 60; ++it) {
               double m = 0.5 * (lo + hi), fm = c_ex(mres, m, thTrue * TMath::DegToRad(), keTrue);
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
            std::vector<double> dC, dZ, dK;
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               if (cmTrue < cmLo[s] || cmTrue >= cmHi[s]) continue;
               dC.push_back(exReco - LEX[l]);
               dK.push_back(keReco - keTrue);
               double a = c_ex(M4, fEb.Eval(zReco), thReco * TMath::DegToRad(), keReco);
               if (!std::isnan(a)) dZ.push_back(a - LEX[l]);
            }
            sConst[l][s][c] = (c_q(dC, .75) - c_q(dC, .25)) / 1.349;
            sZ[l][s][c] = (c_q(dZ, .75) - c_q(dZ, .25)) / 1.349;
            if (l == 0) sKE[s][c] = (c_q(dK, .75) - c_q(dK, .25)) / 1.349;
         }
         fr->Close();
      }

   auto *cv = new TCanvas("dpc", "dpc", 1620, 1020);
   cv->Divide(2, 2);
   const int col[NS] = {kRed + 1, kOrange + 8, kGreen + 3, kAzure + 2};
   auto frame = [&](int pad, const char *yt, double y0, double y1, bool logy) {
      cv->cd(pad);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.15);
      gPad->SetGridy();
      if (logy) gPad->SetLogy();
      auto *h = new TH1F(Form("fr%d", pad), TString(";field [T] / pad pitch [mm];") + yt, NC, 0, NC);
      for (int c = 0; c < NC; ++c) h->GetXaxis()->SetBinLabel(c + 1, CLB[c]);
      h->GetXaxis()->SetLabelSize(0.043);
      h->GetXaxis()->SetTitleSize(0.042);
      h->GetXaxis()->SetTitleOffset(1.3);
      h->GetYaxis()->SetTitleSize(0.045);
      h->GetYaxis()->SetTitleOffset(1.25);
      h->GetYaxis()->SetRangeUser(y0, y1);
      h->SetLineColor(kWhite);
      h->Draw();
      return h;
   };
   auto line = [&](double *v, int s, int sty) {
      auto *g = new TGraph();
      for (int c = 0; c < NC; ++c)
         if (!std::isnan(v[c])) g->SetPoint(g->GetN(), c + 0.5, v[c]);
      g->SetLineColor(col[s]);
      g->SetMarkerColor(col[s]);
      g->SetMarkerStyle(sty == 1 ? 20 : 24);
      g->SetMarkerSize(1.5);
      g->SetLineWidth(3);
      g->SetLineStyle(sty);
      g->Draw("pl same");
      return g;
   };

   // A : sigma(Ex), as analysed
   frame(1, "#sigma(E_{x}) g.s. [MeV]", 0.03, 4.0, kTRUE);
   auto *lgA = new TLegend(0.17, 0.14, 0.72, 0.40);
   lgA->SetBorderSize(0); lgA->SetFillStyle(0); lgA->SetTextSize(0.036);
   for (int s = 0; s < NS; ++s) lgA->AddEntry(line(sConst[0][s], s, 1), SLB[s], "pl");
   lgA->Draw();
   auto *tA = new TLatex(0.16, 0.90, "#bf{A}  constant E_{beam}, as the analysis does it");
   tA->SetNDC(); tA->SetTextSize(0.040); tA->Draw();

   // B : sigma(Ex), corrected
   frame(2, "#sigma(E_{x}) g.s. [MeV]", 0.03, 4.0, kTRUE);
   for (int s = 0; s < NS; ++s) line(sZ[0][s], s, 1);
   auto *tB = new TLatex(0.16, 0.90, "#bf{B}  E_{beam} at the reconstructed vertex");
   tB->SetNDC(); tB->SetTextSize(0.040); tB->Draw();
   auto *tB2 = new TLatex(0.17, 0.17, "same colours as A");
   tB2->SetNDC(); tB2->SetTextSize(0.034); tB2->SetTextColor(kGray + 3); tB2->Draw();

   // C : the mechanism -- sigma(KE)
   frame(3, "#sigma(KE_{reco} - KE_{true}) [MeV]", 0.02, 6.0, kTRUE);
   for (int s = 0; s < NS; ++s) line(sKE[s], s, 1);
   auto *tC = new TLatex(0.16, 0.90, "#bf{C}  why: the ejectile energy measurement");
   tC->SetNDC(); tC->SetTextSize(0.040); tC->Draw();

   // D : separation of the 0.740 doublet, corrected
   frame(4, "#DeltaE / (#sigma_{gs} + #sigma_{0.740})", 0.1, 30, kTRUE);
   double sep[NS][NC];
   for (int s = 0; s < NS; ++s)
      for (int c = 0; c < NC; ++c) sep[s][c] = 0.740 / (sZ[0][s][c] + sZ[1][s][c]);
   for (int s = 0; s < NS; ++s) line(sep[s], s, 1);
   auto *l2 = new TLine(0, 2, NC, 2);
   l2->SetLineStyle(2); l2->SetLineColor(kGray + 2); l2->Draw();
   auto *tD = new TLatex(0.16, 0.90, "#bf{D}  the 15C doublet, 740 keV apart");
   tD->SetNDC(); tD->SetTextSize(0.040); tD->Draw();
   auto *tD2 = new TLatex(0.18, 0.615, "resolvable");
   tD2->SetNDC(); tD2->SetTextSize(0.034); tD2->SetTextColor(kGray + 3); tD2->Draw();

   cv->SaveAs(outDir + "dp_conclusion.png");
   printf("\n  wrote %sdp_conclusion.png\n\n", outDir.Data());
}
