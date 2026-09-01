/// @file tp_kinplots_Be10.C
/// @brief The 10Be(t,p)12Be kinematics, drawn: the analytic curves, what the simulation actually
///        produces, and the two leverage factors that decide what a detector upgrade can buy --
///        with 14C(d,p)15C alongside throughout.
///
///   root -b -q 'tp_kinplots_Be10.C("/mnt/f/Be10_tp","b285_attpc","plots")'
///
/// SIX PANELS:
///   A  KE_p vs theta_lab: the four level curves over the simulated TRUTH distribution
///   B  theta_lab vs theta_cm, the four levels and the (d,p) ground state for scale
///   C  dEx/dKE vs theta_cm -- how much an error in the proton energy costs in excitation energy
///   D  dEx/dE_beam vs theta_cm -- the same for the beam. THIS is the panel that matters here:
///      it is ~3x the (d,p) value, which is why taking the beam energy at the reconstructed
///      vertex is worth a factor 2.6 in resolution and is the difference between seeing the
///      0+_2 and not seeing it.
///   E  the proton-energy separation of the 2.109/2.251 doublet, in keV, vs theta_cm -- the
///      decisive quantity: it is what the detector has to resolve, and it is only 40-80 keV
///   F  acceptance vs theta_cm, per level, from the campaign's own acceptance files
///
/// CM ANGLE CONVENTION, as everywhere in this study: theta_cm is the ANALYSIS/DWBA angle, between
/// the outgoing and the incoming light particle. Small theta_cm is the BACKWARD-lab, low-energy
/// proton, which is where a transfer angular distribution carries its yield.

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <cmath>
#include <cstdio>
#include <vector>

static const double KU = 931.49401;
// 10Be(t,p)12Be
static const double TM1 = 10.0135341 * KU, TM2 = 3.0160493 * KU, TM3 = 1.007825 * KU, TM4 = 12.0269221 * KU;
static const double TEB = 112.20;
// 14C(d,p)15C, for comparison
static const double DM1 = 14.003242 * KU, DM2 = 2.0141018 * KU, DM3 = 1.007825 * KU, DM4 = 15.0105993 * KU;
static const double DEB = 159.75;

static const int NL = 4;
static const double LEX[NL] = {0.0, 2.109, 2.251, 2.715};
static const char *LJP[NL] = {"0^{+} g.s.", "2^{+} 2.109", "0^{+}_{2} 2.251", "1^{-} 2.715"};
static const int LCOL[NL] = {kBlue, kGreen + 2, kMagenta, kOrange + 7};

static double kp_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double kp_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(th) * kp_om2(s, m1 * m1, m2 * m2) * kp_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) / (2 * m2 * m2) + s + uu - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}
/// theta_cm (ANALYSIS convention) -> (theta_lab, KE) of the proton
static bool kp_fwd(double m1, double m2, double m3, double m4, double K, double thcmA, double &thlab, double &Ke)
{
   const double thcm = TMath::Pi() - thcmA;
   double E1 = K + m1;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1;
   double rs = std::sqrt(s);
   if (rs < m3 + m4) return false;
   double pcm = kp_om2(s, m3 * m3, m4 * m4) / (2 * rs);
   double Ecm3 = std::sqrt(pcm * pcm + m3 * m3);
   double plab = std::sqrt(E1 * E1 - m1 * m1);
   double beta = plab / (E1 + m2), gam = 1.0 / std::sqrt(1 - beta * beta);
   double pz = gam * (pcm * std::cos(thcm) + beta * Ecm3);
   double pt = pcm * std::sin(thcm);
   Ke = gam * (Ecm3 + beta * pcm * std::cos(thcm)) - m3;
   thlab = std::atan2(pt, pz);
   return true;
}

void tp_kinplots_Be10(TString root = "/mnt/f/Be10_tp", TString cfg = "b285_attpc", TString outDir = "plots")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gSystem->mkdir(outDir, kTRUE);

   TCanvas *c = new TCanvas("cKin", "kinematics", 1650, 1000);
   c->Divide(3, 2);

   // ---------- A : KE vs theta_lab, curves over the simulated truth ----------
   c->cd(1);
   gPad->SetRightMargin(0.13);
   TH2D *hA = new TH2D("hA", "A: proton KE vs #theta_{lab}, curves + simulated truth;#theta_{lab} [deg];KE_{p} [MeV]",
                       180, 0, 180, 200, 0, 55);
   {
      // the TRUTH sample, all generated reactions (not just the reconstructed ones), from the g.s.
      TString sf = gSystem->GetFromPipe(
         TString::Format("ls %s/sims_b285/gs_s*_sim.root 2>/dev/null | head -1", root.Data()));
      sf = sf.Strip(TString::kBoth);
      if (!sf.IsNull()) {
         gSystem->Load("libAtSimulationData.so");
         TFile *fs = TFile::Open(sf);
         TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
         if (ts) {
            TClonesArray *mc = nullptr;
            ts->SetBranchAddress("MCTrack", &mc);
            for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
               ts->GetEntry(i);
               if (!mc) continue;
               for (int k = 0; k < mc->GetEntriesFast(); ++k) {
                  auto *t = (AtMCTrack *)mc->At(k);
                  if (!t || t->GetPdgCode() != 2212 || t->GetMotherId() != -1) continue;
                  double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000;
                  double p = std::sqrt(px * px + py * py + pz * pz);
                  if (p <= 0) break;
                  hA->Fill(std::acos(pz / p) * TMath::RadToDeg(), std::sqrt(p * p + TM3 * TM3) - TM3);
                  break;
               }
            }
            fs->Close();
         }
      }
   }
   // The truth map is the point of this panel, so give it a log z -- on a linear scale the
   // sparsely populated forward region vanished entirely under the palette maximum.
   gPad->SetLogz();
   hA->Draw("colz");
   TGraph *gKE[NL];
   for (int l = 0; l < NL; ++l) {
      gKE[l] = new TGraph();
      int n = 0;
      for (double a = 1; a <= 179; a += 0.5) {
         double th, ke;
         if (!kp_fwd(TM1, TM2, TM3, TM4 + LEX[l], TEB, a * TMath::DegToRad(), th, ke)) continue;
         if (ke <= 0) continue;
         gKE[l]->SetPoint(n++, th * TMath::RadToDeg(), ke);
      }
      gKE[l]->SetLineColor(LCOL[l]);
      gKE[l]->SetLineWidth(2);
      gKE[l]->Draw("L SAME");
   }
   TLegend *lA = new TLegend(0.45, 0.62, 0.86, 0.88);
   lA->SetFillColorAlpha(kWhite, 0.55);
   lA->SetBorderSize(0);
   for (int l = 0; l < NL; ++l) lA->AddEntry(gKE[l], LJP[l], "l");
   lA->Draw();

   // ---------- B : theta_lab vs theta_cm ----------
   c->cd(2);
   TH2D *fB = new TH2D("fB", "B: #theta_{lab} vs #theta_{cm} (DWBA angle);#theta_{cm} [deg];#theta_{lab} [deg]", 10, 0,
                       180, 10, 0, 180);
   fB->Draw();
   TGraph *gTL[NL];
   for (int l = 0; l < NL; ++l) {
      gTL[l] = new TGraph();
      int n = 0;
      for (double a = 1; a <= 179; a += 0.5) {
         double th, ke;
         if (!kp_fwd(TM1, TM2, TM3, TM4 + LEX[l], TEB, a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
         gTL[l]->SetPoint(n++, a, th * TMath::RadToDeg());
      }
      gTL[l]->SetLineColor(LCOL[l]);
      gTL[l]->SetLineWidth(2);
      gTL[l]->Draw("L SAME");
   }
   TGraph *gDP = new TGraph();
   {
      int n = 0;
      for (double a = 1; a <= 179; a += 0.5) {
         double th, ke;
         if (!kp_fwd(DM1, DM2, DM3, DM4, DEB, a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
         gDP->SetPoint(n++, a, th * TMath::RadToDeg());
      }
      gDP->SetLineColor(kGray + 2);
      gDP->SetLineStyle(2);
      gDP->SetLineWidth(2);
      gDP->Draw("L SAME");
   }
   TLegend *lB = new TLegend(0.45, 0.60, 0.90, 0.88);
   lB->SetFillStyle(0);
   lB->SetBorderSize(0);
   for (int l = 0; l < NL; ++l) lB->AddEntry(gTL[l], LJP[l], "l");
   lB->AddEntry(gDP, "^{14}C(d,p) g.s.", "l");
   lB->Draw();

   // ---------- C, D : the two leverage factors ----------
   auto lever = [&](int pad, bool beam, const char *title, double ymax) {
      c->cd(pad);
      TH2D *fr = new TH2D(Form("fr%d", pad), Form("%s;#theta_{cm} [deg];%s", title, beam ? "|dE_{x}/dE_{beam}|" : "|dE_{x}/dKE_{p}|"),
                          10, 0, 180, 10, 0, ymax);
      fr->Draw();
      TGraph *g1 = new TGraph(), *g2 = new TGraph();
      int n1 = 0, n2 = 0;
      for (double a = 2; a <= 178; a += 1.0) {
         double th, ke, h = 0.01;
         if (kp_fwd(TM1, TM2, TM3, TM4, TEB, a * TMath::DegToRad(), th, ke) && ke > 0.2) {
            double d = beam ? (kp_ex(TM1, TM2, TM3, TM4, TEB + h, th, ke) - kp_ex(TM1, TM2, TM3, TM4, TEB - h, th, ke))
                            : (kp_ex(TM1, TM2, TM3, TM4, TEB, th, ke + h) - kp_ex(TM1, TM2, TM3, TM4, TEB, th, ke - h));
            // Defensive only: the kinematics are finite over the whole range (checked -- dEx/dEbeam
            // runs smoothly to 0.341 at theta_cm 178). This is NOT what made the curve look
            // truncated; that was the legend drawn on top of it. Kept because a silently dropped
            // point is indistinguishable from a physical edge.
            if (std::isfinite(d)) g1->SetPoint(n1++, a, std::fabs(d / (2 * h)));
         }
         if (kp_fwd(DM1, DM2, DM3, DM4, DEB, a * TMath::DegToRad(), th, ke) && ke > 0.2) {
            double d = beam ? (kp_ex(DM1, DM2, DM3, DM4, DEB + h, th, ke) - kp_ex(DM1, DM2, DM3, DM4, DEB - h, th, ke))
                            : (kp_ex(DM1, DM2, DM3, DM4, DEB, th, ke + h) - kp_ex(DM1, DM2, DM3, DM4, DEB, th, ke - h));
            if (std::isfinite(d)) g2->SetPoint(n2++, a, std::fabs(d / (2 * h)));
         }
      }
      g1->SetLineColor(kRed);
      g1->SetLineWidth(3);
      g1->Draw("L SAME");
      g2->SetLineColor(kGray + 2);
      g2->SetLineStyle(2);
      g2->SetLineWidth(3);
      g2->Draw("L SAME");
      // LEGEND PLACEMENT IS NOT COSMETIC HERE. With the legend at upper right, panel D's (t,p)
      // curve -- which rises to 0.341 at theta_cm 178 -- ran underneath an OPAQUE legend box and
      // appeared to stop dead at 137 deg. It read exactly like a channel closing, and the first
      // thing I did was add a NaN guard for a NaN that did not exist. The values were always
      // fine; the box was on top of them. Legends are transparent from here on, so an overlap
      // shows as an overlap instead of as missing data, and D's sits where its curve is not.
      TLegend *lg = beam ? new TLegend(0.16, 0.68, 0.60, 0.86) : new TLegend(0.45, 0.72, 0.90, 0.88);
      lg->SetFillStyle(0);
      lg->SetBorderSize(0);
      lg->AddEntry(g1, "^{10}Be(t,p)^{12}Be", "l");
      lg->AddEntry(g2, "^{14}C(d,p)^{15}C", "l");
      lg->Draw();
      // the window the analysis uses
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.035);
      tx.DrawLatex(0.16, 0.20, "transfer peak: #theta_{cm} 12-45#circ");
   };
   lever(3, false, "C: cost of a proton-energy error", 3.5);
   lever(4, true, "D: cost of a beam-energy error", 0.40);

   // ---------- E : the doublet separation in proton energy ----------
   c->cd(5);
   TH2D *fE = new TH2D("fE", "E: what must be resolved: KE_{p}(2.109) - KE_{p}(2.251);#theta_{cm} [deg];#DeltaKE_{p} [keV]",
                       10, 0, 180, 10, 0, 220);
   fE->Draw();
   TGraph *gSep = new TGraph();
   {
      int n = 0;
      for (double a = 2; a <= 178; a += 1.0) {
         double th1, ke1, th2, ke2;
         if (!kp_fwd(TM1, TM2, TM3, TM4 + LEX[1], TEB, a * TMath::DegToRad(), th1, ke1)) continue;
         if (!kp_fwd(TM1, TM2, TM3, TM4 + LEX[2], TEB, a * TMath::DegToRad(), th2, ke2)) continue;
         if (ke1 <= 0 || ke2 <= 0) continue;
         gSep->SetPoint(n++, a, 1000 * (ke1 - ke2));
      }
      gSep->SetLineColor(kMagenta + 2);
      gSep->SetLineWidth(3);
      gSep->Draw("L SAME");
   }
   {
      TLatex tx;
      tx.SetNDC();
      tx.SetTextSize(0.033);
      tx.DrawLatex(0.15, 0.62, "142 keV apart in excitation energy,");
      tx.DrawLatex(0.15, 0.56, "but only 40-80 keV apart in PROTON");
      tx.DrawLatex(0.15, 0.50, "ENERGY over the useful range.");
   }

   // ---------- F : acceptance vs theta_cm, per level ----------
   c->cd(6);
   TH2D *fF = new TH2D("fF", Form("F: acceptance vs #theta_{cm}, %s;#theta_{cm} [deg];acceptance", cfg.Data()), 10, 0,
                       180, 10, 0, 1.05);
   fF->Draw();
   const char *LVN[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
   TLegend *lF = new TLegend(0.45, 0.16, 0.90, 0.40);
   lF->SetFillStyle(0);
   lF->SetBorderSize(0);
   for (int l = 0; l < NL; ++l) {
      TString fa = gSystem->GetFromPipe(TString::Format("ls %s/%s/acceptance_%s_s*_%s.root 2>/dev/null | head -1",
                                                        root.Data(), cfg.Data(), LVN[l], cfg.Data()));
      fa = fa.Strip(TString::kBoth);
      if (fa.IsNull()) continue;
      TFile *f = TFile::Open(fa);
      if (!f || f->IsZombie()) continue;
      TH1D *hg = nullptr, *hr = nullptr;
      TIter nx(f->GetListOfKeys());
      while (auto *k = (TKey *)nx()) {
         TString nm = k->GetName();
         if (nm.BeginsWith("hGen_")) hg = (TH1D *)f->Get(nm);
         if (nm.BeginsWith("hRec_")) hr = (TH1D *)f->Get(nm);
      }
      if (!hg || !hr) { f->Close(); continue; }
      auto *g = new TGraph();
      int n = 0;
      for (int b = 1; b <= hg->GetNbinsX(); ++b)
         if (hg->GetBinContent(b) > 20)
            g->SetPoint(n++, hg->GetBinCenter(b), hr->GetBinContent(b) / hg->GetBinContent(b));
      g->SetLineColor(LCOL[l]);
      g->SetLineWidth(2);
      g->SetMarkerColor(LCOL[l]);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.6);
      g->DrawClone("LP SAME");
      lF->AddEntry(g, LJP[l], "lp");
      f->Close();
   }
   lF->Draw();

   TString png = outDir + "/kinematics_Be10tp.png";
   c->SaveAs(png);
   printf("\nwrote %s\n", png.Data());

   // the numbers behind panels C, D and E, over the window the analysis uses
   printf("\n  %-9s %10s %10s | %12s %12s | %12s\n", "theta_cm", "theta_lab", "KE_p", "dEx/dKE", "dEx/dEbeam",
          "dKE(doublet)");
   for (double a : {5., 10., 15., 20., 25., 30., 35., 40., 45., 60., 90.}) {
      double th, ke, h = 0.01, th2, ke2;
      if (!kp_fwd(TM1, TM2, TM3, TM4, TEB, a * TMath::DegToRad(), th, ke)) continue;
      double dK = (kp_ex(TM1, TM2, TM3, TM4, TEB, th, ke + h) - kp_ex(TM1, TM2, TM3, TM4, TEB, th, ke - h)) / (2 * h);
      double dE = (kp_ex(TM1, TM2, TM3, TM4, TEB + h, th, ke) - kp_ex(TM1, TM2, TM3, TM4, TEB - h, th, ke)) / (2 * h);
      double s1, s2;
      kp_fwd(TM1, TM2, TM3, TM4 + LEX[1], TEB, a * TMath::DegToRad(), th2, s1);
      kp_fwd(TM1, TM2, TM3, TM4 + LEX[2], TEB, a * TMath::DegToRad(), th2, s2);
      printf("  %7.0f   %10.1f %10.2f | %12.3f %12.4f | %9.0f keV\n", a, th * TMath::RadToDeg(), ke, std::fabs(dK),
             std::fabs(dE), 1000 * (s1 - s2));
   }
   printf("\nkinematics plots done\n\n");
}
