/// @file tp_fb_Be10.C
/// @brief 10Be(t,p)12Be kinematics and excitation energy, split TRANSFER PEAK / FORWARD.
///
///   root -b -q 'tp_fb_Be10.C("/mnt/f/Be10_tp")'
///
/// The (d,p) campaign split forward/backward in the LAB, because there the two regimes divided at
/// theta_lab 90. Here the natural divider is theta_cm 45, which is the same physical split stated
/// in the variable the physics is quoted in: below it the proton is slow (3-10 MeV) and BACKWARD
/// (theta_lab 99-168), spiralling long enough that its momentum is over-determined; above it the
/// proton is fast (15-49 MeV) and forward, on a short arc whose curvature is barely measured.
///
/// Measured, ground state at 2.85 T: sigma(KE)/KE is 0.7 % below the divider and 3-6 % above it,
/// and sigma(Ex) is 0.08-0.10 MeV against 0.6-1.4. They are two different measurements sharing a
/// detector, and any number averaged over both describes neither.
///
/// Writes three figures:
///   tp_fb_kinematics.png  the (theta_lab, KE) plane per configuration, with the two-body loci of
///                         all four levels and the theta_cm = 45 divider drawn in
///   tp_fb_excitation.png  E_x per configuration, the two regions overlaid, each normalised to
///                         itself so the WIDTHS compare rather than the yields
///   tp_fb_summary.png     sigma(E_x) and acceptance, peak vs forward, across the matrix
///
/// Missing samples are skipped and reported, so this runs while a campaign is still going.

#include "tp_common_Be10.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <cstdio>

void tp_fb_Be10(TString root = "/mnt/f/Be10_tp", TString outDir = "plots", Bool_t vertexEbeam = kTRUE)
{
   using namespace tpc;
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gSystem->mkdir(outDir, kTRUE);
   const double DIV = PEAK_HI; // theta_cm divider

   // ---------------- 1. kinematics plane, per configuration ----------------
   auto *c1 = new TCanvas("cFbK", "kinematics", 1650, 950);
   c1->Divide(3, 2);
   for (int c = 0; c < NC; ++c) {
      c1->cd(c + 1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.13);
      auto *h = new TH2D(Form("hk%d", c), Form("%s;#theta_{lab} [deg];KE_{p} [MeV]", CFGL[c]), 180, 0, 180, 160, 0, 55);
      h->SetDirectory(nullptr);
      bool any = false;
      for (int l = 0; l < NL; ++l) {
         TString fn = find(root, CFG[c], LVN[l], "exres");
         if (fn.IsNull()) continue;
         TFile *f = TFile::Open(fn);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (!t) { if (f) f->Close(); continue; }
         double thR, keR;
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keReco", &keR);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); h->Fill(thR, keR); }
         f->Close();
         any = true;
      }
      if (!any) { printf("  kinematics: %s missing\n", CFG[c]); continue; }
      h->Draw("colz");
      for (int l = 0; l < NL; ++l) {
         auto *g = new TGraph();
         int n = 0;
         for (double a = 0.5; a <= 179.5; a += 0.25) {
            double th, ke;
            if (!fwd(M4 + LEX[l], a * TMath::DegToRad(), th, ke) || ke <= 0) continue;
            g->SetPoint(n++, th * TMath::RadToDeg(), ke);
         }
         g->SetLineColor(LCOL[l]);
         g->SetLineWidth(2);
         g->Draw("L SAME");
      }
      double thd, ked;
      if (fwd(M4, DIV * TMath::DegToRad(), thd, ked)) {
         auto *ln = new TLine(thd * TMath::RadToDeg(), 0, thd * TMath::RadToDeg(), 55);
         ln->SetLineColor(kGreen + 2);
         ln->SetLineStyle(7);
         ln->SetLineWidth(2);
         ln->Draw();
         TLatex t;
         t.SetTextSize(0.038);
         t.SetTextColor(kGreen + 2);
         t.DrawLatex(thd * TMath::RadToDeg() - 42, 48, Form("#theta_{cm}=%.0f#circ", DIV));
      }
   }
   c1->SaveAs(outDir + "/tp_fb_kinematics.png");

   // ---------------- 2. Ex, the two regions overlaid, self-normalised ----------------
   auto *c2 = new TCanvas("cFbE", "excitation", 1650, 950);
   c2->Divide(3, 2);
   struct Res { double sPeak{NAN}, sFwd{NAN}, aPeak{NAN}, aFwd{NAN}; long nPeak{0}, nFwd{0}; };
   Res R[NC];
   for (int c = 0; c < NC; ++c) {
      c2->cd(c + 1);
      auto *hp = new TH1D(Form("hp%d", c), Form("%s;E_{x}(^{12}Be) [MeV];normalised", CFGL[c]), 220, -2.5, 5.0);
      auto *hf = new TH1D(Form("hf%d", c), "", 220, -2.5, 5.0);
      hp->SetDirectory(nullptr);
      hf->SetDirectory(nullptr);
      std::vector<double> dp, df;
      bool any = false;
      for (int l = 0; l < NL; ++l) {
         TString fn = find(root, CFG[c], LVN[l], "exres");
         if (fn.IsNull()) continue;
         TFile *f = TFile::Open(fn);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (!t) { if (f) f->Close(); continue; }
         TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
         bool okProf = vertexEbeam ? ebeamProfile(t, LEX[l], fEb) : false;
         double thR, keR, cmT, exR, zR;
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keReco", &keR);
         t->SetBranchAddress("cmTrue", &cmT);
         t->SetBranchAddress("exReco", &exR);
         t->SetBranchAddress("zReco", &zR);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            double e = (vertexEbeam && okProf) ? ex(M4, fEb.Eval(zR), thR * TMath::DegToRad(), keR) : exR;
            if (std::isnan(e)) continue;
            if (cmT < DIV) { hp->Fill(e); dp.push_back(e - LEX[l]); }
            else           { hf->Fill(e); df.push_back(e - LEX[l]); }
         }
         f->Close();
         any = true;
      }
      if (!any) { printf("  excitation: %s missing\n", CFG[c]); continue; }
      R[c].nPeak = dp.size();
      R[c].nFwd = df.size();
      R[c].sPeak = (quant(dp, .75) - quant(dp, .25)) / 1.349;
      R[c].sFwd = (quant(df, .75) - quant(df, .25)) / 1.349;
      if (hp->Integral() > 0) hp->Scale(1. / hp->Integral());
      if (hf->Integral() > 0) hf->Scale(1. / hf->Integral());
      hp->SetLineColor(kBlue + 1);
      hp->SetLineWidth(2);
      hf->SetLineColor(kRed + 1);
      hf->SetLineWidth(2);
      hp->SetMaximum(1.35 * std::max(hp->GetMaximum(), hf->GetMaximum()));
      hp->Draw("HIST");
      hf->Draw("HIST SAME");
      auto *lg = new TLegend(0.55, 0.68, 0.97, 0.88);
      lg->SetFillStyle(0);
      lg->SetBorderSize(0);
      lg->SetTextSize(0.036);
      lg->AddEntry(hp, Form("#theta_{cm}<%.0f  #sigma=%.3f", DIV, R[c].sPeak), "l");
      lg->AddEntry(hf, Form("#theta_{cm}>%.0f  #sigma=%.3f", DIV, R[c].sFwd), "l");
      lg->Draw();
   }
   c2->SaveAs(outDir + "/tp_fb_excitation.png");

   // ---------------- 3. summary across the matrix ----------------
   for (int c = 0; c < NC; ++c) {
      TH1D *hg = nullptr, *hr = nullptr;
      TFile *f = nullptr;
      if (!acceptanceHists(root, CFG[c], "gs", hg, hr, f)) continue;
      double gp = hg->Integral(hg->FindBin(PEAK_LO + .1), hg->FindBin(DIV - .1));
      double rp = hr->Integral(hr->FindBin(PEAK_LO + .1), hr->FindBin(DIV - .1));
      double gf = hg->Integral(hg->FindBin(DIV + .1), hg->GetNbinsX());
      double rf = hr->Integral(hr->FindBin(DIV + .1), hr->GetNbinsX());
      R[c].aPeak = gp > 0 ? rp / gp : NAN;
      R[c].aFwd = gf > 0 ? rf / gf : NAN;
      f->Close();
   }
   auto *c3 = new TCanvas("cFbS", "summary", 1400, 620);
   c3->Divide(2, 1);
   auto bar = [&](int pad, bool res, const char *ttl, const char *yt, double ymax) {
      c3->cd(pad);
      auto *fr = new TH2D(Form("frs%d", pad), Form("%s;configuration;%s", ttl, yt), NC, 0, NC, 10, 0, ymax);
      for (int c = 0; c < NC; ++c) fr->GetXaxis()->SetBinLabel(c + 1, CFGL[c]);
      fr->GetXaxis()->LabelsOption("v");
      fr->Draw();
      auto *gp = new TGraph(), *gf = new TGraph();
      int np = 0, nf = 0;
      for (int c = 0; c < NC; ++c) {
         double vp = res ? R[c].sPeak : R[c].aPeak, vf = res ? R[c].sFwd : R[c].aFwd;
         if (!std::isnan(vp)) gp->SetPoint(np++, c + 0.35, vp);
         if (!std::isnan(vf)) gf->SetPoint(nf++, c + 0.65, vf);
      }
      gp->SetMarkerStyle(20); gp->SetMarkerSize(1.6); gp->SetMarkerColor(kBlue + 1);
      gf->SetMarkerStyle(21); gf->SetMarkerSize(1.6); gf->SetMarkerColor(kRed + 1);
      gp->Draw("P SAME"); gf->Draw("P SAME");
      auto *lg = new TLegend(0.45, 0.76, 0.97, 0.90);
      lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.034);
      lg->AddEntry(gp, Form("#theta_{cm} < %.0f#circ (transfer peak)", DIV), "p");
      lg->AddEntry(gf, Form("#theta_{cm} > %.0f#circ (forward)", DIV), "p");
      lg->Draw();
   };
   bar(1, true, "excitation-energy resolution", "#sigma(E_{x}) [MeV]", 1.6);
   bar(2, false, "acceptance", "accepted / generated", 1.05);
   c3->SaveAs(outDir + "/tp_fb_summary.png");

   printf("\n=== 10Be(t,p)12Be, split at theta_cm = %.0f  (%s beam energy) ===\n", DIV,
          vertexEbeam ? "VERTEX" : "constant");
   printf("  %-16s %10s %10s %10s %10s %8s %8s\n", "config", "sig peak", "sig fwd", "acc peak", "acc fwd", "n peak",
          "n fwd");
   for (int c = 0; c < NC; ++c)
      printf("  %-16s %10.3f %10.3f %10.3f %10.3f %8ld %8ld\n", CFGL[c], R[c].sPeak, R[c].sFwd, R[c].aPeak, R[c].aFwd,
             R[c].nPeak, R[c].nFwd);
   printf("\nwrote %s/tp_fb_{kinematics,excitation,summary}.png\nfb done\n\n", outDir.Data());
}
