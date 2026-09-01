/// @file tp_acceptance_Be10.C
/// @brief Acceptance of 10Be(t,p)12Be vs theta_cm, one curve per configuration.
///
///   root -b -q 'tp_acceptance_Be10.C("/mnt/f/Be10_tp")'
///
/// Same construction as the (d,p) campaign's dp_acceptance_C14.C, and for the same reason: the
/// OVERALL acceptance is 0.63-0.68 in every cell of the matrix, which makes the field look nearly
/// free. It is not. Plotting the ratio to the present detector (right panel) is what makes the
/// cost visible -- the left panel alone would let 7 T pass as harmless, when in fact it keeps only
/// ~60 % of what 2.85 T keeps at the transfer peak while GAINING at forward theta_cm, where the
/// fast protons curl up enough to be measured. Those two effects nearly cancel in the overall
/// number, which is exactly why the overall number must not be quoted on its own.

#include "tp_common_Be10.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <cstdio>

void tp_acceptance_Be10(TString root = "/mnt/f/Be10_tp", TString lvl = "gs", TString outDir = "plots")
{
   using namespace tpc;
   gStyle->SetOptStat(0);
   gSystem->mkdir(outDir, kTRUE);

   TGraph *g[NC] = {}, *r[NC] = {};
   TH1D *hgRef = nullptr, *hrRef = nullptr;
   TFile *fRef = nullptr;
   bool haveRef = acceptanceHists(root, CFG[0], lvl, hgRef, hrRef, fRef);

   printf("\n=== 10Be(t,p)12Be acceptance vs theta_cm, level %s ===\n", lvl.Data());
   printf("  %-16s %8s %8s %8s %8s %8s | %9s\n", "config", "2-20", "20-45", "45-90", "90-135", "135-180", "overall");
   for (int c = 0; c < NC; ++c) {
      TH1D *hg = nullptr, *hr = nullptr;
      TFile *f = nullptr;
      if (!acceptanceHists(root, CFG[c], lvl, hg, hr, f)) {
         printf("  %-16s (missing)\n", CFGL[c]);
         continue;
      }
      g[c] = new TGraph();
      r[c] = new TGraph();
      int n = 0, m = 0;
      for (int b = 1; b <= hg->GetNbinsX(); ++b) {
         if (hg->GetBinContent(b) < 20) continue;
         double a = hr->GetBinContent(b) / hg->GetBinContent(b);
         g[c]->SetPoint(n++, hg->GetBinCenter(b), a);
         if (haveRef && hgRef->GetBinContent(b) >= 20 && hrRef->GetBinContent(b) > 0)
            r[c]->SetPoint(m++, hg->GetBinCenter(b), a / (hrRef->GetBinContent(b) / hgRef->GetBinContent(b)));
      }
      double lo[5] = {2, 20, 45, 90, 135}, hi[5] = {20, 45, 90, 135, 180};
      printf("  %-16s", CFGL[c]);
      for (int s = 0; s < 5; ++s) {
         double gg = hg->Integral(hg->FindBin(lo[s] + .1), hg->FindBin(hi[s] - .1));
         double rr = hr->Integral(hr->FindBin(lo[s] + .1), hr->FindBin(hi[s] - .1));
         printf(" %8.3f", gg > 0 ? rr / gg : 0.);
      }
      printf(" | %9.3f\n", hg->Integral() > 0 ? hr->Integral() / hg->Integral() : 0.);
      for (auto *gr : {g[c], r[c]}) {
         gr->SetLineColor(CCOL[c]);
         gr->SetLineStyle(CSTY[c]);
         gr->SetLineWidth(2);
      }
   }

   auto *c1 = new TCanvas("cAcc", "acceptance", 1500, 620);
   c1->Divide(2, 1);
   c1->cd(1);
   auto *fr1 = new TH2D("fr1", Form("acceptance, %s;#theta_{cm} [deg];accepted / generated", lvl.Data()), 10, 0, 180,
                        10, 0, 1.05);
   fr1->Draw();
   auto *l1 = new TLegend(0.14, 0.14, 0.52, 0.42);
   l1->SetFillStyle(0);
   l1->SetBorderSize(0);
   l1->SetTextSize(0.030);
   for (int c = 0; c < NC; ++c)
      if (g[c]) { g[c]->Draw("L SAME"); l1->AddEntry(g[c], CFGL[c], "l"); }
   l1->Draw();
   // mark the window the analysis uses
   for (double x : {PEAK_LO, PEAK_HI}) {
      auto *ln = new TLine(x, 0, x, 1.05);
      ln->SetLineColor(kGreen + 2);
      ln->SetLineStyle(7);
      ln->Draw();
   }
   {
      TLatex t; t.SetNDC(); t.SetTextSize(0.030); t.SetTextColor(kGreen + 2);
      t.DrawLatex(0.16, 0.90, Form("#theta_{cm} %.0f-%.0f#circ = the transfer peak", PEAK_LO, PEAK_HI));
   }

   c1->cd(2);
   auto *fr2 = new TH2D("fr2", "ratio to the present detector (2.85 T, AT-TPC);#theta_{cm} [deg];acceptance ratio", 10,
                        0, 180, 10, 0, 2.0);
   fr2->Draw();
   auto *one = new TLine(0, 1, 180, 1);
   one->SetLineColor(kBlack);
   one->SetLineStyle(3);
   one->Draw();
   for (int c = 1; c < NC; ++c)
      if (r[c]) r[c]->Draw("L SAME");
   auto *l2 = new TLegend(0.14, 0.68, 0.52, 0.90);
   l2->SetFillStyle(0);
   l2->SetBorderSize(0);
   l2->SetTextSize(0.030);
   for (int c = 1; c < NC; ++c)
      if (r[c]) l2->AddEntry(r[c], CFGL[c], "l");
   l2->Draw();
   for (double x : {PEAK_LO, PEAK_HI}) {
      auto *ln = new TLine(x, 0, x, 2.0);
      ln->SetLineColor(kGreen + 2);
      ln->SetLineStyle(7);
      ln->Draw();
   }
   TString png = outDir + "/tp_acceptance_" + lvl + ".png";
   c1->SaveAs(png);
   printf("\nwrote %s\nacceptance done\n\n", png.Data());
}
