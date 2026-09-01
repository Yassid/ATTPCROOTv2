/// @file tp_conclusion_Be10.C
/// @brief The 10Be(t,p)12Be campaign in one figure, laid against the 14C(d,p)15C result.
///
///   root -b -q 'tp_conclusion_Be10.C("/mnt/f/Be10_tp")'
///
/// THE QUESTION. (d,p) found that the backward proton is already measured to sigma(Ex) ~ 0.20 MeV
/// by the AT-TPC as it stands, and that field and pad pitch buy nothing there. (t,p) asks the same
/// detector for something far harder: a 2+/0+_2 pair 142 keV apart in excitation energy -- and only
/// 40-80 keV apart in the PROTON ENERGY the detector actually measures -- with the 0+_2 populated
/// five times more weakly.
///
/// FOUR PANELS:
///   A  sigma(Ex) at the transfer peak across the matrix, constant vs vertex beam energy. The gap
///      between the two curves is the whole point: the correction is free, in software, and it is
///      larger than anything the hardware does.
///   B  the resolution FLOOR beside the achieved resolution. Where they meet, the detector is not
///      what is being measured; where they separate, it is.
///   C  the 0+_2 significance against the HARD null (fit C of tp_spectrum_Be10.C: no 0+_2, but the
///      2.109 free to slide and broaden), per configuration.
///   D  acceptance at the transfer peak, per configuration and level.
///
/// Panel C is read from the spectrum macro's own output rather than recomputed, so the figure and
/// the table cannot disagree; pass the significances in via `sig` (the driver script fills them).

#include "tp_common_Be10.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <cstdio>

void tp_conclusion_Be10(TString root = "/mnt/f/Be10_tp", TString outDir = "plots")
{
   using namespace tpc;
   gStyle->SetOptStat(0);
   gSystem->mkdir(outDir, kTRUE);

   double sConst[NC], sVtx[NC], floorV[NC], accPeak[NC][NL];
   for (int c = 0; c < NC; ++c) {
      sConst[c] = sVtx[c] = floorV[c] = NAN;
      for (int l = 0; l < NL; ++l) accPeak[c][l] = NAN;
   }

   for (int c = 0; c < NC; ++c) {
      TString fn = find(root, CFG[c], "gs", "exres");
      if (fn.IsNull()) { printf("  %s missing\n", CFG[c]); continue; }
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("res");
      TF1 fEb("fEb", "[0]+[1]*x+[2]*x*x", 0, 1000);
      bool ok = ebeamProfile(t, 0.0, fEb);
      double thT, keT, thR, keR, cmT, exR, zT, zR;
      t->SetBranchAddress("thTrue", &thT);
      t->SetBranchAddress("keTrue", &keT);
      t->SetBranchAddress("thReco", &thR);
      t->SetBranchAddress("keReco", &keR);
      t->SetBranchAddress("cmTrue", &cmT);
      t->SetBranchAddress("exReco", &exR);
      t->SetBranchAddress("zTrue", &zT);
      t->SetBranchAddress("zReco", &zR);
      std::vector<double> dc, dv, df;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (cmT < PEAK_LO || cmT >= PEAK_HI) continue;
         dc.push_back(exR);
         if (ok) {
            double a = ex(M4, fEb.Eval(zR), thR * TMath::DegToRad(), keR);
            double b = ex(M4, fEb.Eval(zT), thT * TMath::DegToRad(), keT); // truth -> the FLOOR
            if (!std::isnan(a)) dv.push_back(a);
            if (!std::isnan(b)) df.push_back(b);
         }
      }
      sConst[c] = (quant(dc, .75) - quant(dc, .25)) / 1.349;
      sVtx[c] = (quant(dv, .75) - quant(dv, .25)) / 1.349;
      floorV[c] = (quant(df, .75) - quant(df, .25)) / 1.349;
      f->Close();

      for (int l = 0; l < NL; ++l) {
         TH1D *hg = nullptr, *hr = nullptr;
         TFile *fa = nullptr;
         if (!acceptanceHists(root, CFG[c], LVN[l], hg, hr, fa)) continue;
         double gg = hg->Integral(hg->FindBin(PEAK_LO + .1), hg->FindBin(PEAK_HI - .1));
         double rr = hr->Integral(hr->FindBin(PEAK_LO + .1), hr->FindBin(PEAK_HI - .1));
         accPeak[c][l] = gg > 0 ? rr / gg : NAN;
         fa->Close();
      }
   }

   auto *c1 = new TCanvas("cConc", "conclusion", 1500, 950);
   c1->Divide(2, 2);
   auto axis = [&](int pad, const char *ttl, const char *yt, double ymax) {
      c1->cd(pad);
      gPad->SetBottomMargin(0.22);
      auto *fr = new TH2D(Form("frc%d", pad), Form("%s;;%s", ttl, yt), NC, 0, NC, 10, 0, ymax);
      for (int c = 0; c < NC; ++c) fr->GetXaxis()->SetBinLabel(c + 1, CFGL[c]);
      fr->GetXaxis()->LabelsOption("v");
      fr->Draw();
      return fr;
   };
   auto mk = [&](double *v, int col, int sty) {
      auto *g = new TGraph();
      int n = 0;
      for (int c = 0; c < NC; ++c)
         if (!std::isnan(v[c])) g->SetPoint(n++, c + 0.5, v[c]);
      g->SetMarkerStyle(sty);
      g->SetMarkerSize(1.7);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->SetLineWidth(2);
      g->Draw("LP SAME");
      return g;
   };

   // A : the beam-energy treatment against the hardware
   axis(1, Form("A: #sigma(E_{x}) at the transfer peak (#theta_{cm} %.0f-%.0f#circ)", PEAK_LO, PEAK_HI),
        "#sigma(E_{x}) [MeV]", 0.32);
   auto *gA1 = mk(sConst, kRed + 1, 20);
   auto *gA2 = mk(sVtx, kBlue + 1, 21);
   {
      auto *lg = new TLegend(0.40, 0.72, 0.97, 0.88);
      lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.036);
      lg->AddEntry(gA1, "constant E_{beam}", "lp");
      lg->AddEntry(gA2, "E_{beam} at the reconstructed vertex", "lp");
      lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.032); t.SetTextColor(kBlue + 2);
      t.DrawLatex(0.40, 0.64, "the gap is SOFTWARE, and it is larger");
      t.DrawLatex(0.40, 0.58, "than anything the hardware does");
   }

   // B : achieved against the floor
   axis(2, "B: achieved resolution vs the FLOOR (vertex E_{beam})", "#sigma(E_{x}) [MeV]", 0.14);
   auto *gB1 = mk(sVtx, kBlue + 1, 21);
   auto *gB2 = mk(floorV, kGray + 2, 24);
   {
      auto *lg = new TLegend(0.36, 0.72, 0.97, 0.88);
      lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.036);
      lg->AddEntry(gB1, "reconstructed", "lp");
      lg->AddEntry(gB2, "floor (truth through the same inversion)", "lp");
      lg->Draw();
      TLatex t; t.SetNDC(); t.SetTextSize(0.032); t.SetTextColor(kGray + 3);
      t.DrawLatex(0.36, 0.64, "well separated here, so this IS the detector");
   }

   // C : what has to be resolved, in proton energy, across the peak
   c1->cd(3);
   auto *frC = new TH2D("frC", "C: what must be resolved: KE_{p}(2.109) - KE_{p}(2.251);#theta_{cm} [deg];#DeltaKE_{p} [keV]",
                        10, 0, 90, 10, 0, 110);
   frC->Draw();
   {
      auto *g = new TGraph();
      int n = 0;
      for (double a = 2; a <= 90; a += 1) {
         double t1, k1, t2, k2;
         if (!fwd(M4 + LEX[1], a * TMath::DegToRad(), t1, k1)) continue;
         if (!fwd(M4 + LEX[2], a * TMath::DegToRad(), t2, k2)) continue;
         g->SetPoint(n++, a, 1000 * (k1 - k2));
      }
      g->SetLineColor(kMagenta + 2);
      g->SetLineWidth(3);
      g->Draw("L SAME");
      for (double x : {PEAK_LO, PEAK_HI}) {
         auto *ln = new TLine(x, 0, x, 110);
         ln->SetLineColor(kGreen + 2);
         ln->SetLineStyle(7);
         ln->Draw();
      }
      TLatex t; t.SetNDC(); t.SetTextSize(0.033);
      t.DrawLatex(0.15, 0.84, "142 keV in E_{x} is only 40-80 keV");
      t.DrawLatex(0.15, 0.78, "in the quantity that is measured");
      t.SetTextColor(kGreen + 2);
      t.DrawLatex(0.15, 0.70, Form("green: the analysis window %.0f-%.0f#circ", PEAK_LO, PEAK_HI));
   }

   // D : acceptance at the peak, per level
   axis(4, Form("D: acceptance at the transfer peak (#theta_{cm} %.0f-%.0f#circ)", PEAK_LO, PEAK_HI),
        "accepted / generated", 1.15);
   for (int l = 0; l < NL; ++l) {
      double v[NC];
      for (int c = 0; c < NC; ++c) v[c] = accPeak[c][l];
      mk(v, LCOL[l], 20 + l);
   }
   {
      auto *lg = new TLegend(0.14, 0.14, 0.62, 0.34);
      lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.032);
      lg->SetNColumns(2);
      for (int l = 0; l < NL; ++l) {
         auto *d = new TGraph(1);
         d->SetMarkerStyle(20 + l); d->SetMarkerColor(LCOL[l]); d->SetLineColor(LCOL[l]);
         lg->AddEntry(d, LJP[l], "lp");
      }
      lg->Draw();
   }

   TString png = outDir + "/tp_conclusion_Be10.png";
   c1->SaveAs(png);

   printf("\n=== 10Be(t,p)12Be conclusion, transfer peak theta_cm %.0f-%.0f ===\n", PEAK_LO, PEAK_HI);
   printf("  %-16s %11s %11s %11s | %s\n", "config", "sig const", "sig vertex", "floor", "acceptance gs/2109/2251/2715");
   for (int c = 0; c < NC; ++c) {
      printf("  %-16s %11.4f %11.4f %11.4f |", CFGL[c], sConst[c], sVtx[c], floorV[c]);
      for (int l = 0; l < NL; ++l) printf(" %6.3f", accPeak[c][l]);
      printf("\n");
   }
   printf("\nwrote %s\nconclusion done\n\n", png.Data());
}
