/// @file plots_matrix_3Hed.C
/// @brief The three figures of the 2026-09-04 campaign: the field x pad-plane resolution matrix,
///        the excitation spectra that show whether the states separate, and the diffusion scan.
///
///   root -b -q 'plots_matrix_3Hed.C()'
///
/// Every number here comes from Ar46::Collect in ex_core_3Hed.h -- the SAME inversion the table
/// macro uses, so a figure cannot quietly disagree with its own table. Do not reimplement the
/// kinematics in this file.
///
/// WIDTHS ARE QUOTED AS THE FWHM AN EQUALLY WIDE GAUSSIAN WOULD HAVE, converted from the IQR
/// (FWHM = 1.7456 x IQR), never measured by walking out from the peak. These distributions are a
/// narrow core on broad tails, and a half-maximum walk returned 0.765 against 1.800 MeV for two
/// histograms differing only by a constant shift. The conversion is exact only for a Gaussian, so
/// it is a COMPARABLE width, not a claim that the lineshape is Gaussian -- which is also why the
/// spectra of figure 2 are drawn: they show the lineshape the number is standing in for.
///
/// FIGURE 2 GIVES THE THREE STATES EQUAL AREA. It answers "does the detector separate them",
/// which is a resolution question; it is NOT a predicted spectrum. The real relative yields need
/// the DWBA (dwba/ar46_dwba_ptolemy.txt), and the l = 3 curve there is on a single-particle scale
/// with no spectroscopic factor, so a weighted spectrum would carry an invented normalisation for
/// exactly one of the three peaks. Equal areas make the assumption visible instead of hiding it.
#include "ex_core_3Hed.h"

static void styleGraph(TGraph *g, int col, int mrk)
{
   g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mrk);
   g->SetLineWidth(3); g->SetMarkerSize(1.4);
}

void plots_matrix_3Hed(TString outDir = "plots", Double_t dThetaMax = 10.0, Double_t driftLength = 100.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->mkdir(outDir, kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

   const TString A = "/mnt/f/ar46_3hed_OLD_2.85T_placeholder"; // 2.85 T generation
   const TString B = "/mnt/f/ar46_3hed_gen_B39";               // 3.9 T generation

   const int NC = 4;
   TString cdir[NC] = {"/mnt/f/ar46_3hed_mx_B285_attpc", "/mnt/f/ar46_3hed_mx_B39_attpc",
                       "/mnt/f/ar46_3hed_mx_B285_2mm",   "/mnt/f/ar46_3hed_mx_B39_2mm"};
   TString csim[NC] = {A, B, A, B};
   TString cnam[NC] = {"2.85 T, AT-TPC pads", "3.9 T, AT-TPC pads", "2.85 T, 2 mm pads", "3.9 T, 2 mm pads"};
   int ccol[NC] = {kGray + 2, kAzure + 2, kRed + 1, kOrange + 8};
   int cmrk[NC] = {24, 25, 20, 21};

   const int NS = 3;
   TString stag[NS] = {"gs_s3001", "360_s3011", "2020_s3021"};
   TString snam[NS] = {"^{47}K g.s.  1/2^{+}  (1s_{1/2}, #font[12]{l} = 0)",
                       "^{47}K 0.36  3/2^{+}  (0d_{3/2}, #font[12]{l} = 2)",
                       "^{47}K 2.02  7/2^{-}  (0f_{7/2}, #font[12]{l} = 3)"};
   double sex[NS] = {0.0, 0.36, 2.02};
   // Short forms for figure 2: the full names above carry the orbital and overflow the legend box,
   // which silently CLIPS the FWHM that is the point of the entry.
   TString sshort[NS] = {"g.s.  1/2^{+}", "0.36  3/2^{+}", "2.02  7/2^{-}"};

   const int NB = 7;
   double edge[NB + 1] = {60, 70, 80, 90, 100, 110, 125, 145};

   // ---------------------------------------------------------------------------------------
   // Collect once per (configuration, state) and keep it: each Collect reads a ~90 MB fit file
   // plus a ~180 MB sim file, so re-reading per figure would triple the cost for nothing.
   // ---------------------------------------------------------------------------------------
   Ar46::Sample S[NC][NS];
   for (int c = 0; c < NC; ++c)
      for (int s = 0; s < NS; ++s) {
         TObjArray *t = stag[s].Tokenize(",");
         S[c][s] = Ar46::Collect(cdir[c], csim[c], t, kTRUE, dThetaMax, -1.0, driftLength, cnam[c]);
         printf("  %-22s %-12s  %ld/%ld fits\n", cnam[c].Data(), stag[s].Data(), S[c][s].nFit, S[c][s].nTruth);
      }

   // =======================================================================================
   // FIGURE 1 -- resolution against theta_lab, one panel per state
   // =======================================================================================
   TCanvas *c1 = new TCanvas("c1", "resolution matrix", 1650, 560);
   c1->Divide(3, 1, 0.001, 0.001);
   TLegend *lg = nullptr;
   for (int s = 0; s < NS; ++s) {
      c1->cd(s + 1);
      gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
      TH1F *fr = gPad->DrawFrame(60, 0.10, 145, 20);
      fr->SetTitle(Form("%s;#theta_{lab} of the deuteron [deg];FWHM(E_{x}) [MeV]", snam[s].Data()));
      fr->GetYaxis()->SetTitleOffset(1.25);

      // the proposal's own target, drawn so every panel is read against it
      TLine *goal = new TLine(60, 0.350, 145, 0.350);
      goal->SetLineColor(kGreen + 2); goal->SetLineStyle(2); goal->SetLineWidth(3);
      goal->Draw();
      if (s == 0) {
         TLatex *gt = new TLatex(62, 0.375, "proposal goal, 350 keV");
         gt->SetTextColor(kGreen + 2); gt->SetTextSize(0.040); gt->Draw();
      }

      for (int c = 0; c < NC; ++c) {
         auto *g = new TGraph();
         for (int b = 0; b < NB; ++b) {
            std::vector<double> v;
            for (size_t j = 0; j < S[c][s].ex.size(); ++j)
               if (S[c][s].thetaTrue[j] >= edge[b] && S[c][s].thetaTrue[j] < edge[b + 1]) v.push_back(S[c][s].ex[j]);
            if (v.size() < 8) continue;
            g->SetPoint(g->GetN(), 0.5 * (edge[b] + edge[b + 1]), Ar46::IqrToFwhm(Ar46::IqrOf(v)));
         }
         styleGraph(g, ccol[c], cmrk[c]);
         g->Draw("PL same");
         if (s == 0) {
            if (!lg) {
               lg = new TLegend(0.34, 0.63, 0.95, 0.89);
               lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
            }
            lg->AddEntry(g, cnam[c], "pl");
         }
      }
      if (s == 0) lg->Draw();
   }
   c1->cd(2);
   TLatex *note = new TLatex();
   note->SetNDC(); note->SetTextSize(0.034); note->SetTextColor(kAzure + 3);
   note->DrawLatex(0.17, 0.20, "3.9 T collapses beyond 110#circ at BOTH pad planes");
   c1->SaveAs(outDir + "/mx_resolution.png");

   // =======================================================================================
   // FIGURE 2 -- the spectra, in the window where the 2 mm plane is best
   // =======================================================================================
   const double wLo = 100, wHi = 125;
   TCanvas *c2 = new TCanvas("c2", "spectra", 1450, 580);
   c2->Divide(2, 1, 0.001, 0.001);
   int show[2] = {0, 2}; // AT-TPC vs 2 mm, both at 2.85 T
   int scol[NS] = {kBlack, kRed + 1, kAzure + 2};
   for (int k = 0; k < 2; ++k) {
      c2->cd(k + 1);
      gPad->SetGridx(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.15); gPad->SetTopMargin(0.09);
      int c = show[k];
      THStack *st = new THStack("st", "");
      TLegend *l2 = new TLegend(0.50, 0.64, 0.97, 0.88);
      l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.033);
      double ymax = 0;
      TH1D *hh[NS];
      for (int s = 0; s < NS; ++s) {
         hh[s] = new TH1D(Form("h%d%d", k, s), "", 120, -2, 4);
         for (size_t j = 0; j < S[c][s].ex.size(); ++j)
            if (S[c][s].thetaTrue[j] >= wLo && S[c][s].thetaTrue[j] < wHi) hh[s]->Fill(S[c][s].ex[j]);
         // EQUAL AREA -- this figure is about separation, not yield. See the header.
         if (hh[s]->Integral() > 0) hh[s]->Scale(1.0 / hh[s]->Integral());
         hh[s]->SetLineColor(scol[s]); hh[s]->SetLineWidth(3);
         ymax = std::max(ymax, hh[s]->GetMaximum());
      }
      TH1F *fr = gPad->DrawFrame(-2, 0, 4, 1.30 * ymax);
      fr->SetTitle(Form("%s,  #theta_{lab} %.0f-%.0f#circ;E_{x}(^{47}K) [MeV];fraction of the state's tracks",
                        cnam[c].Data(), wLo, wHi));
      fr->GetYaxis()->SetTitleOffset(1.35);
      for (int s = 0; s < NS; ++s) {
         hh[s]->Draw("hist same");
         std::vector<double> v;
         for (size_t j = 0; j < S[c][s].ex.size(); ++j)
            if (S[c][s].thetaTrue[j] >= wLo && S[c][s].thetaTrue[j] < wHi) v.push_back(S[c][s].ex[j]);
         l2->AddEntry(hh[s], Form("%s   FWHM %.2f MeV", sshort[s].Data(), Ar46::IqrToFwhm(Ar46::IqrOf(v))), "l");
         TLine *tr = new TLine(sex[s], 0, sex[s], 1.30 * ymax);
         tr->SetLineColor(scol[s]); tr->SetLineStyle(3); tr->SetLineWidth(1); tr->Draw();
      }
      l2->Draw();
      TLatex *eq = new TLatex(); eq->SetNDC(); eq->SetTextSize(0.030); eq->SetTextColor(kGray + 3);
      eq->DrawLatex(0.16, 0.045, "equal area per state: separation, not predicted yield");
   }
   c2->SaveAs(outDir + "/mx_spectra.png");

   // =======================================================================================
   // FIGURE 3 -- the transverse-diffusion scan
   // =======================================================================================
   const int ND = 5;
   TString ddir[ND] = {"/mnt/f/ar46_3hed_ct05", "/mnt/f/ar46_3hed_mb_B285", "/mnt/f/ar46_3hed_ct20",
                       "/mnt/f/ar46_3hed_ct30", "/mnt/f/ar46_3hed_ct40"};
   double dmul[ND] = {0.5, 1.0, 2.0, 3.0, 4.0};        // multiplier on D_T
   double dsig[ND] = {0.63, 1.25, 2.50, 3.75, 5.00};   // sigma_T at the full 100 cm drift [mm]
   int dbin[3] = {4, 5, 6};                            // 100-110, 110-125, 125-145
   int dcol[3] = {kAzure + 2, kRed + 1, kBlack};
   int dmrk[3] = {21, 20, 22};

   TCanvas *c3 = new TCanvas("c3", "diffusion scan", 820, 620);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
   TH1F *fr3 = gPad->DrawFrame(0.3, 0.25, 4.4, 2.0);
   fr3->SetTitle("Transverse diffusion at 2.85 T, AT-TPC pads, ^{47}K g.s.;"
                 "D_{T} / D_{T}(3He+5%CO_{2} at 140 V/cm);FWHM(E_{x}) [MeV]");
   fr3->GetYaxis()->SetTitleOffset(1.25);
   TLine *goal3 = new TLine(0.3, 0.350, 4.4, 0.350);
   goal3->SetLineColor(kGreen + 2); goal3->SetLineStyle(2); goal3->SetLineWidth(3); goal3->Draw();
   TLegend *l3 = new TLegend(0.46, 0.68, 0.92, 0.90);
   l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.035);
   for (int q = 0; q < 3; ++q) {
      auto *g = new TGraph();
      for (int d = 0; d < ND; ++d) {
         TObjArray *t = TString("gs_s3001").Tokenize(",");
         Ar46::Sample Sd = Ar46::Collect(ddir[d], A, t, kTRUE, dThetaMax, -1.0, driftLength, "ct");
         std::vector<double> v;
         for (size_t j = 0; j < Sd.ex.size(); ++j)
            if (Sd.thetaTrue[j] >= edge[dbin[q]] && Sd.thetaTrue[j] < edge[dbin[q] + 1]) v.push_back(Sd.ex[j]);
         if (v.size() < 8) continue;
         g->SetPoint(g->GetN(), dmul[d], Ar46::IqrToFwhm(Ar46::IqrOf(v)));
      }
      styleGraph(g, dcol[q], dmrk[q]);
      g->Draw("PL same");
      l3->AddEntry(g, Form("#theta_{lab} %.0f-%.0f#circ", edge[dbin[q]], edge[dbin[q] + 1]), "pl");
   }
   l3->Draw();
   TLatex *t3 = new TLatex(); t3->SetNDC(); t3->SetTextSize(0.033); t3->SetTextColor(kGray + 3);
   t3->DrawLatex(0.17, 0.235, "#sigma_{T} over the 1 m drift: 0.63 / 1.25 / 2.50 / 3.75 / 5.00 mm,");
   t3->DrawLatex(0.17, 0.190, "against 8 #times 12 mm pads. The optimum is real: it turns over at 4#times.");
   c3->SaveAs(outDir + "/mx_coeft.png");

   printf("\nwrote %s/mx_resolution.png, %s/mx_spectra.png, %s/mx_coeft.png\n", outDir.Data(), outDir.Data(),
          outDir.Data());
}
