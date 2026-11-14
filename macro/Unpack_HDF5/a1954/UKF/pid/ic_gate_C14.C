/// @file ic_gate_C14.C
/// @brief a1954 14C beam gate, step 2 of 2: the ion-chamber amplitude spectrum ALONE, with the
///        accepted window marked. PLAIN ROOT -- no ATTPCROOT, no raw data. It draws the histogram
///        ic_gate_scan_C14.C wrote; the a1954_C14_absnorm repo vendors that as
///        data/ic_gate_C14.root (14 production runs, TB 1050-1250, 390228 IC signals).
///
///        This replaces the left panel of the old two-panel ic_C14.C figure
///        (pid/plots/ic_C14_attpc_win.png, shipped as figures/ic_gate_C14.png in the
///        write-up) with vertical lines at icLo / icHi.
///
///        The defaults are the PRODUCTION gate of pipeline/gate_events_C14.C
///        (amplitude 950-1350, TB 1050-1250) -- note the old figure was made with the
///        wider TB search window 1000-1350, so the two spectra are not byte-identical.
///        The pile-up (exactly-one-pulse) condition of the production gate is NOT applied
///        here: this figure shows the amplitude cut, which is what the lines mark.
///
///   root -b -q 'ic_gate_C14.C'                  // production window on the vendored spectrum
///   root -b -q 'ic_gate_C14.C(900,1400,"_wide")' // another amplitude window, same spectrum
///   root -b -q 'ic_gate_C14.C(950,1350,"_r55","run_0055")'  // after ic_gate_scan_C14.C(...,"_r55")
///
/// The spectrum is looked up as ic_gate_C14<outTag>.root in plots/ beside this macro (where the
/// scan writes it, ATTPCROOT tree), then in ../data/ (the absnorm repo). It is drawn only if its
/// stamp matches the runs and TB window asked for; the amplitude window is applied here and needs
/// no rescan. The png goes to plots/ beside this macro.
void ic_gate_C14(double icLo = 950, double icHi = 1350, TString outTag = "",
                 TString runsCSV = "run_0055,run_0056,run_0057,run_0058,run_0059,run_0060,run_0061,"
                                   "run_0062,run_0063,run_0064,run_0065,run_0066,run_0068,run_0069",
                 Int_t icTbLo = 1050, Int_t icTbHi = 1250)
{
   gStyle->SetOptStat(0);

   TString dir = gSystem->DirName(__FILE__);
   TString name = "ic_gate_C14" + outTag + ".root";
   TString inFile = dir + "/plots/" + name;
   if (gSystem->AccessPathName(inFile))
      inFile = dir + "/../data/" + name;
   if (gSystem->AccessPathName(inFile)) {
      printf("\033[1;31mmissing %s (looked in plots/ and ../data/) -- run ic_gate_scan_C14.C first\033[0m\n",
             name.Data());
      return;
   }

   // ---- the spectrum: drawn only if it was made from the runs and TB window asked for ------
   TString stamp = TString::Format("%s|tb%d-%d", runsCSV.Data(), icTbLo, icTbHi);
   TFile *fc = TFile::Open(inFile);
   auto *st = fc ? (TNamed *)fc->Get("stamp") : nullptr;
   auto *hc = fc ? (TH1D *)fc->Get("hAmp") : nullptr;
   if (!st || !hc || stamp != st->GetTitle()) {
      printf("\033[1;31m%s was not made from these inputs -- rerun ic_gate_scan_C14.C\033[0m\n  want: %s\n  have: %s\n",
             inFile.Data(), stamp.Data(), st ? st->GetTitle() : "(none)");
      return;
   }
   TH1D *hAmp = (TH1D *)hc->Clone("hAmp");
   hAmp->SetDirectory(nullptr);
   fc->Close();

   // ---- what the window keeps ---------------------------------------------------------
   double all = hAmp->Integral(1, hAmp->GetNbinsX());
   double in = hAmp->Integral(hAmp->FindBin(icLo), hAmp->FindBin(icHi) - 1);
   printf("\nIC amplitude (%s): entries=%.0f  inside [%.0f,%.0f]=%.0f (%.1f %%)\n", inFile.Data(), all, icLo, icHi,
          in, all > 0 ? 100. * in / all : 0.);

   // ---- the plot ----------------------------------------------------------------------
   TCanvas *c = new TCanvas("cIcGate", "IC beam gate", 700, 550);
   c->SetLogy();
   c->SetTicks(1, 1);
   hAmp->SetLineColor(kBlack);
   hAmp->SetLineWidth(2);
   hAmp->GetYaxis()->SetTitleOffset(1.25);
   // frame on the data: the 0-4096 ADC range is mostly empty (saturation sits near 2050)
   int last = hAmp->FindLastBinAbove(0.);
   if (last > 0)
      hAmp->GetXaxis()->SetRangeUser(0., hAmp->GetXaxis()->GetBinUpEdge(last) * 1.05);
   hAmp->Draw("hist");
   gPad->Update();

   double y1 = TMath::Power(10, gPad->GetUymin());
   double y2 = TMath::Power(10, gPad->GetUymax());
   for (double x : {icLo, icHi}) {
      auto *l = new TLine(x, y1, x, y2);
      l->SetLineColor(kRed + 1);
      l->SetLineStyle(2);
      l->SetLineWidth(2);
      l->Draw();
   }
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.033);
   tx->SetTextColor(kRed + 1);
   tx->DrawLatex(0.57, 0.87, TString::Format("accepted window %.0f#minus%.0f ADC", icLo, icHi));
   tx->SetTextColor(kBlack);
   tx->DrawLatex(0.57, 0.81, TString::Format("%.1f %% of IC signals", all > 0 ? 100. * in / all : 0.));

   TString plotDir = dir + "/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE); // the absnorm copy has no plots/ dir; SaveAs fails mutely without this
   TString png = plotDir + "ic_gate_C14" + outTag + ".png";
   c->SaveAs(png);
   if (gSystem->AccessPathName(png)) {
      printf("\033[1;31mSaveAs did not write %s\033[0m\n", png.Data());
      return;
   }
   printf("saved %s\n", png.Data());
}
