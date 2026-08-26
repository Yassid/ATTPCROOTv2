/// @file ic_fraction_C14.C
/// @brief Measure f_beam -- the fraction of ion-chamber counts that are 14C -- from the IC spectrum.
///
/// WHY. The scalers give a luminosity that owes nothing to an optical model:
///     L = ic_sca * f_beam * livetime * n_target / 1e27
/// Every term is measured except f_beam, and f_beam multiplies L directly, so it is the whole
/// uncertainty of that route. a1975 measured its own on the H2 block (0.613); assuming that number
/// here would defeat the purpose, since a1954 is a different beam at a different separator setting.
///
/// WHAT IS MEASURED, and what it is not. This is the fraction of TRIGGERED events whose IC pulse
/// falls in the analysis gate -- amplitude in [icLo,icHi] over TB [tbLo,tbHi], exactly one pulse,
/// which is verbatim the gate pipeline/gate_events_C14.C applies. It is the composition of the
/// beam as the IC ADC sees it. It equals the scaler f_beam only if the TPC trigger does not
/// prefer one species over another and if the scaler discriminator counts the same population the
/// ADC digitises. NEITHER IS VERIFIED HERE -- both are stated in the output, not buried.
///
/// MEASURED, runs 55-66,68,69 (the C14 analysis list), 2026-08-27 -- plots/ is gitignored so the
/// result is recorded here. 396102 triggered events with an IC trace; 322116 with exactly one
/// pulse (0.8132); 281714 also inside 950-1350, so
///
///     f_beam = 0.7112
///
/// It is NOT flat across the block: 0.580 on run_0055 rising to 0.804 on run_0069. Taking one run
/// as representative gives 0.58 and a luminosity 22% low -- do not shortcut this.
///
///   root -b -q 'ic_fraction_C14.C("run_0055,run_0058,run_0060")'
///   root -b -q 'ic_fraction_C14.C("")'    // all runs found in inDir
namespace icf {
int CountPulses(const std::vector<Double_t> &adc, double thr, int tbLo, int tbHi)
{
   int n = 0; bool in = false;
   for (int i = tbLo; i <= tbHi && i < (int)adc.size(); ++i) {
      if (!in && adc[i] > thr) { in = true; ++n; }
      else if (in && adc[i] < 0.5 * thr) in = false;
   }
   return n;
}
} // namespace icf

void ic_fraction_C14(TString runsCSV = "", TString inDir = "/mnt/f/a1954_C14_reco_hdb_slim/",
                     double icLo = 950, double icHi = 1350, Int_t tbLo = 1050, Int_t tbHi = 1250,
                     double peakThr = 200, Int_t pkTbLo = 800, Int_t pkTbHi = 1500)
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   std::vector<TString> runs;
   if (runsCSV != "") { TObjArray *a = runsCSV.Tokenize(","); for (auto o : *a) runs.push_back(((TObjString *)o)->String()); }
   else { void *d = gSystem->OpenDirectory(inDir); const char *e;
          while ((e = gSystem->GetDirEntry(d))) { TString s(e);
             if (s.EndsWith("_FRIB.root")) runs.push_back(s(0, s.Length() - 10)); }
          gSystem->FreeDirectory(d); std::sort(runs.begin(), runs.end()); }

   auto *hIC = new TH1D("hIC", "IC amplitude, all triggered events;IC max [ADC];events", 400, 0, 4000);
   auto *hIC1 = new TH1D("hIC1", "", 400, 0, 4000);          // single-pulse only
   long long nTot = 0, nOne = 0, nGate = 0;
   printf("\n  %-10s %10s %10s %10s %8s\n", "run", "events", "1 pulse", "in gate", "frac");
   for (auto &r : runs) {
      TString ff = inDir + r + "_FRIB.root";
      TFile *f = TFile::Open(ff);
      TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t || !t->GetEntries()) { if (f) f->Close(); continue; }
      TClonesArray *ra = nullptr; t->SetBranchAddress("AtRawEvent", &ra);
      long long a = 0, b = 0, c = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         auto *raw = (AtRawEvent *)ra->At(0);
         if (!raw || raw->GetGenTraces().empty()) continue;
         const auto &adc = raw->GetGenTraces()[0]->GetADC();
         double mx = 0;
         for (int k = tbLo; k <= tbHi && k < (int)adc.size(); ++k) mx = std::max(mx, adc[k]);
         ++a; hIC->Fill(mx);
         int np = icf::CountPulses(adc, peakThr, pkTbLo, pkTbHi);
         if (np == 1) { ++b; hIC1->Fill(mx); if (mx >= icLo && mx <= icHi) ++c; }
      }
      printf("  %-10s %10lld %10lld %10lld %8.4f\n", r.Data(), a, b, c, a ? (double)c / a : 0);
      nTot += a; nOne += b; nGate += c;
      f->Close();
   }
   if (!nTot) { printf("\n  no events read from %s\n\n", inDir.Data()); return; }
   printf("\n  ===== TOTALS =====\n");
   printf("    triggered events with an IC trace   %12lld\n", nTot);
   printf("    exactly one pulse                   %12lld   (%.4f)\n", nOne, (double)nOne / nTot);
   printf("    and IC max in [%.0f,%.0f]            %12lld   (%.4f)\n", icLo, icHi, nGate, (double)nGate / nTot);
   printf("\n    f_beam (measured, this definition)  %12.4f\n", (double)nGate / nTot);
   printf("\n  \033[1;33mThis is the fraction of TRIGGERED events, not of scaler counts. It equals the\n"
          "  scaler f_beam only if the trigger is species-blind and the scaler discriminator counts\n"
          "  the same population the ADC digitises. Neither is checked here.\033[0m\n");

   {
      TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
      std::ofstream o((here + "/plots/ic_fraction.txt").Data());
      o << "# f_beam  nGate  nTot  (fraction of TRIGGERED events inside the analysis IC gate)\n";
      o << (double)nGate / nTot << " " << nGate << " " << nTot << "\n";
      printf("  wrote plots/ic_fraction.txt\n");
   }
   auto *c1 = new TCanvas("cicf", "", 1200, 500); c1->Divide(2, 1);
   c1->cd(1); gPad->SetLogy(); gPad->SetGridx();
   hIC->SetLineColor(kGray + 2); hIC->SetLineWidth(2); hIC->Draw("hist");
   hIC1->SetLineColor(kBlue + 1); hIC1->SetLineWidth(2); hIC1->Draw("hist same");
   auto *l1 = new TLine(icLo, 0.5, icLo, hIC->GetMaximum()); l1->SetLineColor(kRed + 1);
   l1->SetLineWidth(2); l1->SetLineStyle(2); l1->Draw();
   auto *l2 = new TLine(icHi, 0.5, icHi, hIC->GetMaximum()); l2->SetLineColor(kRed + 1);
   l2->SetLineWidth(2); l2->SetLineStyle(2); l2->Draw();
   auto *lg = new TLegend(0.45, 0.72, 0.90, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hIC, "all triggered", "l"); lg->AddEntry(hIC1, "exactly one pulse", "l");
   lg->AddEntry(l1, Form("gate %.0f-%.0f", icLo, icHi), "l"); lg->Draw();
   c1->cd(2); gPad->SetGridx(); gPad->SetGridy();
   hIC1->GetXaxis()->SetRangeUser(0, 2500);
   auto *h2 = (TH1D *)hIC1->Clone("h2"); h2->SetTitle("single-pulse events, linear;IC max [ADC];events");
   h2->Draw("hist"); l1->Draw(); l2->Draw();
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
   tx.DrawLatex(0.45, 0.80, Form("f_{beam} = %.3f", (double)nGate / nTot));
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "08_ic_fraction.png");
   printf("  wrote %s08_ic_fraction.png\n\n", out.Data());
}
