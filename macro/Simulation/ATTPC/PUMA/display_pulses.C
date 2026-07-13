/// @file display_pulses.C
/// @brief Example pad pulses (ADC vs time bucket) with a fit to the pulse and
///        the PSA-extracted hit (peak time / amplitude) marked. Illustrates the
///        pulse-shape analysis step.
/// Run: root -b -q 'display_pulses.C("data/reco_pid_base.root")'
void display_pulses(TString file = "data/reco_pid_base.root",
                    TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pulses_psa.png")
{
   gStyle->SetOptStat(0); gStyle->SetOptFit(0); gStyle->SetTextFont(62);
   gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   TFile f(file); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *raw = nullptr;
   t->SetBranchAddress("AtRawEvent", &raw);

   // Find an event with several well-formed pulses (high-max pads).
   Long64_t bestEv = 0; int bestN = 0;
   for (Long64_t i = 0; i < std::min((Long64_t)60, t->GetEntries()); i++) {
      t->GetEntry(i);
      auto *ev = (AtRawEvent *)(raw ? raw->At(0) : nullptr);
      if (!ev) continue;
      int n = 0;
      for (const auto &p : ev->GetPads()) { double mx = 0; for (double a : p->GetADC()) mx = std::max(mx, a); if (mx > 50) n++; }
      if (n > bestN) { bestN = n; bestEv = i; }
   }
   t->GetEntry(bestEv);
   auto *ev = (AtRawEvent *)raw->At(0);

   // Rank pads by peak amplitude, take the top 4.
   std::vector<std::pair<double, const AtPad *>> pads;
   for (const auto &p : ev->GetPads()) { double mx = 0; for (double a : p->GetADC()) mx = std::max(mx, a); pads.emplace_back(mx, p.get()); }
   std::sort(pads.begin(), pads.end(), [](auto &a, auto &b) { return a.first > b.first; });

   auto *c = new TCanvas("c", "", 1100, 850); c->Divide(2, 2);
   int npan = std::min((size_t)4, pads.size());
   for (int k = 0; k < npan; k++) {
      c->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13);
      const AtPad *pad = pads[k].second;
      const auto &adc = pad->GetADC();
      // peak
      int pk = 0; double pkv = 0;
      for (int i = 0; i < 512; i++) if (adc[i] > pkv) { pkv = adc[i]; pk = i; }
      int lo = std::max(0, pk - 40), hi = std::min(511, pk + 60);
      auto *g = new TGraph();
      for (int i = lo; i <= hi; i++) g->SetPoint(g->GetN(), i, adc[i]);
      g->SetTitle(Form("pad %d;time bucket;ADC", pad->GetPadNum()));
      g->SetLineColor(kAzure + 2); g->SetLineWidth(2); g->SetMarkerStyle(20); g->SetMarkerSize(0.4);
      g->SetMarkerColor(kAzure + 3); g->Draw("APL");
      // fit the pulse peak with a Gaussian (illustrates the PSA shape fit)
      auto *fit = new TF1(Form("fit%d", k), "gaus", pk - 15, pk + 15);
      fit->SetParameters(pkv, pk, 6); fit->SetLineColor(kRed + 1); fit->SetLineWidth(3);
      g->Fit(fit, "QR");
      fit->Draw("same");
      // PSA hit: peak time / amplitude
      auto *m = new TMarker(fit->GetParameter(1), fit->GetParameter(0), 29);
      m->SetMarkerColor(kOrange + 8); m->SetMarkerSize(2.4); m->Draw();
      auto *lv = new TLine(fit->GetParameter(1), 0, fit->GetParameter(1), fit->GetParameter(0));
      lv->SetLineColor(kOrange + 8); lv->SetLineStyle(2); lv->SetLineWidth(2); lv->Draw();
      if (k == 0) {
         auto *lg = new TLegend(0.5, 0.72, 0.88, 0.9); lg->SetTextFont(62); lg->SetTextSize(0.05);
         lg->AddEntry(g, "pad ADC pulse", "l"); lg->AddEntry(fit, "PSA pulse fit", "l");
         lg->AddEntry(m, "PSA hit (t, A)", "p"); lg->Draw();
      }
   }
   c->SaveAs(out);
   printf("PULSES_DONE event %lld  (%d signal pads)\n", bestEv, bestN);
}
