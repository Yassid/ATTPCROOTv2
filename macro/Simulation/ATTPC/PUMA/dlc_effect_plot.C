/// @file dlc_effect_plot.C
/// @brief DLC resistive-anode effect on PUMA momentum reconstruction: overlay the
///        PRA circle-fit p_T = 0.3 B R distribution with and without DLC charge
///        dispersion. DLC (1.35 MOhm/sq) spreads charge ~4x but, reconstructed with
///        pad-centre hits, biases the circle radius LOW (no sub-pad gain realised).
/// Run: root -b -q 'dlc_effect_plot.C'
void dlc_effect_plot(TString off = "data/dlc_off.root", TString on = "data/dlc_on.root",
                     TString out = "/Users/quantumlab/fair_install/puma_slides/figs/dlc_effect.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const double B = 4.0, pT0 = 0.3749;
   auto fill = [&](TString f, TH1F *h) {
      TFile ff(f); auto *t = (TTree *)ff.Get("cbmsim"); TClonesArray *pat = nullptr; t->SetBranchAddress("AtPatternEvent", &pat);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); if (!pat->GetEntries()) continue;
         for (auto &tr : ((AtPatternEvent *)pat->At(0))->GetTrackCand()) { double R = tr.GetGeoRadius();
            if (R > 0 && R < 1e4) { double p = 0.3 * B * R / 1000.0; if (p > 0.05 && p < 0.7) h->Fill(p); } } }
   };
   auto *hOff = new TH1F("hOff", "", 60, 0.05, 0.65);
   auto *hOn = new TH1F("hOn", "", 60, 0.05, 0.65);
   fill(off, hOff); fill(on, hOn);

   auto *c = new TCanvas("c", "", 820, 560); c->SetLeftMargin(0.12); c->SetBottomMargin(0.13);
   hOff->SetLineColor(kAzure + 2); hOff->SetLineWidth(3);
   hOn->SetLineColor(kRed + 1); hOn->SetLineWidth(3); hOn->SetLineStyle(2);
   hOff->SetTitle(";reconstructed p_{T} = 0.3 B R  [GeV/c];tracks");
   hOff->SetMaximum(1.3 * std::max(hOff->GetMaximum(), hOn->GetMaximum()));
   hOff->Draw("hist"); hOn->Draw("hist same");
   auto *lt = new TLine(pT0, 0, pT0, hOff->GetMaximum()); lt->SetLineColor(kGreen + 2); lt->SetLineWidth(2); lt->SetLineStyle(9); lt->Draw();
   auto *lg = new TLegend(0.14, 0.72, 0.55, 0.89); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hOff, "no DLC (bias #minus5%, #sigma 17%)", "l");
   lg->AddEntry(hOn, "DLC 1.35 M#Omega/#Box (bias #minus30%, #sigma 23%)", "l");
   lg->AddEntry(lt, "truth p_{T} = 0.375 GeV/c", "l");
   lg->Draw();
   auto *tt = new TLatex(); tt->SetNDC(); tt->SetTextFont(62); tt->SetTextSize(0.042);
   tt->DrawLatex(0.12, 0.93, "DLC charge dispersion biases the reconstructed radius low");
   auto *tn = new TLatex(); tn->SetNDC(); tn->SetTextFont(62); tn->SetTextSize(0.028); tn->SetTextColor(kRed + 1);
   tn->DrawLatex(0.14, 0.66, "pad-centre reco of dispersed charge #Rightarrow no sub-pad gain");
   c->SaveAs(out);
   printf("DLC_EFFECT_DONE\n");
}
