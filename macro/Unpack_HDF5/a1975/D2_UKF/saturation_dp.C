/// @file saturation_dp.C
/// @brief ADC saturation check for the a1975 D2 PID: hit amplitudes (AtHit::GetCharge, the fitted peak of
///        AtPSAMultiFit) for all hits, and for hits on slow-deuteron STREAK tracks (3D Brho 0.11-0.16,
///        sqrt(dE/dx) > 42) and on band tracks at the same Brho. Clipping shows as a pile-up below the
///        12-bit ceiling (~3500-4000 after pedestal subtraction). Also the fraction of each track's
///        charge carried by hits above a threshold.
///   root -l -b -q 'saturation_dp.C("/mnt/f/a1975/reco_d2_dv1104/run_0016_multifit_reco.root","plots/saturation_dp.png")'
void saturation_dp(TString reco, TString png, Long64_t maxEv = 1000000000)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gErrorIgnoreLevel = kError;
   gStyle->SetOptStat(0);
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   AtTools::AtHelix3D h3;
   auto *hAll = new TH1D("hAll", "hit amplitude;amplitude [ADC];hits", 300, 0, 6000);
   auto *hStr = new TH1D("hStr", "hit amplitude;amplitude [ADC];hits", 300, 0, 6000);
   auto *hBnd = new TH1D("hBnd", "hit amplitude;amplitude [ADC];hits", 300, 0, 6000);
   long nStr = 0, nBnd = 0;
   double fracS = 0, fracB = 0;
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      t->GetEntry(ev);
      if (!pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         for (auto &h : H) hAll->Fill(h->GetCharge());
         if (H.size() < 30) continue;
         auto sr = spy.Estimate(const_cast<AtTrack &>(tr));
         if (!sr.valid) continue;
         auto hr = h3.Estimate(tr, sr.direction);
         if (!hr.valid || hr.spacing >= 6 || hr.brho < 0.11 || hr.brho > 0.16) continue;
         bool st = hr.sqrtdEdx > 42, bd = hr.sqrtdEdx > 28 && hr.sqrtdEdx < 38;
         if (!st && !bd) continue;
         double qt = 0, qh = 0;
         for (auto &h : H) {
            (st ? hStr : hBnd)->Fill(h->GetCharge());
            qt += h->GetCharge();
            if (h->GetCharge() > 3000) qh += h->GetCharge();
         }
         if (st) { ++nStr; fracS += qh / qt; } else { ++nBnd; fracB += qh / qt; }
      }
   }
   auto tail = [](TH1D *h, double x) { return h->Integral(h->FindBin(x), h->GetNbinsX() + 1) / h->Integral(0, h->GetNbinsX() + 1); };
   printf("hits above 3000 / 3500 / 4000 ADC:  all %.4f / %.4f / %.4f | streak %.4f / %.4f / %.4f | band %.4f / %.4f / %.4f\n",
          tail(hAll, 3000), tail(hAll, 3500), tail(hAll, 4000), tail(hStr, 3000), tail(hStr, 3500), tail(hStr, 4000),
          tail(hBnd, 3000), tail(hBnd, 3500), tail(hBnd, 4000));
   printf("charge fraction carried by hits > 3000 ADC: streak tracks %.3f (n=%ld), band tracks %.3f (n=%ld)\n",
          nStr ? fracS / nStr : 0, nStr, nBnd ? fracB / nBnd : 0, nBnd);
   printf("highest populated amplitude bin: all %.0f ADC\n", hAll->GetBinCenter(hAll->FindLastBinAbove(0)));
   auto *c = new TCanvas("c", "c", 1400, 600);
   c->Divide(2, 1);
   c->cd(1); gPad->SetLogy(); hAll->SetLineColor(kBlack); hAll->Draw("hist");
   c->cd(2); gPad->SetLogy();
   hBnd->Scale(1.0 / hBnd->Integral()); hStr->Scale(1.0 / hStr->Integral());
   hBnd->SetLineColor(kBlue + 1); hStr->SetLineColor(kRed + 1); hBnd->SetLineWidth(2); hStr->SetLineWidth(2);
   hBnd->SetTitle("hits on band (blue) and streak (red) tracks, Brho 0.11-0.16;amplitude [ADC];fraction");
   hBnd->SetMaximum(1.5 * std::max(hBnd->GetMaximum(), hStr->GetMaximum()));
   hBnd->Draw("hist"); hStr->Draw("hist same");
   c->SaveAs(png);
}
