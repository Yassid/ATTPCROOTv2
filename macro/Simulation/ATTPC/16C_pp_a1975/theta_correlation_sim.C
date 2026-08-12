/// @file theta_correlation_sim.C
/// @brief Generated vs reconstructed polar angle: the correlation, the bias and the resolution.
///
/// WHY THIS IS WORTH CHECKING SEPARATELY. Every efficiency quoted from this simulation is built on
/// a truth match of |theta_rec - theta_true| < 10 deg, and a window that wide would absorb a
/// systematic offset without complaint -- an acceptance can look healthy while the angles it is
/// binned in are shifted. So the residual is shown UNCUT here, and the matched version alongside
/// it, and the bias is plotted against angle rather than quoted as one number.
///
/// The convention is checked rather than assumed: AtSpyralPID reports the polar angle such that
/// the lab angle is 180 - polar, and the raw one is drawn beside it so the two can be told apart
/// on sight instead of by argument.
///
///   root -b -q 'theta_correlation_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void theta_correlation_sim(TString dir = "/mnt/f/a1975_C16_pp_pid",
                           TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                           TString png = "plots/theta_correlation_sim.png", Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   auto *hC = new TH2D("hC", "#theta_{rec} = 180 - polar   vs   truth;#theta_{true} [deg];#theta_{rec} [deg]",
                       90, 0, 90, 180, 0, 180);
   auto *hRaw = new TH2D("hRaw", "raw polar   vs   truth;#theta_{true} [deg];polar [deg]", 90, 0, 90, 180, 0, 180);
   auto *hRes = new TH1D("hRes", "#theta_{rec} - #theta_{true}, ALL valid entries;#Delta#theta [deg];entries",
                         200, -100, 100);
   auto *hResM = new TH1D("hResM", "#Delta#theta for the matched track;#Delta#theta [deg];entries", 100, -10, 10);
   auto *hRvT = new TH2D("hRvT", "#Delta#theta vs #theta_{true};#theta_{true} [deg];#Delta#theta [deg]",
                         45, 0, 90, 200, -20, 20);
   long nAll = 0, nMatch = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fsim = dir + "/" + tg + "_sim.root", fpid = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(fpid))
         continue;
      TFile *fs = TFile::Open(fsim), *fp = TFile::Open(fpid);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tp = fp ? (TTree *)fp->Get("cbmsim") : nullptr;
      if (!ts || !tp || tp->GetEntries() > ts->GetEntries()) {
         if (fs) fs->Close();
         if (fp) fp->Close();
         continue;
      }
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tp->SetBranchAddress("AtPIDEvent", &pe);

      for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
         ts->GetEntry(i);
         tp->GetEntry(i);
         double thTrue = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            thTrue = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thTrue < 0 || !pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         double bd = 1e9, bres = 0;
         bool got = false;
         for (auto &sp : ev->GetSpyral()) {
            if (!sp.valid)
               continue;
            double raw = sp.polar * TMath::RadToDeg();
            double rec = 180.0 - raw;
            ++nAll;
            hC->Fill(thTrue, rec);
            hRaw->Fill(thTrue, raw);
            hRes->Fill(rec - thTrue);
            if (std::fabs(rec - thTrue) < bd) {
               bd = std::fabs(rec - thTrue);
               bres = rec - thTrue;
               got = true;
            }
         }
         if (got && bd < dThetaMax) {
            ++nMatch;
            hResM->Fill(bres);
            hRvT->Fill(thTrue, bres);
         }
      }
      fs->Close();
      fp->Close();
   }
   printf("\n  %ld valid PID entries paired with a truth proton, %ld matched within %.0f deg\n", nAll, nMatch,
          dThetaMax);
   printf("  ALL entries:  mean %+.2f deg, rms %.2f deg\n", hRes->GetMean(), hRes->GetRMS());
   printf("  MATCHED    :  mean %+.2f deg, rms %.2f deg\n\n", hResM->GetMean(), hResM->GetRMS());

   // bias and resolution against angle -- a single mean would hide a rotation
   printf("  theta_true    n      bias [deg]   rms [deg]\n");
   auto *hBias = new TH1D("hBias", "bias;#theta_{true} [deg];<#Delta#theta> [deg]", 45, 0, 90);
   auto *hSig = new TH1D("hSig", "resolution;#theta_{true} [deg];rms(#Delta#theta) [deg]", 45, 0, 90);
   for (int b = 1; b <= hRvT->GetNbinsX(); ++b) {
      TH1D *p = hRvT->ProjectionY("_py", b, b);
      if (p->GetEntries() < 30) {
         delete p;
         continue;
      }
      double m = p->GetMean(), s = p->GetRMS();
      hBias->SetBinContent(b, m);
      hBias->SetBinError(b, s / std::sqrt(p->GetEntries()));
      hSig->SetBinContent(b, s);
      if (b % 3 == 1)
         printf("  %5.1f      %6.0f     %+6.2f      %5.2f\n", hRvT->GetXaxis()->GetBinCenter(b), p->GetEntries(), m, s);
      delete p;
   }

   TCanvas *c = new TCanvas("cT", "theta correlation", 1500, 900);
   c->Divide(3, 2);
   c->cd(1); gPad->SetLogz(); hC->Draw("colz");
   { auto *l = new TLine(0, 0, 90, 90); l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw(); }
   c->cd(2); gPad->SetLogz(); hRaw->Draw("colz");
   { auto *l = new TLine(0, 0, 90, 90); l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw(); }
   c->cd(3); gPad->SetLogy(); hRes->Draw("hist");
   c->cd(4); hResM->Draw("hist");
   c->cd(5);
   hBias->SetLineColor(kRed + 1); hBias->SetLineWidth(2); hBias->SetMinimum(-6); hBias->SetMaximum(6);
   hBias->Draw("e");
   { auto *l = new TLine(0, 0, 90, 0); l->SetLineStyle(2); l->Draw(); }
   c->cd(6);
   hSig->SetLineColor(kBlue + 1); hSig->SetLineWidth(2); hSig->SetMinimum(0);
   hSig->Draw("hist");

   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s\n", png.Data());
}
