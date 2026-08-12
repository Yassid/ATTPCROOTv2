/// @file fit_vs_truth_sim.C
/// @brief Fitted track vs MC truth: angle and energy, bias and resolution, against the PID estimator.
///
/// The Spyral PID angle carries a systematic offset -- flat at about -0.8 deg forward and growing
/// to -3 deg by 80 deg (theta_correlation_sim.C). That estimator is NOT what the physics is binned
/// in: the angular distributions come from the fitted track. So the question this answers is
/// whether the fit inherits that bias, cancels it, or has one of its own, and the two are computed
/// in the SAME event loop so the comparison carries no bookkeeping difference.
///
/// Units are checked rather than assumed. EventFit::AtGenfitter reports kinematics.theta in
/// RADIANS; the legacy AtGenfit path reported degrees, and the two are told apart here by the
/// observed range instead of by which macro was copied from.
///
/// B FIELD. The fit uses -2.85, the SAME sign as the data, which contradicts the note in
/// run_reco_C16_TC.C that simulation needs +1. Measured on 800 events: +2.85 gives dKE -3.53 +-
/// 5.80 MeV and returns polar angles up to 147 deg, which a recoil proton cannot have; -2.85 gives
/// -0.27 +- 2.70 MeV and stops at 84 deg, matching the generated range. The note is wrong for this
/// chain, or at least for this fitter.
///
///   root -b -q 'fit_vs_truth_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void fit_vs_truth_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                      TString fitSuffix = "_genfitter_pp", TString png = "plots/fit_vs_truth_sim.png",
                      Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double mp = 938.272;

   auto *hTh = new TH2D("hTh", "fitted vs true angle;#theta_{true} [deg];#theta_{fit} [deg]", 90, 0, 90, 90, 0, 90);
   auto *hKE = new TH2D("hKE", "fitted vs true KE;KE_{true} [MeV];KE_{fit} [MeV]", 90, 0, 45, 90, 0, 45);
   auto *hRvT = new TH2D("hRvT", ";#theta_{true} [deg];#Delta#theta [deg]", 45, 0, 90, 200, -20, 20);
   auto *hKvK = new TH2D("hKvK", ";KE_{true} [MeV];#DeltaKE [MeV]", 45, 0, 45, 200, -20, 20);
   auto *hSvT = new TH2D("hSvT", ";#theta_{true} [deg];#Delta#theta [deg]", 45, 0, 90, 200, -20, 20); // Spyral
   auto *genTh = new TH1D("genTh", "", 45, 0, 90);
   auto *fitTh = new TH1D("fitTh", "", 45, 0, 90);
   auto *spyTh = new TH1D("spyTh", "", 45, 0, 90);
   long nGen = 0, nFit = 0, nSpy = 0;
   double rawMin = 1e9, rawMax = -1e9;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fsim = dir + "/" + tg + "_sim.root", ffit = dir + "/" + tg + fitSuffix + ".root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(ffit)) {
         printf("  skip %s (missing)\n", tg.Data());
         continue;
      }
      TFile *fs = TFile::Open(fsim), *ff = TFile::Open(ffit);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tf = ff ? (TTree *)ff->Get("cbmsim") : nullptr;
      if (!ts || !tf || tf->GetEntries() > ts->GetEntries()) {
         if (fs) fs->Close();
         if (ff) ff->Close();
         continue;
      }
      TClonesArray *mc = nullptr, *te = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);
      tf->SetBranchAddress("AtPIDEvent", &pe);

      for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
         ts->GetEntry(i);
         tf->GetEntry(i);
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            keT = std::sqrt(pp * pp + mp * mp) - mp;
            thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (keT < 0)
            continue;
         ++nGen;
         genTh->Fill(thT);

         // ---- fitted track, best match in angle ----
         double bd = 1e9, bTh = 0, bKE = 0;
         bool got = false;
         if (te && te->GetEntriesFast()) {
            auto *ev = (AtTrackingEvent *)te->At(0);
            if (ev)
               for (auto &ft : ev->GetFittedTracks()) {
                  if (!ft)
                     continue;
                  auto &k = ft->GetKinematics();
                  rawMin = std::min(rawMin, (double)k.theta);
                  rawMax = std::max(rawMax, (double)k.theta);
                  double th = k.theta * TMath::RadToDeg(); // radians for EventFit::AtGenfitter
                  if (std::fabs(th - thT) < bd) {
                     bd = std::fabs(th - thT);
                     bTh = th;
                     bKE = k.kineticEnergy;
                     got = true;
                  }
               }
         }
         if (got && bd < dThetaMax) {
            ++nFit;
            fitTh->Fill(thT);
            hTh->Fill(thT, bTh);
            hKE->Fill(keT, bKE);
            hRvT->Fill(thT, bTh - thT);
            hKvK->Fill(keT, bKE - keT);
         }

         // ---- Spyral PID angle, same event, for the side-by-side ----
         double sd = 1e9, sTh = 0;
         bool gotS = false;
         if (pe && pe->GetEntriesFast()) {
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (ev)
               for (auto &sp : ev->GetSpyral()) {
                  if (!sp.valid)
                     continue;
                  double th = 180.0 - sp.polar * TMath::RadToDeg();
                  if (std::fabs(th - thT) < sd) {
                     sd = std::fabs(th - thT);
                     sTh = th;
                     gotS = true;
                  }
               }
         }
         if (gotS && sd < dThetaMax) {
            ++nSpy;
            spyTh->Fill(thT);
            hSvT->Fill(thT, sTh - thT);
         }
      }
      fs->Close();
      ff->Close();
      printf("  %-8s done\n", tg.Data());
   }
   printf("\n  raw kinematics.theta range %.3f .. %.3f  => %s\n", rawMin, rawMax,
          rawMax > 6.5 ? "DEGREES (units assumption WRONG, fix the macro)" : "radians, as assumed");
   printf("  %ld generated,  %ld matched to a FITTED track (%.1f %%),  %ld to a Spyral entry (%.1f %%)\n\n", nGen, nFit,
          nGen ? 100.0 * nFit / nGen : 0., nSpy, nGen ? 100.0 * nSpy / nGen : 0.);

   auto band = [](TH2D *h, TH1D *&bias, TH1D *&sig, const char *nm) {
      bias = new TH1D(Form("%s_b", nm), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
      sig = new TH1D(Form("%s_s", nm), "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         TH1D *p = h->ProjectionY("_p", b, b);
         if (p->GetEntries() >= 30) {
            bias->SetBinContent(b, p->GetMean());
            bias->SetBinError(b, p->GetRMS() / std::sqrt(p->GetEntries()));
            sig->SetBinContent(b, p->GetRMS());
         }
         delete p;
      }
   };
   TH1D *bFit, *sFit, *bSpy, *sSpy, *bKE, *sKE;
   band(hRvT, bFit, sFit, "fit");
   band(hSvT, bSpy, sSpy, "spy");
   band(hKvK, bKE, sKE, "ke");

   printf("  theta_true   fit bias   fit rms   |   Spyral bias   Spyral rms\n");
   for (int b = 1; b <= bFit->GetNbinsX(); b += 3)
      if (sFit->GetBinContent(b) > 0)
         printf("  %6.1f      %+6.2f     %5.2f    |    %+6.2f        %5.2f\n", bFit->GetBinCenter(b),
                bFit->GetBinContent(b), sFit->GetBinContent(b), bSpy->GetBinContent(b), sSpy->GetBinContent(b));

   TCanvas *c = new TCanvas("cF", "fit vs truth", 1500, 900);
   c->Divide(3, 2);
   c->cd(1); gPad->SetLogz(); hTh->Draw("colz");
   { auto *l = new TLine(0, 0, 90, 90); l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw(); }
   c->cd(2); gPad->SetLogz(); hKE->Draw("colz");
   { auto *l = new TLine(0, 0, 45, 45); l->SetLineColor(kRed); l->SetLineWidth(2); l->Draw(); }
   c->cd(3);
   bFit->SetTitle("angle bias: fit (red) vs Spyral PID (blue);#theta_{true} [deg];<#Delta#theta> [deg]");
   bFit->SetLineColor(kRed + 1); bFit->SetLineWidth(2); bFit->SetMinimum(-5); bFit->SetMaximum(5);
   bSpy->SetLineColor(kAzure + 2); bSpy->SetLineWidth(2);
   bFit->Draw("e"); bSpy->Draw("e same");
   { auto *l = new TLine(0, 0, 90, 0); l->SetLineStyle(2); l->Draw(); }
   c->cd(4);
   sFit->SetTitle("angle resolution: fit (red) vs Spyral (blue);#theta_{true} [deg];rms [deg]");
   sFit->SetLineColor(kRed + 1); sFit->SetLineWidth(2); sFit->SetMinimum(0);
   sSpy->SetLineColor(kAzure + 2); sSpy->SetLineWidth(2);
   sFit->Draw("hist"); sSpy->Draw("hist same");
   c->cd(5);
   bKE->SetTitle("KE bias;KE_{true} [MeV];<#DeltaKE> [MeV]");
   bKE->SetLineColor(kRed + 1); bKE->SetLineWidth(2); bKE->SetMinimum(-6); bKE->SetMaximum(6);
   bKE->Draw("e");
   { auto *l = new TLine(0, 0, 45, 0); l->SetLineStyle(2); l->Draw(); }
   c->cd(6);
   auto *eF = (TH1D *)fitTh->Clone("eF"); eF->Divide(fitTh, genTh, 1, 1, "B");
   auto *eS = (TH1D *)spyTh->Clone("eS"); eS->Divide(spyTh, genTh, 1, 1, "B");
   eF->SetTitle("efficiency: fit (red) vs Spyral (blue);#theta_{true} [deg];matched / generated");
   eF->SetLineColor(kRed + 1); eF->SetLineWidth(2); eF->SetMinimum(0); eF->SetMaximum(1.05);
   eS->SetLineColor(kAzure + 2); eS->SetLineWidth(2);
   eF->Draw("hist"); eS->Draw("hist same");

   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s\n", png.Data());
}
