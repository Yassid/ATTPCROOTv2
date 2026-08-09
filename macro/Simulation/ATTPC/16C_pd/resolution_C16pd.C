/// @file resolution_C16pd.C
/// @brief What the reconstruction does to the deuteron, measured against MC truth.
///
/// Pairs the fitted tracks with the truth by entry index, so the fit MUST be of the ungated
/// reconstruction of the SAME simulation file; the macro refuses to run if the entry counts
/// disagree rather than pairing the wrong events together.
///
/// THE POINT OF THE PLOT is that a single resolution number is misleading here. The deuteron
/// energy resolution degrades by a factor ten across the energy range the ground state covers,
/// so an excitation-energy resolution quoted for the whole angular range is an average over
/// conditions that differ enormously. The bottom row separates it by angle, which is the form
/// the analysis can actually use: there is a region of large lab angle -- low deuteron energy --
/// where the ground state and the 0.740 MeV state should be separable, and a forward region where
/// they cannot be.
///
/// The cause is geometric. A high-energy deuteron has a large helix radius and its arc does not
/// close inside the chamber, so the curvature is poorly constrained. At 52 MeV the radius is about
/// 30 cm against a 25 cm chamber, which is the same failure that produced the energy bias in the
/// a1954 analysis above theta_lab 62 deg.
///
/// Resolutions are quoted as the CORE sigma of a gaussian fitted over a restricted range, not the
/// rms: the residual distributions have tails from failed fits that make the rms four times the
/// core width and describe neither.
///
///   root -b -q 'resolution_C16pd.C()'

void resolution_C16pd(TString simFile = "/mnt/f/a1975_C16_pd_sim/gs_s1001_sim.root",
                      TString fitFile = "/mnt/f/a1975_C16_pd_sim/gs_s1001_genfitter.root", Double_t resEx = 0.0,
                      Double_t Ebeam = 192.0, Double_t dEdrift = 28.0, Double_t chi2Cut = 10.0, TString tag = "gs")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const double u = 931.49401;
   const double mb = 16.0147 * u, mt = 1.0078250322 * u, md = 2.0141017781 * u, mr = 15.0105993 * u;
   const double mres = mr + resEx;

   TFile *fs = TFile::Open(simFile);
   TFile *ff = TFile::Open(fitFile);
   if (!fs || fs->IsZombie() || !ff || ff->IsZombie()) {
      printf("\033[1;31mcannot open %s or %s\033[0m\n", simFile.Data(), fitFile.Data());
      return;
   }
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TTree *tf = (TTree *)ff->Get("cbmsim");
   if (!ts || !tf)
      return;
   if (ts->GetEntries() != tf->GetEntries()) {
      printf("\033[1;31mENTRY MISMATCH: sim %lld vs fit %lld -- pairing is by index, so this would\n"
             "compare unrelated events. Fit the UNGATED reco of this same sim file.\033[0m\n",
             ts->GetEntries(), tf->GetEntries());
      return;
   }

   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   auto exOf = [&](double Eb, double ke, double th) {
      double E1 = Eb + mb, p1 = std::sqrt(Eb * Eb + 2 * Eb * mb);
      double E3 = ke + md, p3 = std::sqrt(E3 * E3 - md * md);
      double E4 = E1 + mt - E3, p4 = p1 * p1 + p3 * p3 - 2 * p1 * p3 * std::cos(th);
      return std::sqrt(std::max(E4 * E4 - p4, 0.0)) - mres;
   };

   auto *hcorr = new TH2D("hcorr", "deuteron energy: reconstructed vs true;KE_{true} [MeV];KE_{fit} [MeV]", 60, 0, 60,
                          60, 0, 60);
   auto *hdke = new TH1D("hdke", "KE_{fit} - KE_{true};#DeltaKE [MeV];tracks", 160, -16, 16);
   auto *k2 = new TH2D("k2", "", 12, 0, 60, 200, -16, 16);
   auto *t2 = new TH2D("t2", "", 10, 0, 40, 200, -16, 16);
   auto *tex = new TH2D("tex", "", 10, 0, 40, 200, -4, 8);
   auto *exFwd = new TH1D("exFwd", "", 120, -4, 8);
   auto *exBwd = new TH1D("exBwd", "", 120, -4, 8);

   long nfit = 0, ngood = 0;
   for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
      ts->GetEntry(i);
      tf->GetEntry(i);
      double keT = -1, thT = -1, zT = 0;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
            continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pp = std::sqrt(px * px + py * py + pz * pz);
         if (pp <= 0)
            continue;
         keT = std::sqrt(pp * pp + md * md) - md;
         thT = std::acos(pz / pp);
         zT = p->GetStartZ() * 10;
         break;
      }
      if (keT < 0 || !te || te->GetEntriesFast() == 0)
         continue;
      auto *e = (AtTrackingEvent *)te->At(0);
      if (!e)
         continue;
      for (auto &ft : e->GetFittedTracks()) {
         if (!ft)
            continue;
         ++nfit;
         const auto &m = ft->GetTrackMetadata();
         double ndf = m ? m->GetNdf() : 0, c2 = m ? m->GetChi2() : 0;
         if (!(ndf > 0 && c2 / ndf < chi2Cut))
            break;
         ++ngood;
         double ke = ft->GetKinematics().kineticEnergy, th = ft->GetKinematics().theta;
         double zR = ft->GetVertex(0).Z(), thDeg = th * TMath::RadToDeg();
         hcorr->Fill(keT, ke);
         hdke->Fill(ke - keT);
         k2->Fill(keT, ke - keT);
         t2->Fill(thDeg, ke - keT);
         double ex = exOf(Ebeam - dEdrift * zR / 1000.0, ke, th);
         tex->Fill(thDeg, ex);
         if (thDeg > 30)
            exBwd->Fill(ex);
         else
            exFwd->Fill(ex);
         break;
      }
   }
   printf("\n  %ld fitted tracks, %ld with chi2/ndf < %.0f\n", nfit, ngood, chi2Cut);

   // core sigma of a slice, or -1 if the slice is too thin to fit
   auto coreOf = [](TH2D *h, int b, double lo, double hi, double &bias) {
      TH1D *s = h->ProjectionY("s_tmp", b, b);
      if (s->Integral() < 40) {
         delete s;
         bias = 0;
         return -1.0;
      }
      TF1 g("g", "gaus", lo, hi);
      s->Fit(&g, "QNR");
      bias = g.GetParameter(1);
      double sg = g.GetParameter(2);
      delete s;
      return sg;
   };

   auto *gKE = new TGraph();
   auto *gTH = new TGraph();
   auto *gEX = new TGraph();
   int n1 = 0, n2 = 0, n3 = 0;
   double bias;
   printf("\n   KE[MeV]  sigma(KE)      theta[deg]  sigma(KE)   sigma(Ex)\n");
   for (int b = 1; b <= 12; ++b) {
      double s = coreOf(k2, b, -6, 6, bias);
      if (s > 0)
         gKE->SetPoint(n1++, k2->GetXaxis()->GetBinCenter(b), s);
   }
   for (int b = 1; b <= 10; ++b) {
      double s = coreOf(t2, b, -6, 6, bias);
      double sx = coreOf(tex, b, -3, 5, bias);
      double c = t2->GetXaxis()->GetBinCenter(b);
      if (s > 0)
         gTH->SetPoint(n2++, c, s);
      if (sx > 0)
         gEX->SetPoint(n3++, c, sx);
      if (s > 0)
         printf("   %5.0f    %6.2f          %5.0f      %6.2f      %6.2f\n", k2->GetXaxis()->GetBinCenter(b),
                n1 > b - 1 ? gKE->GetPointY(std::min(b - 1, n1 - 1)) : 0, c, s, sx > 0 ? sx : 0);
   }

   TCanvas *c1 = new TCanvas("cres", "resolution", 1500, 950);
   c1->Divide(2, 2);

   c1->cd(1);
   gPad->SetLogz();
   hcorr->Draw("colz");
   auto *diag = new TLine(0, 0, 60, 60);
   diag->SetLineColor(kRed + 1);
   diag->SetLineWidth(2);
   diag->SetLineStyle(2);
   diag->Draw();

   c1->cd(2);
   auto *fr2 = new TH1D("fr2", "energy resolution;KE_{true} [MeV]   /   #theta_{lab} [deg];core #sigma(KE) [MeV]", 1,
                        0, 60);
   fr2->SetMinimum(0);
   fr2->SetMaximum(5);
   fr2->Draw();
   gKE->SetMarkerStyle(20);
   gKE->SetMarkerColor(kAzure + 2);
   gKE->SetLineColor(kAzure + 2);
   gKE->SetLineWidth(3);
   gKE->SetMarkerSize(1.3);
   gKE->Draw("LP");
   gTH->SetMarkerStyle(21);
   gTH->SetMarkerColor(kRed + 1);
   gTH->SetLineColor(kRed + 1);
   gTH->SetLineWidth(3);
   gTH->SetMarkerSize(1.3);
   gTH->Draw("LP");
   auto *lg = new TLegend(0.16, 0.70, 0.60, 0.87);
   lg->AddEntry(gKE, "vs true KE [MeV]", "lp");
   lg->AddEntry(gTH, "vs #theta_{lab} [deg]", "lp");
   lg->SetTextSize(0.04);
   lg->Draw();

   c1->cd(3);
   double mx = std::max(exFwd->GetMaximum(), exBwd->GetMaximum());
   auto *fr3 = new TH1D("fr3", TString::Format("E_{x} (generated %.2f MeV), split by lab angle;E_{x} [MeV];tracks",
                                               resEx),
                        1, -4, 8);
   fr3->SetMaximum(mx * 1.3);
   fr3->Draw();
   exBwd->SetLineColor(kRed + 1);
   exBwd->SetLineWidth(2);
   exBwd->Draw("hist same");
   exFwd->SetLineColor(kAzure + 2);
   exFwd->SetLineWidth(2);
   exFwd->Draw("hist same");
   auto *lv = new TLine(resEx, 0, resEx, mx * 1.3);
   lv->SetLineStyle(2);
   lv->SetLineWidth(2);
   lv->Draw();
   auto *lg3 = new TLegend(0.55, 0.68, 0.89, 0.87);
   lg3->AddEntry(exBwd, "#theta_{lab} > 30#circ", "l");
   lg3->AddEntry(exFwd, "#theta_{lab} < 30#circ", "l");
   lg3->SetTextSize(0.04);
   lg3->Draw();

   c1->cd(4);
   auto *fr4 = new TH1D("fr4", "excitation-energy resolution vs angle;#theta_{lab} [deg];core #sigma(E_{x}) [MeV]", 1,
                        0, 40);
   fr4->SetMinimum(0);
   fr4->SetMaximum(3.5);
   fr4->Draw();
   gEX->SetMarkerStyle(20);
   gEX->SetMarkerColor(kGreen + 3);
   gEX->SetLineColor(kGreen + 3);
   gEX->SetLineWidth(3);
   gEX->SetMarkerSize(1.4);
   gEX->Draw("LP");
   // the level spacing this has to beat to separate the first two states
   auto *lsp = new TLine(0, 0.74, 40, 0.74);
   lsp->SetLineColor(kRed + 1);
   lsp->SetLineWidth(3);
   lsp->SetLineStyle(2);
   lsp->Draw();
   auto *tx = new TLatex(2, 0.85, "0.740 MeV level spacing");
   tx->SetTextColor(kRed + 1);
   tx->SetTextSize(0.038);
   tx->Draw();

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/resolution_" + tag + ".png";
   c1->SaveAs(png);
   printf("\n  wrote %s\n\n", png.Data());
}
