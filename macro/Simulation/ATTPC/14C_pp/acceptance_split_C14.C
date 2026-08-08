/// @file acceptance_split_C14.C
/// @brief Same acceptance as acceptance_C14.C, but split by WHICH cut removed the event.
///
/// acceptance_C14.C answers "how many generated reactions become a good fitted proton". When the
/// answer has a hole in it -- the GENFIT g.s. acceptance drops to 0.42 at theta_cm 25-35 deg while
/// the UKF one is 0.68-0.91 there -- the number alone cannot say whether the tracks are missing,
/// unconverged, or merely failing the truth match. This macro keeps the identical denominator and
/// numerator logic and additionally classifies every generated reaction by how far it got:
///
///   NOTRK    no fitted track in the event at all
///   KEBAD    a fit exists but its kinetic energy is not in (0, 1000) MeV
///   CHI2     survives that, fails chi2/ndf < chi2Cut
///   DTH      survives chi2, fails |theta_fit - theta_true| < dThetaMax
///   RATIO    survives chi2 and dtheta, fails KE_fit/KE_true in [keRatioMin, keRatioMax]
///   BOTH     survives chi2, fails dtheta AND ratio
///   GOOD     passes everything (this is acceptance_C14.C's numerator)
///
/// An event is labelled by the BEST outcome over its fitted tracks, so the categories partition
/// the denominator exactly and the fractions in each theta_cm bin sum to 1.
///
/// It also builds the counterfactual acceptances -- what the curve would be with one cut removed:
/// noDth, noRatio, noMatch (neither), noChi2 (no chi2 cut, match kept), any (any converged fit,
/// which is what a numerator with no truth match at all would give). Comparing those to the
/// nominal curve says directly which cut owns the hole.
///
///   root -b -q 'acceptance_split_C14.C("/mnt/f/a1954_C14_acc_gf/gs_s1001_sim.root",
///                                      "/mnt/f/a1954_C14_acc_gf/gs_s1001_genfit.root",
///                                      "gs",0.0,161.0,5.0,36,180.0,10.0,0.5,2.0,kTRUE,"gf_s1001")'

#include <tuple>

static double spl_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics, identical to acceptance_C14.C / pp/ex_C14.C
static std::tuple<double, double> spl_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                           double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * spl_omega2(s, m1 * m1, m2 * m2) * spl_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (spl_omega2(s, m1 * m1, m2 * m2) * spl_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

void acceptance_split_C14(TString simFile, TString fitFile, TString tag = "gs", Double_t resEx = 0.0,
                          Double_t Ebeam = 161.0, Double_t chi2Cut = 5.0, Int_t nBins = 36, Double_t cmMax = 180.0,
                          Double_t dThetaMax = 10.0, Double_t keRatioMin = 0.5, Double_t keRatioMax = 2.0,
                          Bool_t useXtr = kFALSE, TString outTag = "", TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;
   const double m_resid = m_C14 + resEx;

   TFile *fs = TFile::Open(simFile);
   TFile *ff = TFile::Open(fitFile);
   if (!fs || fs->IsZombie() || !ff || ff->IsZombie()) {
      printf("\033[1;31mcannot open %s or %s\033[0m\n", simFile.Data(), fitFile.Data());
      return;
   }
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TTree *tf = (TTree *)ff->Get("cbmsim");
   if (!ts || !tf) {
      printf("\033[1;31mmissing cbmsim\033[0m\n");
      return;
   }
   if (ts->GetEntries() != tf->GetEntries()) {
      printf("\033[1;31mENTRY MISMATCH: sim %lld vs fit %lld\033[0m\n", ts->GetEntries(), tf->GetEntries());
      return;
   }

   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   // category codes, best-outcome-wins (higher = got further)
   enum { CAT_NOTRK = 0, CAT_KEBAD, CAT_CHI2, CAT_BOTH, CAT_DTH, CAT_RATIO, CAT_GOOD, NCAT };
   const char *catName[NCAT] = {"NOTRK", "KEBAD", "CHI2", "BOTH", "DTH", "RATIO", "GOOD"};
   const int catColor[NCAT] = {kGray + 1, kGray + 2, kOrange + 7, kMagenta + 1, kRed + 1, kGreen + 2, kAzure + 2};

   TH1D *hGen = new TH1D("hGen_" + tag, "generated;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
   hGen->Sumw2();
   TH1D *hCat[NCAT];
   for (int c = 0; c < NCAT; ++c) {
      hCat[c] = new TH1D(TString::Format("hCat_%s_%s", tag.Data(), catName[c]), catName[c], nBins, 0, cmMax);
      hCat[c]->Sumw2();
   }
   // counterfactual numerators: same loop, one cut dropped
   TH1D *hNoDth = new TH1D("hNoDth_" + tag, "", nBins, 0, cmMax);
   TH1D *hNoRat = new TH1D("hNoRat_" + tag, "", nBins, 0, cmMax);
   TH1D *hNoMatch = new TH1D("hNoMatch_" + tag, "", nBins, 0, cmMax);
   TH1D *hNoChi2 = new TH1D("hNoChi2_" + tag, "", nBins, 0, cmMax);
   TH1D *hAny = new TH1D("hAny_" + tag, "", nBins, 0, cmMax);
   for (auto *h : {hNoDth, hNoRat, hNoMatch, hNoChi2, hAny})
      h->Sumw2();

   // how badly the survivors miss, in the dip region and on the plateau
   TH1D *hDth25 = new TH1D("hDth25_" + tag, "|#Delta#theta| of chi2-passing fits;|#theta_{fit}-#theta_{true}| [deg];fits",
                           90, 0, 180);
   TH1D *hDth45 = new TH1D("hDth45_" + tag, "", 90, 0, 180);
   TH1D *hRat25 = new TH1D("hRat25_" + tag, "KE_{fit}/KE_{true} of chi2-passing fits;ratio;fits", 100, 0, 5);
   TH1D *hRat45 = new TH1D("hRat45_" + tag, "", 100, 0, 5);
   TH2D *hKE = new TH2D("hKE_" + tag, "chi2-passing fits;KE_{true} [MeV];KE_{fit} [MeV]", 100, 0, 40, 100, 0, 40);
   TH1D *hNtrk25 = new TH1D("hNtrk25_" + tag, "fitted tracks per reaction;n tracks;events", 8, -0.5, 7.5);
   TH1D *hNtrk45 = new TH1D("hNtrk45_" + tag, "", 8, -0.5, 7.5);
   TH1D *hKEtrue = new TH1D("hKEtrue_" + tag, "true recoil KE;KE_{true} [MeV];reactions", 100, 0, 40);

   Long64_t N = ts->GetEntries();
   long nGen = 0;
   long catTot[NCAT] = {0};
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc)
         continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *t = (AtMCTrack *)mc->At(k);
         if (!t || t->GetPdgCode() != 2212 || t->GetMotherId() != -1)
            continue;
         double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000;
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0)
            continue;
         keT = std::sqrt(p * p + m_p * m_p) - m_p;
         thT = std::acos(pz / p);
         break;
      }
      if (keT <= 0)
         continue;
      auto [exT, cmT] = spl_kine(m_C14, m_p, m_p, m_resid, Ebeam, thT, keT);
      if (std::isnan(cmT))
         continue;
      ++nGen;
      hGen->Fill(cmT);
      hKEtrue->Fill(keT);
      const bool dip = (cmT >= 25 && cmT < 35);      // the GENFIT hole
      const bool plateau = (cmT >= 45 && cmT < 55);  // where both fitters are fine

      tf->GetEntry(i);
      int best = CAT_NOTRK;
      bool anyFit = false, okNoDth = false, okNoRat = false, okNoMatch = false, okNoChi2 = false;
      int nTrk = 0;
      if (te && te->GetEntriesFast() > 0) {
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (ev)
            for (auto &ft : ev->GetFittedTracks()) {
               if (!ft)
                  continue;
               ++nTrk;
               const auto &md = ft->GetTrackMetadata();
               double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
               double c2n = ndf > 0 ? chi2 / ndf : 1e9;
               const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
               double ke = kin.kineticEnergy;
               if (!(ke > 0 && ke < 1000)) {
                  if (best < CAT_KEBAD)
                     best = CAT_KEBAD;
                  continue;
               }
               double dth = std::fabs(kin.theta * TMath::RadToDeg() - thT * TMath::RadToDeg());
               double r = ke / keT;
               bool dthOK = (dThetaMax <= 0) || dth <= dThetaMax;
               bool ratOK = r >= keRatioMin && r <= keRatioMax;
               anyFit = true;
               if (dthOK && ratOK)
                  okNoChi2 = true;
               if (c2n >= chi2Cut) {
                  if (best < CAT_CHI2)
                     best = CAT_CHI2;
                  continue;
               }
               // survives chi2 -- record how far off it is
               if (dip) {
                  hDth25->Fill(dth);
                  hRat25->Fill(r);
               } else if (plateau) {
                  hDth45->Fill(dth);
                  hRat45->Fill(r);
               }
               hKE->Fill(keT, ke);
               okNoMatch = true;
               if (ratOK)
                  okNoDth = true;
               if (dthOK)
                  okNoRat = true;
               int cat = (dthOK && ratOK) ? CAT_GOOD : (dthOK ? CAT_RATIO : (ratOK ? CAT_DTH : CAT_BOTH));
               if (cat > best)
                  best = cat;
            }
      }
      if (dip)
         hNtrk25->Fill(std::min(nTrk, 7));
      else if (plateau)
         hNtrk45->Fill(std::min(nTrk, 7));
      hCat[best]->Fill(cmT);
      ++catTot[best];
      if (okNoDth)
         hNoDth->Fill(cmT);
      if (okNoRat)
         hNoRat->Fill(cmT);
      if (okNoMatch)
         hNoMatch->Fill(cmT);
      if (okNoChi2)
         hNoChi2->Fill(cmT);
      if (anyFit)
         hAny->Fill(cmT);
   }

   printf("\n===== %s : %lld entries, %ld generated reactions =====\n", fitFile.Data(), N, nGen);
   printf("cuts: chi2/ndf < %g, |dtheta| <= %g deg, KE ratio in [%g, %g], useXtr = %d\n", chi2Cut, dThetaMax,
          keRatioMin, keRatioMax, (int)useXtr);
   printf("\noverall outcome of every generated reaction:\n");
   for (int c = NCAT - 1; c >= 0; --c)
      printf("   %-6s %7ld  %5.1f %%\n", catName[c], catTot[c], nGen ? 100.0 * catTot[c] / nGen : 0.0);

   auto acc = [&](TH1D *num, const char *nm) {
      auto *h = (TH1D *)num->Clone(TString::Format("hAcc_%s_%s", tag.Data(), nm));
      h->Divide(num, hGen, 1, 1, "B");
      return h;
   };
   TH1D *aGood = acc(hCat[CAT_GOOD], "nominal");
   TH1D *aNoDth = acc(hNoDth, "noDth");
   TH1D *aNoRat = acc(hNoRat, "noRatio");
   TH1D *aNoMatch = acc(hNoMatch, "noMatch");
   TH1D *aNoChi2 = acc(hNoChi2, "noChi2");
   TH1D *aAny = acc(hAny, "anyFit");

   printf("\n  theta_cm     gen |");
   for (int c = NCAT - 1; c >= 0; --c)
      printf(" %6s", catName[c]);
   printf(" |  nominal  noDth noRatio noMatch  noChi2  anyFit\n");
   for (int b = 1; b <= nBins; ++b) {
      double g = hGen->GetBinContent(b);
      double ctr = hGen->GetBinCenter(b);
      if (g < 1 || ctr > 155)
         continue;
      printf("  %3.0f-%3.0f %6.0f |", hGen->GetBinLowEdge(b), hGen->GetBinLowEdge(b) + hGen->GetBinWidth(b), g);
      for (int c = NCAT - 1; c >= 0; --c)
         printf(" %5.1f%%", 100.0 * hCat[c]->GetBinContent(b) / g);
      printf(" |  %6.3f %6.3f  %6.3f  %6.3f  %6.3f  %6.3f\n", aGood->GetBinContent(b), aNoDth->GetBinContent(b),
             aNoRat->GetBinContent(b), aNoMatch->GetBinContent(b), aNoChi2->GetBinContent(b),
             aAny->GetBinContent(b));
   }

   auto med = [](TH1D *h) {
      if (h->Integral() <= 0)
         return -1.0;
      double q = 0.5, x = 0;
      h->GetQuantiles(1, &x, &q);
      return x;
   };
   printf("\nchi2-passing fits, dip (theta_cm 25-35) vs plateau (45-55):\n");
   printf("   median |dtheta|   dip %6.2f deg   plateau %6.2f deg   (cut at %g)\n", med(hDth25), med(hDth45),
          dThetaMax);
   printf("   median KE ratio   dip %6.3f       plateau %6.3f       (cut [%g, %g])\n", med(hRat25), med(hRat45),
          keRatioMin, keRatioMax);
   printf("   mean fitted tracks/reaction   dip %.2f   plateau %.2f\n", hNtrk25->GetMean(), hNtrk45->GetMean());
   printf("   events with zero fitted track dip %5.1f %%   plateau %5.1f %%\n",
          hNtrk25->Integral() > 0 ? 100.0 * hNtrk25->GetBinContent(1) / hNtrk25->Integral() : 0.0,
          hNtrk45->Integral() > 0 ? 100.0 * hNtrk45->GetBinContent(1) / hNtrk45->Integral() : 0.0);

   // ---------- figure ----------
   TCanvas *c1 = new TCanvas("c1", "acceptance split", 1500, 950);
   c1->Divide(2, 2);

   c1->cd(1);
   auto *st = new THStack("st", TString::Format("where each generated reaction is lost (%s);#theta_{cm} [deg];fraction",
                                                outTag.Data()));
   for (int c = 0; c < NCAT; ++c) {
      auto *h = (TH1D *)hCat[c]->Clone(TString::Format("f_%s", catName[c]));
      h->Divide(h, hGen, 1, 1, "B");
      h->SetFillColor(catColor[c]);
      h->SetLineColor(catColor[c]);
      st->Add(h);
   }
   st->Draw("hist");
   st->GetXaxis()->SetRangeUser(15, 155);
   st->SetMaximum(1.05);
   auto *lg = new TLegend(0.62, 0.15, 0.89, 0.45);
   for (int c = NCAT - 1; c >= 0; --c) {
      auto *d = new TH1D(TString::Format("d_%s", catName[c]), "", 1, 0, 1);
      d->SetFillColor(catColor[c]);
      lg->AddEntry(d, catName[c], "f");
   }
   lg->Draw();

   c1->cd(2);
   auto style = [](TH1D *h, int col, int mk) {
      h->SetMarkerStyle(mk);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
   };
   style(aGood, kBlack, 20);
   style(aNoDth, kRed + 1, 21);
   style(aNoRat, kGreen + 2, 22);
   style(aNoMatch, kAzure + 2, 23);
   style(aNoChi2, kOrange + 7, 33);
   style(aAny, kMagenta + 1, 34);
   aGood->SetTitle("acceptance with one cut dropped;#theta_{cm} [deg];acceptance");
   aGood->GetXaxis()->SetRangeUser(15, 155);
   aGood->GetYaxis()->SetRangeUser(0, 1.15);
   aGood->Draw("E1");
   aNoDth->Draw("E1 same");
   aNoRat->Draw("E1 same");
   aNoMatch->Draw("E1 same");
   aNoChi2->Draw("E1 same");
   aAny->Draw("E1 same");
   auto *lg2 = new TLegend(0.55, 0.13, 0.89, 0.42);
   lg2->AddEntry(aGood, "nominal (all cuts)", "lp");
   lg2->AddEntry(aNoDth, "no #Delta#theta cut", "lp");
   lg2->AddEntry(aNoRat, "no KE-ratio cut", "lp");
   lg2->AddEntry(aNoMatch, "no truth match", "lp");
   lg2->AddEntry(aNoChi2, "no #chi^{2} cut", "lp");
   lg2->AddEntry(aAny, "any converged fit", "lp");
   lg2->Draw();

   c1->cd(3);
   gPad->SetLogy();
   hDth25->SetLineColor(kRed + 1);
   hDth25->SetLineWidth(2);
   hDth45->SetLineColor(kAzure + 2);
   hDth45->SetLineWidth(2);
   hDth25->DrawNormalized("hist");
   hDth45->DrawNormalized("hist same");
   auto *cutl = new TLine(dThetaMax, 0, dThetaMax, 1);
   cutl->SetLineStyle(2);
   cutl->SetLineColor(kGray + 2);
   auto *lg3 = new TLegend(0.55, 0.7, 0.89, 0.87);
   lg3->AddEntry(hDth25, "#theta_{cm} 25-35 (dip)", "l");
   lg3->AddEntry(hDth45, "#theta_{cm} 45-55 (plateau)", "l");
   lg3->Draw();

   c1->cd(4);
   gPad->SetLogy();
   hRat25->SetLineColor(kRed + 1);
   hRat25->SetLineWidth(2);
   hRat45->SetLineColor(kAzure + 2);
   hRat45->SetLineWidth(2);
   hRat25->DrawNormalized("hist");
   hRat45->DrawNormalized("hist same");
   lg3->Draw();

   TString dir = outDir.Length() ? outDir : TString(gSystem->DirName(gInterpreter->GetCurrentMacroName()));
   TString base = dir + "/acc_split_" + tag + (outTag.Length() ? "_" + outTag : "");
   c1->SaveAs(base + ".png");

   TFile fo(base + ".root", "RECREATE");
   hGen->Write();
   for (int c = 0; c < NCAT; ++c)
      hCat[c]->Write();
   // raw counterfactual numerators too, so seeds can be summed before dividing
   hNoDth->Write();
   hNoRat->Write();
   hNoMatch->Write();
   hNoChi2->Write();
   hAny->Write();
   aGood->Write();
   aNoDth->Write();
   aNoRat->Write();
   aNoMatch->Write();
   aNoChi2->Write();
   aAny->Write();
   hDth25->Write();
   hDth45->Write();
   hRat25->Write();
   hRat45->Write();
   hKE->Write();
   hKEtrue->Write();
   fo.Close();
   printf("\nwrote %s.png and .root\n\n", base.Data());
}
