/// @file acceptance_C14.C
/// @brief Detector acceptance vs theta_cm for 14C(p,p'), from MC truth.
///
/// acceptance(theta_cm) = (reactions that end up as a good fitted proton) / (reactions generated)
///
/// The denominator is MC TRUTH, not anything reconstructed: for every generated reaction the
/// true recoil proton is taken from the MCTrack list and its (KE, theta_lab) pushed through the
/// same two-body expressions the analysis uses, giving the true theta_cm. So the denominator is
/// unaffected by any reconstruction inefficiency -- which is the whole point.
///
/// The numerator walks the FIT output entry by entry. That requires the fit file to have one
/// entry per generated event, i.e. it must be the fit of the UNGATED reco. With the 3 cm beam
/// hole in place the beam produces no hits at all, so the ungated reco is already essentially
/// pure recoil protons and gate_truth_C14.C is no longer needed to keep it clean. The macro
/// CHECKS the entry counts match and refuses to guess if they do not.
///
/// Acceptance depends on the level: an inelastic channel puts the same theta_cm at a different
/// theta_lab and a different recoil energy, so a 6.09 MeV run does not have the elastic run's
/// acceptance. Run it per level.
///
///   root -b -q 'acceptance_C14.C("./diagnostics/acc/gs_sim.root","./diagnostics/acc/gs_ukf.root","gs",0.0)'
///   root -b -q 'acceptance_C14.C("./diagnostics/acc/ex1_sim.root","./diagnostics/acc/ex1_ukf.root","ex1",6.094)'

#include <tuple>

static double acc_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics, identical to pp/ex_C14.C: returns {Ex, theta_cm[deg]}
static std::tuple<double, double> acc_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                           double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * acc_omega2(s, m1 * m1, m2 * m2) * acc_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (acc_omega2(s, m1 * m1, m2 * m2) * acc_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

void acceptance_C14(TString simFile, TString fitFile, TString tag = "gs", Double_t resEx = 0.0, Double_t Ebeam = 161.0,
                    Double_t chi2Cut = 5.0, Int_t nBins = 36, Double_t cmMax = 180.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   const double m_C14 = 14.003242 * u, m_p = 1.007825 * u;
   const double m_resid = m_C14 + resEx; // the residual is left excited

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
      printf("\033[1;31mENTRY MISMATCH: sim %lld vs fit %lld -- the fit must be of the UNGATED reco,\n"
             "otherwise entries do not correspond and the acceptance would be nonsense.\033[0m\n",
             ts->GetEntries(), tf->GetEntries());
      return;
   }

   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   TH1D *hGen = new TH1D("hGen_" + tag, "generated;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
   TH1D *hRec = new TH1D("hRec_" + tag, "reconstructed", nBins, 0, cmMax);
   TH1D *hLab = new TH1D("hLab_" + tag, ";#theta_{lab} [deg];reactions", 90, 0, 180);
   TH1D *hLabR = new TH1D("hLabR_" + tag, "", 90, 0, 180);
   hGen->Sumw2();
   hRec->Sumw2();

   Long64_t N = ts->GetEntries();
   long nGen = 0, nRec = 0;
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc)
         continue;
      // the true recoil proton: PDG 2212, primary (mother -1)
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *t = (AtMCTrack *)mc->At(k);
         if (!t || t->GetPdgCode() != 2212 || t->GetMotherId() != -1)
            continue;
         double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000; // GeV -> MeV
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0)
            continue;
         keT = std::sqrt(p * p + m_p * m_p) - m_p;
         thT = std::acos(pz / p);
         break;
      }
      if (keT <= 0)
         continue; // beam-only event: no reaction happened here
      auto [exT, cmT] = acc_kine(m_C14, m_p, m_p, m_resid, Ebeam, thT, keT);
      if (std::isnan(cmT))
         continue;
      ++nGen;
      hGen->Fill(cmT);
      hLab->Fill(thT * TMath::RadToDeg());

      // numerator: a converged fitted proton in the SAME entry
      tf->GetEntry(i);
      bool good = false;
      if (te && te->GetEntriesFast() > 0) {
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (ev)
            for (auto &ft : ev->GetFittedTracks()) {
               if (!ft)
                  continue;
               // GetTrackMetadata() returns a unique_ptr reference, not a raw pointer
               const auto &md = ft->GetTrackMetadata();
               double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
               double c2n = ndf > 0 ? chi2 / ndf : 1e9;
               double ke = ft->GetKinematics().kineticEnergy;
               if (ke > 0 && ke < 1000 && c2n < chi2Cut) {
                  good = true;
                  break;
               }
            }
      }
      if (good) {
         ++nRec;
         hRec->Fill(cmT);
         hLabR->Fill(thT * TMath::RadToDeg());
      }
   }

   printf("\n===== acceptance %s (residual Ex = %.3f MeV) =====\n", tag.Data(), resEx);
   printf("generated reactions %ld   reconstructed %ld   overall acceptance %.3f\n", nGen, nRec,
          nGen ? double(nRec) / nGen : 0.0);
   if (!nGen) {
      printf("\033[1;31mno truth protons found -- wrong branch or beam-only file?\033[0m\n");
      return;
   }

   TH1D *hAcc = (TH1D *)hRec->Clone("hAcc_" + tag);
   hAcc->Divide(hRec, hGen, 1, 1, "B"); // binomial errors
   hAcc->SetTitle(TString::Format("14C(p,p') acceptance, Ex = %.2f MeV;#theta_{cm} [deg];acceptance", resEx));
   hAcc->SetMinimum(0);
   hAcc->SetMaximum(1.05);

   printf("\n  theta_cm      gen     reco   acceptance\n");
   for (int b = 1; b <= hGen->GetNbinsX(); ++b) {
      double g = hGen->GetBinContent(b);
      if (g < 1)
         continue;
      printf("  %3.0f-%3.0f  %7.0f  %7.0f   %.3f +- %.3f\n", hGen->GetBinLowEdge(b),
             hGen->GetBinLowEdge(b) + hGen->GetBinWidth(b), g, hRec->GetBinContent(b), hAcc->GetBinContent(b),
             hAcc->GetBinError(b));
   }

   TString outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/diagnostics/";
   TFile fo(outDir + "acceptance_" + tag + ".root", "RECREATE");
   hGen->Write(); hRec->Write(); hAcc->Write(); hLab->Write(); hLabR->Write();
   fo.Close();

   TCanvas *c = new TCanvas("c", "acc", 1200, 480);
   c->Divide(2, 1);
   c->cd(1);
   hAcc->SetMarkerStyle(20);
   hAcc->SetLineColor(kAzure + 2);
   hAcc->SetMarkerColor(kAzure + 2);
   hAcc->Draw("E1");
   c->cd(2);
   hGen->SetLineColor(kGray + 2);
   hGen->Draw("hist");
   hRec->SetLineColor(kOrange + 7);
   hRec->Draw("hist same");
   TLegend *l = new TLegend(0.5, 0.75, 0.88, 0.88);
   l->AddEntry(hGen, "generated (truth)", "l");
   l->AddEntry(hRec, "reconstructed", "l");
   l->Draw();
   c->SaveAs(outDir + "acceptance_" + tag + ".png");
   printf("\nwrote %sacceptance_%s.{root,png}\n\n", outDir.Data(), tag.Data());
}
