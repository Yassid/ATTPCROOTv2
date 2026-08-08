/// @file apply_acceptance_C14.C
/// @brief Correct the a1954 14C(p,p') angular distribution with the simulated acceptance.
///
///   N_corrected(theta_cm) = N_measured(theta_cm) / acceptance(theta_cm)
///
/// The acceptance comes from macro/Simulation/ATTPC/14C_pp (5 seeds, truth-matched numerator,
/// 3 cm beam hole, B = -28.5 kG). It is PER LEVEL: pass the g.s. curve for the elastic peak and
/// the 6.094 MeV curve for the inelastic one -- the two differ by more than a factor 20 at some
/// angles, so using the elastic acceptance for an inelastic yield is simply wrong.
///
/// RULES BAKED IN HERE (see project_a1954_c14_acceptance):
///   * usable range theta_cm 20-150 deg. Below 20 the acceptance is near zero or steeply
///     varying; above 140 the systematic is as large as the correction. Bins outside are
///     dropped, NOT corrected, because dividing by a near-zero acceptance manufactures signal.
///   * 5 deg bins, matching the acceptance binning exactly. The acceptance has a 34 % step at
///     30 deg (g.s.) and doubles at 20 deg (ex1); a coarser bin would need the DATA's within-bin
///     angular distribution to equal the simulation's, which cannot be assumed.
///   * per-region systematic, applied on top of the statistical error:
///       20-30 deg  g.s. +6/-3.5 % one-sided (sim under-produces charge and pads, and both can
///                  only LOSE marginal tracks, so the truth is at or above the simulated value)
///       20-30 deg  ex1  +-4 % symmetric (the two gain settings bracket zero)
///       30-110     +-2 %      110-140  +-7 %
///
/// dsigma/dOmega: a histogram of counts in equal theta bins carries the 2*pi*sin(theta) solid
/// angle factor, so the last panel divides it out. Shape only -- no luminosity or target
/// thickness here, so the normalisation is arbitrary.
///
/// The data cache and the acceptance must come from the SAME fitter in the SAME configuration:
/// UKF data goes with /mnt/f/a1954_C14_acc/, GENFIT data (proton_kin_300gfx.root, refit with
/// matEffects OFF + backExtrap) goes with /mnt/f/a1954_C14_acc_gf/. Mixing them corrects the
/// data with an efficiency it does not have. GENFIT's g.s. sits at +0.115 rather than UKF's
/// -0.059, so its Ex windows must be shifted by +0.174 to select the same physics.
///
/// outTag keeps the two fitters' outputs apart: without it a GENFIT run silently overwrites
/// plots/angdist_gs.{png,root} from the UKF run.
///
///   root -b -q 'apply_acceptance_C14.C("plots/proton_kin_300_ukf.root",-0.6,0.6,"gs","elastic")'
///   root -b -q 'apply_acceptance_C14.C("plots/proton_kin_300_ukf.root",5.5,7.0,"ex1","Ex 5.5-7.0")'
///   root -b -q 'apply_acceptance_C14.C("plots/proton_kin_300gfx.root",-0.485,0.715,"gs",
///                                      "elastic (GENFIT)",5,"/mnt/f/a1954_C14_acc_gf/",20,150,"_gf")'

void apply_acceptance_C14(TString dataCache = "plots/proton_kin_300_ukf.root", Double_t exLo = -0.6,
                          Double_t exHi = 0.6, TString level = "gs", TString label = "elastic",
                          Double_t chi2Cut = 5.0, TString accDir = "/mnt/f/a1954_C14_acc/",
                          Double_t cmMin = 20.0, Double_t cmMax = 150.0, TString outTag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + dataCache);
   if (!fd || fd->IsZombie()) {
      printf("\033[1;31mcannot open data cache %s\033[0m\n", dataCache.Data());
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   TFile *fa = TFile::Open(accDir + "acceptance_merged_" + level + ".root");
   if (!t || !fa || fa->IsZombie()) {
      printf("\033[1;31mmissing ntuple or acceptance file for level %s\033[0m\n", level.Data());
      return;
   }
   auto *acc = (TH1D *)fa->Get("hAcc_" + level + "_sum");
   if (!acc) {
      printf("\033[1;31mno hAcc_%s_sum\033[0m\n", level.Data());
      return;
   }

   // data yield on the acceptance binning -- identical binning is what makes a bin-by-bin
   // division meaningful
   const int nb = acc->GetNbinsX();
   auto *raw = new TH1D("raw", "", nb, acc->GetXaxis()->GetXmin(), acc->GetXaxis()->GetXmax());
   raw->Sumw2();
   TString cut = TString::Format("chi2ndf<%g&&ex>%g&&ex<%g", chi2Cut, exLo, exHi);
   t->Draw("thcm>>raw", cut, "goff");

   auto *cor = (TH1D *)raw->Clone("cor");
   cor->Reset();
   auto *dsdo = (TH1D *)raw->Clone("dsdo");
   dsdo->Reset();

   printf("\n===== %s : %s, Ex in [%.2f, %.2f], level '%s' =====\n", dataCache.Data(), label.Data(), exLo, exHi,
          level.Data());
   printf("  theta_cm     raw    acc            corrected           syst\n");
   double sumRaw = 0, sumCor = 0;
   for (int b = 1; b <= nb; ++b) {
      double lo = raw->GetBinLowEdge(b), wid = raw->GetBinWidth(b), ctr = raw->GetBinCenter(b);
      double N = raw->GetBinContent(b), A = acc->GetBinContent(b), dA = acc->GetBinError(b);
      if (ctr < cmMin || ctr > cmMax || N <= 0)
         continue;
      if (A <= 0.05) { // never divide by a vanishing acceptance
         printf("  %3.0f-%3.0f  %6.0f   acceptance %.3f -- DROPPED\n", lo, lo + wid, N, A);
         continue;
      }
      // per-region systematic on the acceptance
      double sysUp = 0.02, sysDn = 0.02;
      if (ctr < 30) {
         if (level == "gs") { sysUp = 0.060; sysDn = 0.035; }
         else               { sysUp = 0.040; sysDn = 0.040; }
      } else if (ctr > 110) {
         sysUp = sysDn = 0.07;
      }
      double C = N / A;
      double relStat = std::sqrt(1.0 / N + (dA / A) * (dA / A));
      cor->SetBinContent(b, C);
      cor->SetBinError(b, C * relStat);
      // dsigma/dOmega shape: strip the 2*pi*sin(theta) solid-angle factor of an equal-theta bin
      double s = std::sin(ctr * TMath::DegToRad());
      if (s > 1e-3) {
         dsdo->SetBinContent(b, C / s);
         dsdo->SetBinError(b, C * relStat / s);
      }
      sumRaw += N;
      sumCor += C;
      printf("  %3.0f-%3.0f  %6.0f  %.3f+-%.3f  %8.0f +- %-6.0f  +%.0f%%/-%.0f%%\n", lo, lo + wid, N, A, dA, C,
             C * relStat, 100 * sysUp, 100 * sysDn);
   }
   printf("  ---- totals in %.0f-%.0f deg: raw %.0f -> corrected %.0f (factor %.3f)\n", cmMin, cmMax, sumRaw, sumCor,
          sumRaw > 0 ? sumCor / sumRaw : 0);

   TCanvas *c = new TCanvas("c", "acc-corrected", 1400, 470);
   c->Divide(3, 1);
   c->cd(1);
   raw->SetTitle(TString::Format("%s: measured yield;#theta_{cm} [deg];counts", label.Data()));
   raw->GetXaxis()->SetRangeUser(cmMin - 10, cmMax + 10);
   raw->SetLineColor(kGray + 2);
   raw->Draw("hist");
   c->cd(2);
   acc->SetTitle(TString::Format("acceptance (%s);#theta_{cm} [deg];acceptance", level.Data()));
   acc->GetXaxis()->SetRangeUser(cmMin - 10, cmMax + 10);
   acc->SetMarkerStyle(20);
   acc->SetMarkerColor(kAzure + 2);
   acc->SetLineColor(kAzure + 2);
   acc->Draw("E1");
   c->cd(3);
   gPad->SetLogy();
   dsdo->SetTitle(TString::Format("%s: d#sigma/d#Omega (shape);#theta_{cm} [deg];corrected / sin#theta", label.Data()));
   dsdo->GetXaxis()->SetRangeUser(cmMin - 10, cmMax + 10);
   dsdo->SetMarkerStyle(20);
   dsdo->SetMarkerColor(kOrange + 7);
   dsdo->SetLineColor(kOrange + 7);
   dsdo->Draw("E1");
   TString png = here + "/plots/angdist_" + level + outTag + ".png";
   c->SaveAs(png);

   TFile fo(here + "/plots/angdist_" + level + outTag + ".root", "RECREATE");
   raw->Write("raw");
   cor->Write("corrected");
   dsdo->Write("dsigma_dOmega");
   acc->Write("acceptance");
   fo.Close();
   printf("wrote %s and .root\n\n", png.Data());
}
