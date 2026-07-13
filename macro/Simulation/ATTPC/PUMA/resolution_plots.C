/// @file resolution_plots.C
/// @brief Proper resolution DISTRIBUTIONS (not just summary bars) for UKF vs
///        GENFIT on the PUMA test channel: momentum (dp/p), polar angle, and
///        vertex-z residuals, overlaid with Gaussian fits. Same reading logic
///        as compare_ukf_genfit_test8.C.
/// Run: root -b -q 'resolution_plots.C("pi")'   (or "K")
struct Res {
   std::vector<double> dpFrac, theta, dz;
};

static double medianOf(std::vector<double> v)
{
   if (v.empty()) return 0;
   std::sort(v.begin(), v.end());
   size_t n = v.size();
   return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
static double iqrSigma(std::vector<double> v)
{
   if (v.size() < 4) return 0;
   std::sort(v.begin(), v.end());
   return (v[3 * v.size() / 4] - v[v.size() / 4]) / 1.349; // Gaussian-equivalent sigma
}

// Overlay UKF vs GENFIT with robust median/sigma_IQR in the legend and a
// Gaussian drawn over the robust core (fit tails are non-Gaussian, so we quote
// the IQR sigma, matching compare_ukf_genfit_test8.C).
static void fitDraw(std::vector<double> &vU, std::vector<double> &vG, int nb, double lo, double hi, const char *xt,
                    TString out, TString head, const char *unit)
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(0);
   auto *hU = new TH1F(Form("hU_%s", out.Data()), "", nb, lo, hi);
   auto *hG = new TH1F(Form("hG_%s", out.Data()), "", nb, lo, hi);
   for (double v : vU) hU->Fill(v);
   for (double v : vG) hG->Fill(v);
   double mU = medianOf(vU), sU = iqrSigma(vU), mG = medianOf(vG), sG = iqrSigma(vG);

   auto *c = new TCanvas("c", "", 640, 520);
   c->SetLeftMargin(0.14);
   c->SetBottomMargin(0.13);
   hU->SetLineColor(kAzure + 2);
   hU->SetLineWidth(3);
   hG->SetLineColor(kRed + 1);
   hG->SetLineWidth(3);
   hG->SetLineStyle(2);
   hU->SetTitle(Form("%s;%s;tracks", head.Data(), xt));
   hU->SetMaximum(1.3 * std::max(hU->GetMaximum(), hG->GetMaximum()));
   hU->Draw("hist");
   hG->Draw("hist same");
   // Gaussian over the robust core for visual guidance only
   auto *gU = new TF1("gU", "gaus", mU - 1.5 * sU, mU + 1.5 * sU);
   gU->SetParameters(hU->GetMaximum(), mU, sU);
   hU->Fit(gU, "QR0");
   gU->SetLineColor(kAzure + 2); gU->SetLineWidth(2); gU->Draw("same");
   auto *gG = new TF1("gG", "gaus", mG - 1.5 * sG, mG + 1.5 * sG);
   gG->SetParameters(hG->GetMaximum(), mG, sG);
   hG->Fit(gG, "QR0");
   gG->SetLineColor(kRed + 1); gG->SetLineWidth(2); gG->SetLineStyle(2); gG->Draw("same");

   auto *lg = new TLegend(0.15, 0.70, 0.60, 0.89);
   lg->SetTextFont(62);
   lg->SetTextSize(0.038);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);
   lg->AddEntry(hU, Form("UKF: med %+.1f, #sigma %.1f %s", mU, sU, unit), "l");
   lg->AddEntry(hG, Form("GENFIT: med %+.1f, #sigma %.1f %s", mG, sG, unit), "l");
   lg->Draw();
   c->SaveAs(out);
   printf("  %s : UKF med %+.2f sig %.2f | GENFIT med %+.2f sig %.2f\n", xt, mU, sU, mG, sG);
   delete c;
}

void resolution_plots(TString species = "pi", Double_t testEnergy = -1,
                      TString digiFile = "./data/output_digi_both8.root",
                      TString simFile = "./data/attpcsim.root",
                      TString outdir = "/Users/quantumlab/fair_install/puma_slides/figs")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetTextFont(62);
   gStyle->SetLabelFont(62, "xyz");
   gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   const bool isK = (species == "K" || species == "kaon");
   const double m = isK ? 493.677 : 139.57039;
   const double E0 = (testEnergy > 0) ? testEnergy : (isK ? 0.777 : 0.4);
   const double p0 = std::sqrt(E0 * E0 - (m / 1000) * (m / 1000)) * 1000;
   const double vz0 = 75.0;
   const char *sp = isK ? "K^{+}K^{+}" : "#pi^{+}#pi^{-}";

   TFile fD(digiFile);
   auto *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);
   auto *tS = (TTree *)fS.Get("cbmsim");
   if (!tD || !tS) { printf("missing tree\n"); return; }

   auto *ukfArr = new TClonesArray("AtTrackingEvent");
   auto *gfArr = new TClonesArray("AtTrackingEvent");
   tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);

   Res rU, rG;
   auto grab = [&](TClonesArray *arr, Res &r) {
      if (arr->GetEntries() == 0) return;
      auto *te = (AtTrackingEvent *)arr->At(0);
      for (const auto &ft : te->GetFittedTracks()) {
         const auto &kin = ft->GetKinematics(0);
         double KE = kin.kineticEnergy;
         if (!(KE > 0)) continue;
         double p = std::sqrt(KE * KE + 2 * KE * m);
         r.dpFrac.push_back(100.0 * (p - p0) / p0);
         r.theta.push_back(kin.theta * 180.0 / M_PI - 90.0);
         const auto &props = ft->GetTrackPropertiesStruct();
         double vz = props.initialPositionXtr.Z();
         if (std::abs(vz) < 1e-9) vz = ft->GetVertex(0).Z();
         r.dz.push_back(vz - vz0);
      }
   };

   Long64_t nE = tD->GetEntries();
   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e);
      grab(ukfArr, rU);
      grab(gfArr, rG);
   }
   printf("RES %s: UKF n=%zu  GENFIT n=%zu  (truth |p|=%.1f MeV/c)\n", species.Data(), rU.dpFrac.size(),
          rG.dpFrac.size(), p0);

   fitDraw(rU.dpFrac, rG.dpFrac, 40, -120, 120, "(p_{fit}-p_{truth})/p_{truth} [%]",
           outdir + "/res_" + species + "_dp.png", Form("%s momentum resolution", sp), "%");
   fitDraw(rU.theta, rG.theta, 40, -45, 45, "#theta_{fit}-#theta_{truth} [deg]",
           outdir + "/res_" + species + "_theta.png", Form("%s polar-angle residual", sp), "deg");
   fitDraw(rU.dz, rG.dz, 40, -10, 10, "z_{fit}-z_{truth} [mm]", outdir + "/res_" + species + "_vz.png",
           Form("%s vertex-z residual", sp), "mm");

   printf("RES_DONE %s\n", species.Data());
}
