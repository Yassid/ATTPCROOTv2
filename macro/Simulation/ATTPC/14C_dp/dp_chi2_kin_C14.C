/// @file dp_chi2_kin_C14.C
/// @brief The reconstructed kinematic plane at a series of chi2/ndf cuts.
///
///   root -b -q 'dp_chi2_kin_C14.C()'
///
/// Read straight from the genfit output, so NO truth match and no pre-selection: every fitted
/// track appears, and the only thing changing between panels is the quality cut. The typical
/// chi2/ndf of a good fit in this configuration is 0.02 -- measSigma is 4 mm while the residuals
/// are ~0.6 mm -- so the production cut at 5 sits far above anything meaningful, and what it does
/// or does not remove is worth looking at rather than assuming.

#include <vector>

static const double U = 931.49401;
static const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;
static double ck_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static bool ck_fwd(double m4, double K, double tcm, double &ke, double &thl)
{
   double E1 = K + M1, s = M1 * M1 + M2 * M2 + 2 * M2 * E1, rs = std::sqrt(s);
   if (rs < M3 + m4) return false;
   double pcm = ck_om2(s, M3 * M3, m4 * m4) / (2 * rs), E3cm = std::sqrt(pcm * pcm + M3 * M3);
   double p1 = std::sqrt(E1 * E1 - M1 * M1), beta = p1 / (E1 + M2), g = 1.0 / std::sqrt(1 - beta * beta);
   double th = tcm * TMath::DegToRad();
   ke = g * (E3cm + beta * pcm * std::cos(th)) - M3;
   thl = std::atan2(pcm * std::sin(th), g * (pcm * std::cos(th) + beta * E3cm)) * TMath::RadToDeg();
   return ke > 0;
}

void dp_chi2_kin_C14(TString fitFile = "/mnt/f/a1954_C14dp/b285_attpc/gs_s9001_b285_attpc_genfit.root",
                     Double_t Ebeam = 155.9, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const int NC = 6;
   const double cut[NC] = {1e9, 5.0, 1.0, 0.1, 0.05, 0.02};
   const char *lab[NC] = {"no cut", "#chi^{2}/ndf < 5  (production)", "< 1", "< 0.1", "< 0.05", "< 0.02"};

   TFile *ff = TFile::Open(fitFile);
   TTree *tf = (TTree *)ff->Get("cbmsim");
   TClonesArray *te = nullptr;
   tf->SetBranchAddress("AtTrackingEvent", &te);

   std::vector<TH2D *> h(NC);
   for (int c = 0; c < NC; ++c)
      h[c] = new TH2D(Form("hk%d", c), Form("%s;#theta_{lab} [deg];proton KE [MeV]", lab[c]), 180, 0, 180, 240, 0, 60);
   long n[NC] = {0};

   for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
      tf->GetEntry(i);
      if (!te || te->GetEntriesFast() == 0) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue;
         const auto &md = ft->GetTrackMetadata();
         double ndf = md ? md->GetNdf() : 0, c2 = md ? md->GetChi2() : 0;
         if (ndf <= 0) continue;
         double r = c2 / ndf;
         const auto &kin = ft->GetKinematicsXtr();
         double ke = kin.kineticEnergy, th = kin.theta * TMath::RadToDeg();
         if (ke <= 0 || ke > 1000) continue;
         for (int c = 0; c < NC; ++c)
            if (r < cut[c]) { h[c]->Fill(th, ke); ++n[c]; }
         break;
      }
   }

   auto *loc = new TGraph();
   for (double tcm = 1; tcm <= 179; tcm += 0.4) {
      double ke, thl;
      if (ck_fwd(M4, Ebeam, tcm, ke, thl) && ke < 60) loc->SetPoint(loc->GetN(), thl, ke);
   }
   loc->SetLineColor(kRed + 1);
   loc->SetLineWidth(2);

   auto *cv = new TCanvas("ck", "ck", 1700, 1000);
   cv->Divide(3, 2);
   printf("\n  fitted tracks surviving each cut (no truth match)\n");
   for (int c = 0; c < NC; ++c) {
      cv->cd(c + 1);
      gPad->SetLogz();
      gPad->SetLeftMargin(0.13);
      h[c]->Draw("colz");
      loc->Draw("l same");
      auto *l90 = new TLine(90, 0, 90, 60);
      l90->SetLineStyle(2); l90->SetLineColor(kGray + 2); l90->Draw();
      auto *tx = new TLatex(0.16, 0.84, TString::Format("%ld tracks  (%.0f %% of all)", n[c], 100.0 * n[c] / n[0]));
      tx->SetNDC(); tx->SetTextSize(0.042); tx->Draw();
      printf("  %-32s %6ld  %5.1f %%\n", lab[c], n[c], 100.0 * n[c] / n[0]);
   }
   cv->SaveAs(outDir + "dp_chi2_kin.png");
   printf("\n  wrote %sdp_chi2_kin.png\n\n", outDir.Data());
   ff->Close();
}
