/// @file thetacm_res_C14.C
/// @brief theta_cm resolution of the reconstruction, from MC truth.
///
/// The FRESCO/DWBA curve is a zero-resolution calculation with a sharp diffraction minimum. The
/// measured distribution is that curve convolved with the experimental theta_cm resolution, which
/// fills the minimum in and flattens the secondary maximum. Comparing raw data to a raw DWBA
/// therefore overstates the disagreement exactly where the curve varies fastest -- which is the
/// 60-85 deg region where the sideband extraction currently overshoots by 1.4-2.8x.
///
/// theta_cm here is NOT a measured angle: the analysis computes it from (KE, theta_lab) through
/// two-body kinematics, so its resolution inherits the KE resolution as well as the angular one,
/// and is much worse than theta_lab's ~1 deg. This macro measures it properly: reconstruct
/// theta_cm from the fitted (KE, theta) and from the true (KE, theta) of the same proton, and
/// histogram the difference in bins of true theta_cm.
///
///   root -b -q 'thetacm_res_C14.C("/mnt/f/a1954_C14_acc/","ukf",kFALSE)'
///   root -b -q 'thetacm_res_C14.C("/mnt/f/a1954_C14_acc_gf/","genfit",kTRUE)'

#include <tuple>

static double tr_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static std::tuple<double, double> tr_kine(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * tr_om2(s, m1 * m1, m2 * m2) * tr_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   if (a <= 0)
      return {NAN, NAN};
   double m4e = std::sqrt(a);
   double t = m2 * m2 + m4e * m4e - 2 * m2 * Et4;
   double ct = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4e * m4e) +
                (m1 * m1 - m2 * m2) * (m3 * m3 - m4e * m4e)) /
               (tr_om2(s, m1 * m1, m2 * m2) * tr_om2(s, m3 * m3, m4e * m4e));
   ct = std::max(-1.0, std::min(1.0, ct));
   return {m4e - m4, (TMath::Pi() - std::acos(ct)) * TMath::RadToDeg()};
}

void thetacm_res_C14(TString accDir = "/mnt/f/a1954_C14_acc/", TString fitter = "ukf", Bool_t useXtr = kFALSE,
                     Double_t Eb = 161.0, TString tag = "gs", Int_t s0 = 1001, Int_t s1 = 1005)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   const int NB = 28;
   auto *pRes = new TProfile("pRes", "", NB, 10, 150, "s");
   auto *h2 = new TH2D("h2", "", NB, 10, 150, 120, -30, 30);
   std::vector<std::vector<double>> dl(NB);

   long n = 0;
   for (int s = s0; s <= s1; ++s) {
      TFile *fs = TFile::Open(TString::Format("%s%s_s%d_sim.root", accDir.Data(), tag.Data(), s));
      TFile *ff = TFile::Open(TString::Format("%s%s_s%d_%s.root", accDir.Data(), tag.Data(), s, fitter.Data()));
      if (!fs || fs->IsZombie() || !ff || ff->IsZombie())
         continue;
      TTree *ts = (TTree *)fs->Get("cbmsim"), *tf = (TTree *)ff->Get("cbmsim");
      if (!ts || !tf || ts->GetEntries() != tf->GetEntries())
         continue;
      TClonesArray *mc = nullptr, *te = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);
      for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
         ts->GetEntry(i);
         double keT = -1, thT = -1;
         if (mc)
            for (int k = 0; k < mc->GetEntriesFast(); ++k) {
               auto *tr = (AtMCTrack *)mc->At(k);
               if (!tr || tr->GetPdgCode() != 2212 || tr->GetMotherId() != -1)
                  continue;
               double px = tr->GetPx() * 1000, py = tr->GetPy() * 1000, pz = tr->GetPz() * 1000;
               double p = std::sqrt(px * px + py * py + pz * pz);
               if (p <= 0)
                  continue;
               keT = std::sqrt(p * p + m2 * m2) - m2;
               thT = std::acos(pz / p);
               break;
            }
         if (keT <= 0)
            continue;
         auto [exT, cmT] = tr_kine(m1, m2, m2, m1, Eb, thT, keT);
         if (!std::isfinite(cmT))
            continue;
         tf->GetEntry(i);
         if (!te || te->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         double bestD = 1e9, bKE = -1, bTh = -1;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
            double ke = kin.kineticEnergy, th = kin.theta;
            if (!(ke > 0 && ke < 1000))
               continue;
            double d = std::fabs(th - thT);
            if (d < bestD) {
               bestD = d;
               bKE = ke;
               bTh = th;
            }
         }
         if (bKE < 0 || bestD > 10 * TMath::DegToRad())
            continue;
         auto [exR, cmR] = tr_kine(m1, m2, m2, m1, Eb, bTh, bKE);
         if (!std::isfinite(cmR))
            continue;
         pRes->Fill(cmT, cmR - cmT);
         h2->Fill(cmT, cmR - cmT);
         int b = (int)((cmT - 10) / ((150.0 - 10) / NB));
         if (b >= 0 && b < NB)
            dl[b].push_back(cmR - cmT);
         ++n;
      }
      fs->Close();
      ff->Close();
   }

   printf("\n===== %s : theta_cm resolution from truth, %ld protons =====\n", fitter.Data(), n);
   printf("  theta_cm_true |   bias    rms   |  robust sigma (IQR/1.349)   N\n");
   for (int b = 0; b < NB; ++b) {
      if (dl[b].size() < 40)
         continue;
      auto v = dl[b];
      std::sort(v.begin(), v.end());
      double q1 = v[(size_t)(0.25 * v.size())], q3 = v[(size_t)(0.75 * v.size())];
      double med = v[v.size() / 2];
      printf("  %5.0f-%-5.0f   | %+6.2f %7.2f  |         %8.2f          %6zu\n", 10 + b * 5.0, 15 + b * 5.0, med,
             pRes->GetBinError(b + 1), (q3 - q1) / 1.349, v.size());
   }

   TCanvas *c1 = new TCanvas("c1", "thetacm res", 1200, 500);
   c1->Divide(2, 1);
   c1->cd(1);
   gPad->SetLogz();
   h2->SetTitle(TString::Format("%s: #theta_{cm} reconstruction;#theta_{cm}^{true} [deg];"
                                "#theta_{cm}^{reco} - #theta_{cm}^{true} [deg]",
                                fitter.Data()));
   h2->Draw("colz");
   c1->cd(2);
   pRes->SetTitle(TString::Format("%s: bias and spread;#theta_{cm}^{true} [deg];#Delta#theta_{cm} [deg]",
                                  fitter.Data()));
   pRes->SetMarkerStyle(20);
   pRes->SetLineWidth(2);
   pRes->GetYaxis()->SetRangeUser(-20, 20);
   pRes->Draw();
   auto *z = new TLine(10, 0, 150, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString png = here + "/diagnostics/thetacm_res_" + fitter + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/diagnostics/thetacm_res_" + fitter + ".root", "RECREATE");
   pRes->Write("dThetaCm");
   h2->Write("dThetaCm2D");
   fo.Close();
   printf("\nwrote %s\n\n", png.Data());
}
