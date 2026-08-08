/// @file response_matrix_C14.C
/// @brief Build the theta_cm migration matrix from the sim and store it for folding into DWBA.
///
/// theta_cm is not measured, it is computed from (KE, theta_lab), so it inherits the KE resolution.
/// thetacm_res_C14.C showed the CORE is sharp (robust sigma 0.2-1.4 deg) but the RMS is 2-11 deg:
/// heavy tails. That distinction matters enormously here, because the elastic cross section spans
/// a factor ~350 between the forward peak (~500 mb/sr at 20-25 deg) and the diffraction minimum
/// (~1.5 mb/sr at 62 deg). A 0.3 % tail leaking out of the forward peak doubles the yield at the
/// minimum. So comparing an acceptance-corrected measurement to a RAW DWBA curve is not a fair
/// test anywhere the curve varies fast.
///
/// This writes P(theta_reco | theta_true) as a normalised TH2D, taken directly from the truth-
/// matched sim (which was generated flat in CM, so each true-angle row is an unbiased sample of
/// the response). fold_dwba_C14.C then applies it.
///
///   root -b -q 'response_matrix_C14.C("/mnt/f/a1954_C14_acc/","ukf",kFALSE)'
///   root -b -q 'response_matrix_C14.C("/mnt/f/a1954_C14_acc_gf/","genfit",kTRUE)'

#include <tuple>

static double rm_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static std::tuple<double, double> rm_kine(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * rm_om2(s, m1 * m1, m2 * m2) * rm_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   if (a <= 0)
      return {NAN, NAN};
   double m4e = std::sqrt(a);
   double t = m2 * m2 + m4e * m4e - 2 * m2 * Et4;
   double ct = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4e * m4e) +
                (m1 * m1 - m2 * m2) * (m3 * m3 - m4e * m4e)) /
               (rm_om2(s, m1 * m1, m2 * m2) * rm_om2(s, m3 * m3, m4e * m4e));
   ct = std::max(-1.0, std::min(1.0, ct));
   return {m4e - m4, (TMath::Pi() - std::acos(ct)) * TMath::RadToDeg()};
}

void response_matrix_C14(TString accDir = "/mnt/f/a1954_C14_acc/", TString fitter = "ukf", Bool_t useXtr = kFALSE,
                         Double_t Eb = 161.0, TString tag = "gs", Int_t s0 = 1001, Int_t s1 = 1005, Int_t nb = 36)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   auto *R = new TH2D("R", "P(#theta_{reco}|#theta_{true});#theta_{true} [deg];#theta_{reco} [deg]", nb, 0, 180, nb, 0,
                      180);
   auto *gen = new TH1D("genTrue", "", nb, 0, 180);
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
         auto [exT, cmT] = rm_kine(m1, m2, m2, m1, Eb, thT, keT);
         if (!std::isfinite(cmT))
            continue;
         gen->Fill(cmT);
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
         auto [exR, cmR] = rm_kine(m1, m2, m2, m1, Eb, bTh, bKE);
         if (!std::isfinite(cmR))
            continue;
         R->Fill(cmT, cmR);
         ++n;
      }
      fs->Close();
      ff->Close();
   }

   // normalise each TRUE column to unit probability (over the reco axis)
   for (int bx = 1; bx <= nb; ++bx) {
      double s = 0;
      for (int by = 1; by <= nb; ++by)
         s += R->GetBinContent(bx, by);
      if (s <= 0)
         continue;
      for (int by = 1; by <= nb; ++by)
         R->SetBinContent(bx, by, R->GetBinContent(bx, by) / s);
   }

   printf("\n===== %s : theta_cm response, %ld truth-matched protons =====\n", fitter.Data(), n);
   printf("  theta_true |  diagonal   within +-1 bin  |  leakage beyond +-2 bins\n");
   for (int bx = 1; bx <= nb; ++bx) {
      if (gen->GetBinContent(bx) < 40)
         continue;
      double d = R->GetBinContent(bx, bx), near = 0, far = 0;
      for (int by = 1; by <= nb; ++by) {
         double v = R->GetBinContent(bx, by);
         if (std::abs(by - bx) <= 1)
            near += v;
         if (std::abs(by - bx) > 2)
            far += v;
      }
      printf("  %4.0f-%-4.0f  |  %7.3f   %11.3f  |  %12.4f\n", R->GetXaxis()->GetBinLowEdge(bx),
             R->GetXaxis()->GetBinUpEdge(bx), d, near, far);
   }

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TCanvas *c1 = new TCanvas("c1", "response", 700, 620);
   gPad->SetLogz();
   R->Draw("colz");
   c1->SaveAs(here + "/diagnostics/response_" + fitter + ".png");
   TFile fo(here + "/diagnostics/response_" + fitter + ".root", "RECREATE");
   R->Write("response");
   gen->Write("genTrue");
   fo.Close();
   printf("\nwrote diagnostics/response_%s.{png,root}\n\n", fitter.Data());
}
