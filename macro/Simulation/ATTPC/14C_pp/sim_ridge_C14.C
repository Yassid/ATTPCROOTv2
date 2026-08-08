/// @file sim_ridge_C14.C
/// @brief Control: apply the DATA's ridge estimator to the simulation.
///
/// kine_plane_C14.C measures the elastic ridge in the data as the MODE of the KE distribution in
/// each theta_lab slice, and finds it +1 to +3 MeV above the elastic kinematic line for
/// theta_lab 31-61 deg. fit_truth_bias_C14.C measures a truth-matched mean bias of only +1-2 % on
/// the simulation. Those two numbers are not directly comparable: one is a mode against a line,
/// the other a mean against truth, and with a KE resolution of several MeV a mode can be pulled.
///
/// So this macro runs the DATA estimator on the SIMULATION -- same slices, same mode-finding,
/// same kinematic line at the same nominal beam energy -- turning the comparison into a like-for-
/// like one. If the sim ridge also sits high, the offset is an estimator/resolution effect and
/// says nothing about calibration; if it sits on the line, the data offset is real.
///
///   root -b -q 'sim_ridge_C14.C("/mnt/f/a1954_C14_acc/","ukf",kFALSE)'
///   root -b -q 'sim_ridge_C14.C("/mnt/f/a1954_C14_acc_gf/","genfit",kTRUE)'

#include <tuple>

static double sr_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double sr_ex(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * sr_om2(s, m1 * m1, m2 * m2) * sr_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

static double sr_keFor(double m1, double m2, double Eb, double thl, double exW)
{
   double lo = 0.02, hi = 0.999 * Eb;
   auto g = [&](double k) { return sr_ex(m1, m2, m2, m1, Eb, thl, k) - exW; };
   double flo = g(lo), fhi = g(hi);
   if (!std::isfinite(flo) || !std::isfinite(fhi) || flo * fhi > 0)
      return NAN;
   for (int i = 0; i < 120; ++i) {
      double mid = 0.5 * (lo + hi), fm = g(mid);
      if (!std::isfinite(fm))
         return NAN;
      if (flo * fm <= 0) {
         hi = mid;
      } else {
         lo = mid;
         flo = fm;
      }
   }
   return 0.5 * (lo + hi);
}

void sim_ridge_C14(TString accDir = "/mnt/f/a1954_C14_acc/", TString fitter = "ukf", Bool_t useXtr = kFALSE,
                   Double_t Eb = 161.0, TString tag = "gs", Int_t s0 = 1001, Int_t s1 = 1005)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   // flat (ke, theta) of every fitted track, exactly what the data ntuple holds
   std::vector<float> vKE, vTh, vKEt, vTht;
   for (int s = s0; s <= s1; ++s) {
      TString sf = TString::Format("%s%s_s%d_sim.root", accDir.Data(), tag.Data(), s);
      TString ff = TString::Format("%s%s_s%d_%s.root", accDir.Data(), tag.Data(), s, fitter.Data());
      TFile *fs = TFile::Open(sf), *ff2 = TFile::Open(ff);
      if (!fs || fs->IsZombie() || !ff2 || ff2->IsZombie())
         continue;
      TTree *ts = (TTree *)fs->Get("cbmsim"), *tf = (TTree *)ff2->Get("cbmsim");
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
               thT = std::acos(pz / p) * TMath::RadToDeg();
               break;
            }
         tf->GetEntry(i);
         if (!te || te->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
            double ke = kin.kineticEnergy, th = kin.theta * TMath::RadToDeg();
            if (!(ke > 0.2 && ke < 200 && th > 5 && th < 100))
               continue;
            vKE.push_back(ke);
            vTh.push_back(th);
            vKEt.push_back(keT);
            vTht.push_back(thT);
         }
      }
      fs->Close();
      ff2->Close();
   }
   printf("\n===== %s SIM ridge vs the elastic line (Ebeam = %.1f), %zu fitted tracks =====\n", fitter.Data(), Eb,
          vKE.size());
   printf("  theta_lab | KE_line |  RECO ridge   diff    rel  |  TRUTH ridge  diff   rel   |  N\n");
   for (double th = 22; th < 84; th += 3) {
      double keL = sr_keFor(m1, m2, Eb, (th + 1.5) * TMath::DegToRad(), 0.0);
      if (!std::isfinite(keL))
         continue;
      TH1D hR("hR", "", 160, std::max(0.0, 0.55 * keL), 1.55 * keL);
      TH1D hT("hT", "", 160, std::max(0.0, 0.55 * keL), 1.55 * keL);
      hR.SetDirectory(nullptr);
      hT.SetDirectory(nullptr);
      for (size_t e = 0; e < vKE.size(); ++e) {
         if (vTh[e] >= th && vTh[e] < th + 3)
            hR.Fill(vKE[e]);
         if (vTht[e] >= th && vTht[e] < th + 3)
            hT.Fill(vKEt[e]);
      }
      if (hR.Integral() < 60)
         continue;
      hR.Smooth(1);
      hT.Smooth(1);
      double pR = hR.GetBinCenter(hR.GetMaximumBin());
      double pT = hT.Integral() > 30 ? hT.GetBinCenter(hT.GetMaximumBin()) : NAN;
      printf("  %4.0f-%-4.0f | %7.2f | %10.2f %+7.2f %+6.1f%% | %10.2f %+6.2f %+6.1f%% | %5.0f\n", th, th + 3, keL, pR,
             pR - keL, 100 * (pR - keL) / keL, pT, pT - keL, 100 * (pT - keL) / keL, hR.Integral());
   }
   printf("\n");
}
