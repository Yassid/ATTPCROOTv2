/// @file ex_res_C14_hf.C
/// @brief Excitation-energy resolution of a 14C(p,p') simulation sample, per theta_lab slice.
///
/// The acceptance macro answers "was the proton reconstructed at all". This answers "and how
/// well", which is the question the field x pad-pitch matrix exists to settle: 2 mm pads and a
/// 4 or 7 T field are only worth anything if they narrow Ex enough to separate the 6.09 / 6.59 /
/// 6.73 / 7.01 MeV multiplet.
///
/// WHAT IS QUOTED, AND WHAT IS NOT. Median and interquartile range, never a walk-outward FWHM.
/// On the 46Ar campaign that estimator returned 0.765 vs 1.800 MeV for two histograms differing
/// by a CONSTANT shift of every entry -- it is noise on distributions with a narrow core and
/// broad tails, which is exactly what these are. IQR/1.349 is quoted alongside as the gaussian-
/// equivalent sigma so the numbers can be read as a resolution, but the IQR is the measurement.
///
/// The residual is Ex(reconstructed) - Ex(generated), so it carries the CENTROID BIAS as well as
/// the width. Both matter: a 7 T configuration that halves the width but moves the centroid by
/// 0.5 MeV has not resolved anything.
///
/// Every accepted event is also written to a flat TTree ("res"), so slices, cuts and any other
/// statistic can be recomputed later without re-running the chain.
///
///   root -b -q 'ex_res_C14_hf.C("sims_b700/gs_s7021_sim.root","b700_2mm/gs_s7021_b700_2mm_genfit.root","gs_s7021_b700_2mm",0.0)'

#include <algorithm>
#include <tuple>
#include <vector>

static double exr_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics, byte for byte the one in acceptance_C14.C and pp/ex_C14.C:
/// returns {Ex, theta_cm[deg]} from the EJECTILE (KE, theta_lab)
static std::tuple<double, double> exr_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                           double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * exr_omega2(s, m1 * m1, m2 * m2) * exr_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (exr_omega2(s, m1 * m1, m2 * m2) * exr_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

static double exr_quantile(std::vector<double> &v, double q)
{
   if (v.empty())
      return std::numeric_limits<double>::quiet_NaN();
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, q * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

/// @param useXtr  kTRUE for genfit (GetKinematicsXtr is the vertex-corrected energy; plain
///                GetKinematics is the RAW first-cluster value and is the wrong quantity).
/// @param mTargetAmu / mEjectAmu / mResidAmu  the reaction, in amu. The defaults are 14C(p,p'), so
///        every existing caller is unchanged; 14C(d,p)15C is (2.0141018, 1.007825, 15.0105993).
///        The beam is always 14C here. resEx is added to mResidAmu to make the generated residual,
///        and the RECONSTRUCTED excitation is always referred to the ground-state residual -- that
///        is what an experiment does, so the residual below is the error.
void ex_res_C14_hf(TString simFile, TString fitFile, TString tag, Double_t resEx = 0.0, Double_t Ebeam = 159.75,
                   Double_t chi2Cut = 5.0, Bool_t useXtr = kTRUE, TString outDir = "./", Double_t dThetaMax = 10.0,
                   Double_t keRatioMin = 0.5, Double_t keRatioMax = 2.0, Double_t mTargetAmu = 1.007825,
                   Double_t mEjectAmu = 1.007825, Double_t mResidAmu = 14.003242,
                   // BEAM mass in amu; trailing with the old hard-coded 14C value as default, so
                   // existing callers are unchanged. 10Be(t,p)12Be passes 10.0135341.
                   Double_t mBeamAmu = 14.003242)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   const double m_C14 = mBeamAmu * u;             // beam (name kept; it is whatever mBeamAmu says)
   const double m_tgt = mTargetAmu * u;           // target
   const double m_ej = mEjectAmu * u;             // ejectile (the track that is fitted)
   const double m_r0 = mResidAmu * u;             // residual, ground state
   const double m_p = m_ej;                       // the truth-track mass, for the KE from momentum
   const double m_resid = m_r0 + resEx;
   printf("channel: 14C + %.4f u -> %.4f u + %.4f u (residual Ex = %.3f MeV), Ebeam = %.2f MeV\n", mTargetAmu,
          mEjectAmu, mResidAmu, resEx, Ebeam);

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
      printf("\033[1;31mENTRY MISMATCH: sim %lld vs fit %lld -- the fit must be of the UNGATED "
             "reco, otherwise entry i of one file is not entry i of the other.\033[0m\n",
             ts->GetEntries(), tf->GetEntries());
      return;
   }

   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   TFile fo(outDir + "exres_" + tag + ".root", "RECREATE");
   double bExR, bExT, bThT, bThR, bKeT, bKeR, bCmT, bC2n, bZT, bZR;
   TTree out("res", "Ex residuals");
   out.Branch("exReco", &bExR);
   out.Branch("exTrue", &bExT);
   out.Branch("thTrue", &bThT);   // deg, lab
   out.Branch("thReco", &bThR);   // deg, lab
   out.Branch("keTrue", &bKeT);   // MeV
   out.Branch("keReco", &bKeR);   // MeV
   out.Branch("cmTrue", &bCmT);   // deg
   out.Branch("chi2ndf", &bC2n);
   out.Branch("zTrue", &bZT);     // mm
   out.Branch("zReco", &bZR);     // mm

   TH1D *hEx = new TH1D("hEx_" + tag, tag + ";E_{x} reconstructed [MeV];counts", 400, -5, 15);
   TH2D *hExTh = new TH2D("hExTh_" + tag, tag + ";#theta_{lab} true [deg];E_{x} reconstructed [MeV]", 45, 0, 90, 200,
                          -5, 15);

   // theta_lab slices. The recoil proton of 14C(p,p') comes out at theta_lab ~ (180-theta_cm)/2,
   // so the whole physical range is 0-90 deg and the forward end is the low-energy one.
   const int nSl = 7;
   const double slLo[nSl] = {20, 30, 40, 50, 60, 70, 80};
   const double slHi[nSl] = {30, 40, 50, 60, 70, 80, 90};
   std::vector<double> res[nSl], all;

   Long64_t N = ts->GetEntries();
   long nGen = 0, nGood = 0;
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc)
         continue;
      double keT = -1, thT = -1, zT = -1e9;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *t = (AtMCTrack *)mc->At(k);
         if (!t || t->GetPdgCode() != 2212 || t->GetMotherId() != -1)
            continue;
         zT = t->GetStartZ() * 10.0;
         double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000;
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0)
            continue;
         keT = std::sqrt(p * p + m_p * m_p) - m_p;
         thT = std::acos(pz / p);
         break;
      }
      if (keT <= 0)
         continue; // beam-only event
      ++nGen;
      auto [exT, cmT] = exr_kine(m_C14, m_tgt, m_ej, m_resid, Ebeam, thT, keT);

      tf->GetEntry(i);
      if (!te || te->GetEntriesFast() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev)
         continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         const auto &md = ft->GetTrackMetadata();
         double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
         double c2n = ndf > 0 ? chi2 / ndf : 1e9;
         const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
         double ke = kin.kineticEnergy, th = kin.theta;
         if (!(ke > 0 && ke < 1000 && c2n < chi2Cut))
            continue;
         // SAME truth match as the acceptance, and for the same reason: without it the scattered
         // 14C supplies fits that are not this proton, and a resolution measured on those is a
         // measurement of nothing.
         double dth = std::fabs(th * TMath::RadToDeg() - thT * TMath::RadToDeg());
         double r = ke / keT;
         if (dThetaMax > 0 && (dth > dThetaMax || r < keRatioMin || r > keRatioMax))
            continue;

         // Ex is reconstructed with the RESIDUAL MASS OF THE GROUND STATE (m_C14), not m_resid:
         // that is what an experiment does -- it does not know the level in advance -- so the
         // reconstructed value should come out at resEx, and the residual below is the error.
         auto [exR, cmR] = exr_kine(m_C14, m_tgt, m_ej, m_r0, Ebeam, th, ke);
         if (std::isnan(exR))
            continue;

         bExR = exR;
         bExT = exT;
         bThT = thT * TMath::RadToDeg();
         bThR = th * TMath::RadToDeg();
         bKeT = keT;
         bKeR = ke;
         bCmT = cmT;
         bC2n = c2n;
         bZT = zT;
         bZR = ft->GetVertex(0).Z();
         out.Fill();
         hEx->Fill(exR);
         hExTh->Fill(bThT, exR);
         double d = exR - resEx;
         all.push_back(d);
         for (int s = 0; s < nSl; ++s)
            if (bThT >= slLo[s] && bThT < slHi[s])
               res[s].push_back(d);
         ++nGood;
         break; // one fitted proton per event
      }
   }

   printf("\n===== Ex resolution %s  (generated Ex = %.3f MeV, chi2/ndf < %.1f) =====\n", tag.Data(), resEx, chi2Cut);
   printf("truth reactions %ld   truth-matched good fits %ld   (%.1f %%)\n", nGen, nGood,
          nGen ? 100.0 * nGood / nGen : 0.0);
   printf("\n  theta_lab      n   median [MeV]    IQR [MeV]   IQR/1.349   |resid|<0.5\n");
   auto report = [&](const char *label, std::vector<double> &v) {
      if (v.size() < 20) {
         printf("  %-10s %6zu        --            --           --          --\n", label, v.size());
         return;
      }
      double q25 = exr_quantile(v, 0.25), q50 = exr_quantile(v, 0.50), q75 = exr_quantile(v, 0.75);
      long nIn = std::count_if(v.begin(), v.end(), [](double d) { return std::fabs(d) < 0.5; });
      printf("  %-10s %6zu   %+8.3f      %8.3f     %8.3f     %6.1f %%\n", label, v.size(), q50, q75 - q25,
             (q75 - q25) / 1.349, 100.0 * nIn / v.size());
   };
   for (int s = 0; s < nSl; ++s) {
      TString lab = TString::Format("%.0f-%.0f", slLo[s], slHi[s]);
      report(lab.Data(), res[s]);
   }
   report("ALL", all);

   fo.cd();
   out.Write();
   hEx->Write();
   hExTh->Write();
   // the slice summary as a histogram too, so compare_hf_C14.C does not have to re-read the tree
   TH1D *hMed = new TH1D("hMedian_" + tag, tag + ";#theta_{lab} [deg];median E_{x} residual [MeV]", nSl, 20, 90);
   TH1D *hIqr = new TH1D("hIQR_" + tag, tag + ";#theta_{lab} [deg];E_{x} IQR [MeV]", nSl, 20, 90);
   for (int s = 0; s < nSl; ++s) {
      if (res[s].size() < 20)
         continue;
      double q25 = exr_quantile(res[s], 0.25), q50 = exr_quantile(res[s], 0.50), q75 = exr_quantile(res[s], 0.75);
      hMed->SetBinContent(s + 1, q50);
      hMed->SetBinError(s + 1, 1.253 * (q75 - q25) / 1.349 / std::sqrt((double)res[s].size()));
      hIqr->SetBinContent(s + 1, q75 - q25);
      hIqr->SetBinError(s + 1, (q75 - q25) / std::sqrt(2.0 * res[s].size()));
   }
   hMed->Write();
   hIqr->Write();
   fo.Close();

   printf("\nwrote %sexres_%s.root\nex res done\n\n", outDir.Data(), tag.Data());
}
