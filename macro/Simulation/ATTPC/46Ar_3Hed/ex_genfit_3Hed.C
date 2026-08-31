/// @file ex_genfit_3Hed.C
/// @brief Excitation-energy resolution of 46Ar(3He,d)47K from the GENFIT fit, by theta_lab band.
///
/// WHY THIS EXISTS. `compare_configs.C`, which produced the field x pad matrix quoted through
/// August, reads <tag>_reco.root and runs AtTools::AtSpyralPID -- the PRE-FIT first-arc circle
/// estimate. No genfit fit entered that table at all. Its headline conclusions, "2 mm pads yes,
/// 3.8 T no" and "the 3.8 T backward residual goes bimodal", are therefore statements about a
/// circle estimator, and the standing open item was to test them against a real fit. This macro
/// is that test: same kinematics, same binning, same median/IQR statistics, genfit input.
///
/// READS GetKinematicsXtr, NOT GetKinematics. With SetBackExtrapToAxis on -- the default in
/// fitGenfitter_Ar46.C -- the vertex-gap-corrected momentum lands in the Xtr slot and the raw
/// genfit momentum stays in the other one. Reading the wrong slot silently loses the correction.
/// Pass useXtr = kFALSE to read the raw slot instead; the macro prints which it used.
///
/// MEDIAN AND IQR, NOT MEAN AND FWHM. These distributions are a narrow core on broad tails. A
/// half-maximum walk out from the peak returned 0.765 against 1.800 MeV for two histograms that
/// differed only by a CONSTANT shift of every entry, and the earlier "matEffects halves the width"
/// claim came from exactly that. Do not reintroduce it.
///
/// DRIFT-Z IS MIRRORED in this simulation (r = -1.000 against truth, z_true + z_reco = 101.2 cm),
/// so the beam energy at the vertex needs zUse = driftLength - z_reco. The sign is not assumed:
/// it is measured per configuration from the correlation and printed.
///
/// The defaults are the CURRENT three-field ground-state comparison (2.0 / 2.85 / 3.8 T on real
/// per-field Magboltz transport). Note simDirs: the 2.85 and 3.8 T arms re-reconstruct from sims
/// that live in the OLD placeholder directories, because generation depends on the B field and
/// never on the .par, so those sims stayed valid when the transport was corrected.
///
///   root -b -q 'ex_genfit_3Hed.C()'                                  // the adopted comparison
#include <algorithm>
#include <vector>

static double medOf(std::vector<double> v)
{
   if (v.empty()) return -999;
   std::sort(v.begin(), v.end());
   size_t n = v.size();
   return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
static double iqrOf(std::vector<double> v)
{
   if (v.size() < 8) return -999;
   std::sort(v.begin(), v.end());
   return v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())];
}

void ex_genfit_3Hed(TString tags = "gs_s3001,gs_s3002",
                    TString dirs = "/mnt/f/ar46_3hed_mb_B20,/mnt/f/ar46_3hed_mb_B285,/mnt/f/ar46_3hed_mb_B38",
                    TString names = "2.0 T,2.85 T,3.8 T",
                    Double_t dThetaMax = 10.0, Double_t driftLength = 100.0, Bool_t useXtr = kTRUE,
                    Double_t chi2Max = -1.0,
                    // The _mb_ arms at 2.85 and 3.8 T re-reconstruct from sims that live in the
                    // OLD directories (generation depends on the field, not the .par, so it is
                    // reused). Their _sim.root is therefore NOT next to their fit. Empty = same
                    // directory as the fit.
                    TString simDirs = "/mnt/f/ar46_3hed_mb_B20,/mnt/f/ar46_3hed_OLD_2.85T_placeholder,"
                                      "/mnt/f/ar46_3hed_OLD_3.8T_placeholder")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");

   // 46Ar + 3He -> d + 47K, masses in MeV
   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957; // beam KE at z=0 and its loss per cm of drift

   TObjArray *aT = tags.Tokenize(","), *aD = dirs.Tokenize(","), *aN = names.Tokenize(",");
   TObjArray *aS = simDirs.Length() ? simDirs.Tokenize(",") : nullptr;
   const int NB = 8;
   double edge[NB + 1] = {50, 60, 70, 80, 90, 100, 110, 125, 145};

   printf("\nExcitation energy from GENFIT, %s slot%s\n", useXtr ? "GetKinematicsXtr (back-extrapolated)" : "GetKinematics (raw)",
          chi2Max > 0 ? Form(", chi2/ndf < %.2f", chi2Max) : ", no chi2 cut");

   for (int ic = 0; ic < aD->GetEntries(); ++ic) {
      TString dir = ((TObjString *)aD->At(ic))->GetString();
      TString nm = (ic < aN->GetEntries()) ? ((TObjString *)aN->At(ic))->GetString() : dir;

      std::vector<double> band[NB], all;
      // vertex handedness, measured not assumed
      double cx = 0, cy = 0, cxy = 0, cxx = 0, cyy = 0;
      long nv = 0, nTruth = 0, nFit = 0, nCut = 0;
      // (z_reco, theta, KE) kept until the handedness is known
      std::vector<double> kz, kth, kT, ktrue;

      for (int it = 0; it < aT->GetEntries(); ++it) {
         TString tg = ((TObjString *)aT->At(it))->GetString();
         TString sdir = (aS && ic < aS->GetEntries()) ? ((TObjString *)aS->At(ic))->GetString() : dir;
         TString fs = sdir + "/" + tg + "_sim.root", ff = dir + "/" + tg + "_genfitter_d.root";
         if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(ff)) {
            printf("  [%s] MISSING %s or its fit -- skipped\n", nm.Data(), tg.Data());
            continue;
         }
         TFile *Fs = TFile::Open(fs), *Ff = TFile::Open(ff);
         TTree *ts = (TTree *)Fs->Get("cbmsim"), *tf = (TTree *)Ff->Get("cbmsim");
         TClonesArray *mc = nullptr, *te = nullptr;
         ts->SetBranchAddress("MCTrack", &mc);
         tf->SetBranchAddress("AtTrackingEvent", &te);

         Long64_t N = std::min(ts->GetEntries(), tf->GetEntries());
         for (Long64_t i = 0; i < N; ++i) {
            ts->GetEntry(i);
            tf->GetEntry(i);
            double thTrue = -1, zTrue = -1;
            for (int k = 0; k < mc->GetEntriesFast(); ++k) {
               auto *p = (AtMCTrack *)mc->At(k);
               if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020) continue;
               double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
               double pp = std::sqrt(px * px + py * py + pz * pz);
               if (pp <= 0) break;
               thTrue = std::acos(pz / pp) * TMath::RadToDeg();
               zTrue = p->GetStartZ();
               break;
            }
            if (thTrue < 0) continue;
            ++nTruth;
            if (!te || !te->GetEntriesFast()) continue;
            auto *ev = (AtTrackingEvent *)te->At(0);
            if (!ev) continue;

            double bd = 1e9, bTh = 0, bT = 0, bZ = 0, bChi = -1;
            bool got = false;
            for (auto &ft : ev->GetFittedTracks()) {
               if (!ft) continue;
               auto &kk = ft->GetKinematics();
               if (kk.kineticEnergy <= 0) continue;
               double th = kk.theta * TMath::RadToDeg();
               double d = std::fabs(th - thTrue);
               if (d >= bd) continue;
               double ke = useXtr ? ft->GetKinematicsXtr().kineticEnergy : kk.kineticEnergy;
               if (!(ke > 0)) continue;
               bd = d; bTh = th; bT = ke;
               bZ = ft->GetVertex().Z() / 10.0; // mm -> cm
               // chi2/ndf lives on the fit metadata, not on TrackProperties
               auto &md = ft->GetTrackMetadata();
               bChi = (md && md->GetNdf() > 0) ? md->GetChi2() / md->GetNdf() : -1;
               got = true;
            }
            if (!got) continue;
            // Handedness of the fitted vertex, MEASURED against truth. This must not be inferred
            // from the Ex spread: with the three states mixed, that spread is ~2 MeV of level
            // separation and is nearly blind to a 100 cm z flip.
            cx += zTrue; cy += bZ; cxy += zTrue * bZ; cxx += zTrue * zTrue; cyy += bZ * bZ; ++nv;
            if (bd > dThetaMax) continue;
            if (chi2Max > 0 && bChi > 0 && bChi > chi2Max) { ++nCut; continue; }
            ++nFit;
            kz.push_back(bZ);
            kth.push_back(bTh);
            kT.push_back(bT);
            ktrue.push_back(thTrue);
         }
         Fs->Close();
         Ff->Close();
      }

      double rV = (nv > 2) ? (nv * cxy - cx * cy) / std::sqrt((nv * cxx - cx * cx) * (nv * cyy - cy * cy)) : 0;
      bool useMirror = (rV < 0);

      for (size_t j = 0; j < kz.size(); ++j) {
         double zUse = useMirror ? (driftLength - kz[j]) : kz[j];
         double Tb = Tb0 - dEdz * zUse;
         if (Tb < 50 || Tb > Tb0 + 20) continue;
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double Ed = kT[j] + M_e, pd = std::sqrt(kT[j] * (kT[j] + 2 * M_e));
         double th = kth[j] * TMath::DegToRad();
         double ER = Eb + M_t - Ed;
         double pRz = pb - pd * std::cos(th), pRt = pd * std::sin(th);
         double m2 = ER * ER - pRz * pRz - pRt * pRt;
         if (m2 <= 0) continue;
         double ex = std::sqrt(m2) - M_R;
         all.push_back(ex);
         for (int b = 0; b < NB; ++b)
            if (ktrue[j] >= edge[b] && ktrue[j] < edge[b + 1]) band[b].push_back(ex);
      }

      printf("\n===== %s  (%s) =====\n", nm.Data(), dir.Data());
      printf("  truth deuterons %ld | matched fits %ld (%.1f%%)%s | drift-z mirror %s\n", nTruth, nFit,
             100.0 * nFit / std::max(1L, nTruth),
             chi2Max > 0 ? Form(" | chi2-cut %ld", nCut) : "", useMirror ? "APPLIED" : "not applied");
      printf("  vertex-z correlation against truth r = %+.3f  (negative => the sim mirrors drift z)\n", rV);
      printf("  %-12s %7s %10s %10s\n", "theta_lab", "n", "med Ex", "IQR");
      for (int b = 0; b < NB; ++b) {
         if (band[b].size() < 8) continue;
         printf("  %4.0f-%-7.0f %7zu %10.3f %10.3f\n", edge[b], edge[b + 1], band[b].size(), medOf(band[b]),
                iqrOf(band[b]));
      }
      printf("  %-12s %7zu %10.3f %10.3f\n", "ALL", all.size(), medOf(all), iqrOf(all));
   }
   printf("\nNOTE: these tags mix the three states (0, 0.36, 2.02 MeV), so med Ex is a mixture\n"
          "centroid and only the IQR is a resolution. Run one state's tags alone for a centroid.\n");
}
