/// @file ex_core_3Hed.h
/// @brief The one 46Ar(3He,d)47K excitation-energy inversion, shared by the table macro
///        (ex_genfit_3Hed.C) and the plot macro (plots_matrix_3Hed.C).
///
/// WHY A SHARED HEADER RATHER THAN A COPY. The tables and the figures have to be the same
/// measurement. Two implementations of this loop would agree on the day they were written and
/// then drift -- and a figure that disagrees with its own table is worse than no figure, because
/// nothing announces it. Everything physical lives here; the callers only bin and draw.
///
/// Extracted verbatim from ex_genfit_3Hed.C on 2026-09-05 and checked by re-running the full
/// 3-state x 4-configuration matrix and diffing against the pre-refactor output: identical.
///
/// The three things in here that are easy to get wrong, all of them silent:
///
///   * READS GetKinematicsXtr, NOT GetKinematics. With SetBackExtrapToAxis on (the default in
///     fitGenfitter_Ar46.C) the vertex-gap-corrected momentum lands in the Xtr slot and the raw
///     genfit momentum stays in the other. Reading the wrong slot loses the correction quietly.
///   * DRIFT-Z IS MIRRORED in this simulation, so the beam energy at the vertex needs
///     zUse = driftLength - z_reco. The sign is MEASURED per configuration from the correlation
///     against truth, never assumed -- and it must be measured on the vertex, not inferred from
///     the Ex spread, which with three states mixed is dominated by the level separation and is
///     nearly blind to a 100 cm flip.
///   * The track is matched to truth by ANGLE (dThetaMax), so a configuration that fits fewer
///     tracks reports a smaller n rather than a worse resolution. Always read nFit/nTruth beside
///     any width, or a survivorship-biased cell reads as a good one.
#ifndef EX_CORE_3HED_H
#define EX_CORE_3HED_H

#include "ar46_masses.h"   // kMb/kMt/kMR/kMe/kTb0/kdEdz -- the ONE definition, never retyped

#include <algorithm>
#include <cmath>
#include <vector>

namespace Ar46 {

inline double MedOf(std::vector<double> v)
{
   if (v.empty()) return -999;
   std::sort(v.begin(), v.end());
   size_t n = v.size();
   return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
inline double IqrOf(std::vector<double> v)
{
   if (v.size() < 8) return -999;
   std::sort(v.begin(), v.end());
   return v[(size_t)(0.75 * v.size())] - v[(size_t)(0.25 * v.size())];
}
/// IQR -> the FWHM an equally wide Gaussian would have. IQR = 1.349 sigma, FWHM = 2.3548 sigma.
/// Quoted rather than a measured FWHM on purpose: a half-maximum walk out from the peak returned
/// 0.765 against 1.800 MeV for two histograms differing only by a CONSTANT shift of every entry.
inline double IqrToFwhm(double iqr) { return 1.7456 * iqr; }

struct Sample {
   std::vector<double> ex;        ///< excitation energy per accepted track [MeV]
   std::vector<double> thetaTrue; ///< the MC-truth deuteron lab angle of the same track [deg]
   // The RECONSTRUCTED quantities behind each ex, kept so a viewer can recompute Ex live under a
   // different beam energy or vertex treatment instead of being stuck with this file's choice.
   // Same order and same length as ex/thetaTrue -- one entry per accepted track.
   std::vector<double> ke;        ///< fitted deuteron KE [MeV] (the Xtr slot when useXtr)
   std::vector<double> theta;     ///< fitted deuteron lab angle [deg]
   std::vector<double> vz;        ///< fitted vertex z [cm], ALREADY un-mirrored (see mirror)
   std::vector<double> chi2;      ///< chi2/ndf, or -1 when ndf <= 0
   long nTruth = 0, nFit = 0, nCut = 0;
   double rV = 0;        ///< vertex-z correlation against truth
   bool mirror = false;  ///< whether the drift-z mirror was applied (rV < 0)
   bool ok = false;
};

/// One configuration = one fit directory, over every tag given.
/// @param dir   directory holding <tag>_genfitter_d.root
/// @param simDir directory holding <tag>_sim.root (they differ when generation is shared)
inline Sample Collect(TString dir, TString simDir, TObjArray *tags, bool useXtr = true, double dThetaMax = 10.0,
                      double chi2Max = -1.0, double driftLength = 100.0, TString label = "")
{
   Sample S;
   std::vector<double> kz, kth, kT, ktrue, kc;
   double cx = 0, cy = 0, cxy = 0, cxx = 0, cyy = 0;
   long nv = 0;

   for (int it = 0; it < tags->GetEntries(); ++it) {
      TString tg = ((TObjString *)tags->At(it))->GetString();
      TString fs = simDir + "/" + tg + "_sim.root", ff = dir + "/" + tg + "_genfitter_d.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(ff)) {
         printf("  [%s] MISSING %s or its fit -- skipped\n", label.Data(), tg.Data());
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
         ++S.nTruth;
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
            auto &md = ft->GetTrackMetadata();
            bChi = (md && md->GetNdf() > 0) ? md->GetChi2() / md->GetNdf() : -1;
            got = true;
         }
         if (!got) continue;
         cx += zTrue; cy += bZ; cxy += zTrue * bZ; cxx += zTrue * zTrue; cyy += bZ * bZ; ++nv;
         if (bd > dThetaMax) continue;
         if (chi2Max > 0 && bChi > 0 && bChi > chi2Max) { ++S.nCut; continue; }
         ++S.nFit;
         kz.push_back(bZ);
         kth.push_back(bTh);
         kT.push_back(bT);
         ktrue.push_back(thTrue);
         kc.push_back(bChi);
      }
      Fs->Close();
      Ff->Close();
   }

   S.rV = (nv > 2) ? (nv * cxy - cx * cy) / std::sqrt((nv * cxx - cx * cx) * (nv * cyy - cy * cy)) : 0;
   S.mirror = (S.rV < 0);

   for (size_t j = 0; j < kz.size(); ++j) {
      double zUse = S.mirror ? (driftLength - kz[j]) : kz[j];
      double Tb = kTb0 - kdEdz * zUse;
      if (Tb < 50 || Tb > kTb0 + 20) continue;
      double Eb = Tb + kMb, pb = std::sqrt(Tb * (Tb + 2 * kMb));
      double Ed = kT[j] + kMe, pd = std::sqrt(kT[j] * (kT[j] + 2 * kMe));
      double th = kth[j] * TMath::DegToRad();
      double ER = Eb + kMt - Ed;
      double pRz = pb - pd * std::cos(th), pRt = pd * std::sin(th);
      double m2 = ER * ER - pRz * pRz - pRt * pRt;
      if (m2 <= 0) continue;
      S.ex.push_back(std::sqrt(m2) - kMR);
      S.thetaTrue.push_back(ktrue[j]);
      S.ke.push_back(kT[j]);
      S.theta.push_back(kth[j]);
      S.vz.push_back(zUse);   // un-mirrored, so a consumer never has to redo the handedness test
      S.chi2.push_back(kc[j]);
   }
   S.ok = !S.ex.empty();
   return S;
}

} // namespace Ar46
#endif
