/// @file response_C15dd.C
/// @brief The Ex RESPONSE FUNCTION of the C15d chain, per true theta_CM bin and per Ex.
///
/// The data's kinematics tree takes ke and theta from the genfit-FITTED sample. This macro uses
/// the AtPIDTask output instead (AtSpyralResult: brho, polar, vertex), which is a port of
/// Spyral's own estimator -- so what comes out here is the response of the SPYRAL-LIKE chain, and
/// is the right comparison for the Spyral column of LEVELS_15C.md. The genfit response needs the
/// fit stage and is not covered here.
///
///   Brho -> p = Brho * 299.792458 MeV/c for z = 1, then KE from the deuteron mass.
///   polar -> theta_lab.  vertex.Z() -> vz, which feeds the analysis' own Ebeam(vz).
///
/// Output: Ex_reco - Ex_true per (Ex, theta_CM bin), fitted with the SAME shape the data is
/// fitted with (bifurcated Gaussian + exponential tails, lineshape.C), so the simulated and
/// measured response are described by the same function and can be compared parameter by
/// parameter.
#include "/home/yassid/Data/C15d_statefit/excm.h"
#include "/home/yassid/Data/C15d_statefit/lineshape.C"

void response_C15dd(TString simDir = "/home/yassid/Data/C15dd_sim",
                    TString outFile = "/home/yassid/Data/C15dd_sim/response_C15dd.root",
                    double Ebeam = 191.25, bool applyThCorr = false,
                    double cmLo = 10, double cmStep = 10, int nCM = 6)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   const double md = 2.0141017781 * 931.49401, CKIN = 299.792458;
   Pars P; P.ebeam = Ebeam; P.mBeam = 15.0105993; P.mTarg = 2.014102;
   P.mEject = 2.014102;     P.mResid = 15.0105993;
   P.ebzOn = true; P.ebzSlope = 11.0; P.ebzPivot = 500;
   P.kcOn = applyThCorr; P.kcDenom = 2535.2; P.kcPivot = 14.5;   // the Spyral-column values

   std::vector<TString> tags;
   void *d = gSystem->OpenDirectory(simDir); const char *e;
   while ((e = gSystem->GetDirEntry(d))) { TString s(e);
      if (s.BeginsWith("ex") && s.EndsWith("_reco.root")) tags.push_back(s(0, s.Length() - 10)); }
   gSystem->FreeDirectory(d); std::sort(tags.begin(), tags.end());

   TFile fo(outFile, "RECREATE");
   printf("\n  response of the Spyral-like estimator:  Ex_reco - Ex_true\n");
   printf("  Ex      N      bias     FWHM_L  FWHM_R   alpha_L alpha_R   (tailed fit)\n");
   for (auto &tag : tags) {
      TFile fm(simDir + "/" + tag + "_sim.root"), fr(simDir + "/" + tag + "_reco.root");
      if (fm.IsZombie() || fr.IsZombie()) continue;
      auto *tm = (TTree *)fm.Get("cbmsim"); auto *tr = (TTree *)fr.Get("cbmsim");
      if (!tm || !tr) continue;
      TClonesArray *mc = nullptr, *pi = nullptr;
      if (!tr->GetBranch("AtPIDEvent")) { printf("  %-8s no AtPIDEvent\n", tag.Data()); continue; }
      // ~7 GB per reco file (AtRawEvent and AtEventH are persisted); read only what is used.
      tm->SetBranchStatus("*", 0); tm->SetBranchStatus("MCTrack*", 1);
      tr->SetBranchStatus("*", 0); tr->SetBranchStatus("AtPIDEvent*", 1);
      tm->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPIDEvent", &pi);
      double exTrue = TString(tag(2, tag.Length() - 2)).Atof();

      TH1D *hR = new TH1D("res_" + tag, ";E_{x}^{reco} - E_{x}^{true} (MeV);", 300, -4, 4);
      std::vector<TH1D *> hCM(nCM);
      for (int b = 0; b < nCM; b++)
         hCM[b] = new TH1D(Form("res_%s_cm%d", tag.Data(), b), "", 200, -4, 4);

      Long64_t n = std::min(tm->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < n; i++) {
         tm->GetEntry(i); tr->GetEntry(i);
         // truth theta_CM
         double tcm = -1, vzT = -1;
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *m = (FairMCTrack *)mc->At(j);
            if (m->GetPdgCode() != 1000010020 || m->GetMotherId() != -1) continue;
            double px = m->GetPx() * 1000, py = m->GetPy() * 1000, pz = m->GetPz() * 1000;
            double p = std::sqrt(px * px + py * py + pz * pz), ke = std::sqrt(p * p + md * md) - md;
            vzT = m->GetStartZ() * 10.0;
            tcm = thcmProj(P, ebAt(P, vzT, true), std::acos(pz / p) * 180.0 / M_PI, ke);
            break;
         }
         if (tcm < 0) continue;
         if (!pi->GetEntriesFast()) continue;
         auto *pe = (AtPIDEvent *)pi->At(0);
         for (auto &sr : pe->GetSpyral()) {
            if (!sr.valid) continue;
            double pMeV = sr.brho * CKIN;                       // z = 1
            double ke = std::sqrt(pMeV * pMeV + md * md) - md;
            double th = sr.polar * 180.0 / M_PI; if (th > 90) th = 180 - th;
            double vz = sr.vertex.Z();
            double thc = P.kcOn ? thCorr(P, th, ke) : th;
            double ex = exOf(P, ebAt(P, vz, true), thc, ke);
            if (!std::isfinite(ex)) continue;
            hR->Fill(ex - exTrue);
            int b = (int)((tcm - cmLo) / cmStep);
            if (b >= 0 && b < nCM) hCM[b]->Fill(ex - exTrue);
            break;                                              // one per event
         }
      }
      // fit with the SAME tailed shape used on the data
      double pk = hR->GetBinContent(hR->GetMaximumBin());
      TF1 gt("gt_" + tag, f_GT, -3, 3, 8);
      gt.SetParameters(pk, 0, 0.3, 0.3, 1.5, 1.5, 1, 0);
      gt.SetParLimits(1, -2, 2); gt.SetParLimits(2, 0.02, 2.0); gt.SetParLimits(3, 0.02, 2.0);
      gt.SetParLimits(4, 0.3, 6.0); gt.SetParLimits(5, 0.3, 6.0);
      hR->Fit(&gt, "RQ0NL");
      printf("  %5.2f %7.0f  %+7.3f  %7.3f %7.3f  %7.2f %7.2f\n", exTrue, hR->GetEntries(),
             gt.GetParameter(1), 2.355 * gt.GetParameter(2), 2.355 * gt.GetParameter(3),
             gt.GetParameter(4), gt.GetParameter(5));
      fo.cd(); hR->Write(); for (auto *h : hCM) h->Write();
   }
   fo.Close();
   printf("\n  wrote %s\n", outFile.Data());
   printf("  compare with the DATA response (lineshape.C on dd_ex.root):\n");
   printf("    ATTPCROOT  FWHM_L/R 0.752 / 1.142   alpha_L 1.19  alpha_R 1.50\n");
   printf("    Spyral     FWHM_L/R 0.497 / 0.523   alpha_L 2.07  alpha_R 1.10\n");
}
