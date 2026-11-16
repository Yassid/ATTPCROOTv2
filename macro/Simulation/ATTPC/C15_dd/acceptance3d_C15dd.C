/// @file acceptance3d_C15dd.C
/// @brief Efficiency of the C15d chain, differential in true theta_CM AND true vz, per Ex.
///
/// WHY: the analysis' accBin() has NO vz dependence, while the real acceptance plausibly does --
/// how much of a track is seen depends on where in the chamber it started. Every spectral shape
/// extracted from the data moves with the vz window (the continuum exponent swings over ~1.2
/// units against a +-0.1 statistical error), and a vz-dependent acceptance is the obvious
/// suspect. This macro MEASURES that dependence instead of assuming it. If the efficiency comes
/// out flat in vz, the diagnosis is wrong and the cause is elsewhere.
///
/// Needs no new simulation: the existing ex*_sim.root carry the truth vertex already.
///
/// DENOMINATOR: every generated event whose MC deuteron falls in the theta_CM bin, with the
/// vertex inside the analysis' own vz window. theta_CM is computed from the MC deuteron with the
/// SAME expression the data analysis uses (excm.h::thcmProj at Ebeam 191.25), so the simulated
/// and measured angular variables are the same quantity, biases included.
///
/// NUMERATOR: the event has a reconstructed track (AtPatternEvent) with at least minHits hits.
/// This is the PRA-level efficiency; it does NOT yet include the genfit stage, so it is an UPPER
/// bound on the efficiency the yields actually see.
///
/// ⚠ Requires the reconstruction to have been run at the VALIDATED gain
/// (ATTPC.C15d_sim_g150000.par). At the data par's Gain = 10000 the efficiency comes out ~2x too
/// small and nothing in the output looks wrong -- see README.md.
#include "/home/yassid/Data/C15d_statefit/excm.h"

void acceptance3d_C15dd(TString simDir = "/home/yassid/Data/C15dd_sim",
                      TString outFile = "/home/yassid/Data/C15dd_sim/acceptance3d_C15dd.root",
                      double Ebeam = 191.25, int minHits = 20,
                      double vzLo = 0, double vzHi = 1100,
                      double cmLo = 10, double cmStep = 5, int nCM = 10,
                        double vzBinLo = 0, double vzBinHi = 1000, int nVZ = 10)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   const double md = 2.0141017781 * 931.49401;
   Pars P; P.ebeam = Ebeam; P.mBeam = 15.0105993; P.mTarg = 2.014102;
   P.mEject = 2.014102;     P.mResid = 15.0105993;
   P.kcOn = false;          P.ebzOn = true; P.ebzSlope = 11.0; P.ebzPivot = 500;

   // collect the Ex points present on disk
   std::vector<TString> tags;
   void *d = gSystem->OpenDirectory(simDir); const char *e;
   while ((e = gSystem->GetDirEntry(d))) { TString s(e);
      if (s.BeginsWith("ex") && s.EndsWith("_sim.root")) tags.push_back(s(0, s.Length() - 9)); }
   gSystem->FreeDirectory(d); std::sort(tags.begin(), tags.end());

   TFile fo(outFile, "RECREATE");
   printf("\n  Ex      N gen   N reco   eff(all)  |  efficiency per vz slice [mm]\n");
   printf("                                      |");
   for (int v = 0; v < nVZ; v++) printf("  %3.0f-%3.0f", vzBinLo + v * (vzBinHi - vzBinLo) / nVZ,
                                        vzBinLo + (v + 1) * (vzBinHi - vzBinLo) / nVZ);
   printf("\n");
   for (auto &tag : tags) {
      TFile fm(simDir + "/" + tag + "_sim.root"), fr(simDir + "/" + tag + "_reco.root");
      if (fm.IsZombie() || fr.IsZombie()) { printf("  %-8s (missing reco)\n", tag.Data()); continue; }
      auto *tm = (TTree *)fm.Get("cbmsim"); auto *tr = (TTree *)fr.Get("cbmsim");
      if (!tm || !tr) { printf("  %-8s (no tree)\n", tag.Data()); continue; }
      TClonesArray *mc = nullptr, *pa = nullptr;
      // The reco files carry AtRawEvent and AtEventH too and are ~7 GB each; reading every
      // branch makes this macro I/O bound for no reason.  Enable only what is used.
      tm->SetBranchStatus("*", 0); tm->SetBranchStatus("MCTrack*", 1);
      tr->SetBranchStatus("*", 0); tr->SetBranchStatus("AtPatternEvent*", 1);
      tm->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);
      double exTrue = TString(tag(2, tag.Length() - 2)).Atof();

      TH2D *hG = new TH2D("g_" + tag, ";#theta_{CM};v_{z} [mm]", nCM, cmLo, cmLo + nCM * cmStep,
                          nVZ, vzBinLo, vzBinHi);
      TH2D *hN = new TH2D("n_" + tag, ";#theta_{CM};v_{z} [mm]", nCM, cmLo, cmLo + nCM * cmStep,
                          nVZ, vzBinLo, vzBinHi);
      hG->Sumw2(); hN->Sumw2();
      Long64_t n = std::min(tm->GetEntries(), tr->GetEntries());
      long ng = 0, nr = 0;
      for (Long64_t i = 0; i < n; i++) {
         tm->GetEntry(i); tr->GetEntry(i);
         // --- truth ---
         double tcm = -1, vz = -1;
         for (int j = 0; j < mc->GetEntriesFast(); j++) {
            auto *m = (FairMCTrack *)mc->At(j);
            if (m->GetPdgCode() != 1000010020 || m->GetMotherId() != -1) continue;
            double px = m->GetPx() * 1000, py = m->GetPy() * 1000, pz = m->GetPz() * 1000;
            double p = std::sqrt(px * px + py * py + pz * pz), ke = std::sqrt(p * p + md * md) - md;
            double th = std::acos(pz / p) * 180.0 / M_PI;
            vz = m->GetStartZ() * 10.0;
            tcm = thcmProj(P, ebAt(P, vz, true), th, ke);
            break;
         }
         if (tcm < 0 || vz < vzLo || vz > vzHi) continue;
         ng++; hG->Fill(tcm, vz);
         // --- reconstructed ---
         bool got = false;
         if (pa->GetEntriesFast()) {
            auto *pe = (AtPatternEvent *)pa->At(0);
            for (auto &t : pe->GetTrackCand())
               if ((int)t.GetHitArray().size() >= minHits) { got = true; break; }
         }
         if (got) { nr++; hN->Fill(tcm, vz); }
      }
      TH2D *hE = (TH2D *)hN->Clone("eff_" + tag);
      hE->Divide(hN, hG, 1, 1, "B");
      printf("  %5.2f  %7ld  %7ld   %6.3f   |", exTrue, ng, nr, ng ? (double)nr / ng : 0);
      // efficiency projected on vz, summed over theta_CM -- the number the diagnosis turns on
      for (int v = 1; v <= nVZ; v++) {
         double g = hG->Integral(1, nCM, v, v), nn = hN->Integral(1, nCM, v, v);
         printf("  %6.3f", g > 0 ? nn / g : 0.0);
      }
      printf("\n");
      fo.cd(); hG->Write(); hN->Write(); hE->Write();
   }
   fo.Close();
   printf("\n  wrote %s\n", outFile.Data());
}
