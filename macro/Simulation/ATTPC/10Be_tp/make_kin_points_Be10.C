/// @file make_kin_points_Be10.C
/// @brief Cache the (theta_lab, KE) points the kinematics GUI needs, so it opens instantly.
///
///   root -b -q 'make_kin_points_Be10.C("/mnt/f/Be10_tp")'
///
/// The GENERATED truth has to come from the sim files, and walking 4 x 16000 events of MCTrack
/// takes long enough that doing it on every GUI redraw would make the thing unusable. It is also
/// field-independent -- the generator does not know about B -- so one pass over the 2.85 T sims
/// serves every configuration.
///
/// The ACCEPTED and RECONSTRUCTED samples are already small (a few thousand rows per level) and
/// live in the campaign's exres trees, so the GUI reads those directly and they are not cached
/// here; that way a rerun of a configuration is picked up without rebuilding anything.
///
/// Writes plots/kin_truth_Be10.root with one TTree per level: thTrue [deg], keTrue [MeV].

#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"

#include <cmath>
#include <cstdio>

void make_kin_points_Be10(TString root = "/mnt/f/Be10_tp", TString outDir = "plots")
{
   gSystem->Load("libAtSimulationData.so");
   gSystem->mkdir(outDir, kTRUE);
   const int NL = 4;
   const char *LVN[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
   const double MP = 1.007825 * 931.49401;

   TFile *fo = TFile::Open(outDir + "/kin_truth_Be10.root", "RECREATE");
   for (int l = 0; l < NL; ++l) {
      TString sf = gSystem->GetFromPipe(
         TString::Format("ls %s/sims_b285/%s_s*_sim.root 2>/dev/null | head -1", root.Data(), LVN[l]));
      sf = sf.Strip(TString::kBoth);
      if (sf.IsNull()) { printf("\033[1;31mno sim file for %s\033[0m\n", LVN[l]); continue; }
      TFile *fs = TFile::Open(sf);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      if (!ts) { printf("\033[1;31mno cbmsim in %s\033[0m\n", sf.Data()); continue; }
      TClonesArray *mc = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      fo->cd();
      TTree *to = new TTree(Form("truth_%s", LVN[l]), Form("generated truth, %s", LVN[l]));
      double th, ke, z;
      to->Branch("th", &th);
      to->Branch("ke", &ke);
      to->Branch("z", &z);
      long n = 0;
      for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
         ts->GetEntry(i);
         if (!mc) continue;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *t = (AtMCTrack *)mc->At(k);
            if (!t || t->GetPdgCode() != 2212 || t->GetMotherId() != -1) continue;
            double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000;
            double p = std::sqrt(px * px + py * py + pz * pz);
            if (p <= 0) break;
            th = std::acos(pz / p) * 180.0 / M_PI;
            ke = std::sqrt(p * p + MP * MP) - MP;
            z = t->GetStartZ() * 10.0;
            to->Fill();
            ++n;
            break;
         }
      }
      to->Write();
      printf("  %-8s %7ld generated protons  (from %s)\n", LVN[l], n, gSystem->BaseName(sf));
      fs->Close();
      fo->cd();
   }
   fo->Close();
   printf("\nwrote %s/kin_truth_Be10.root\ncache done\n", outDir.Data());
}
