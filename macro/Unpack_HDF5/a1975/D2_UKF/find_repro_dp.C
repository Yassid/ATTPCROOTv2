/// @file find_repro_dp.C
/// @brief Fit ONE track at a ladder of prefixes, INCLUDING below FIND's 5 % floor, with a fresh
///        fitter each time and then with one reused fitter, twice.
///
/// TWO HYPOTHESES AT ONCE, on run_0016 e1808 (984 clusters), which gf_dp_find805 stores as
/// ndf = -3, chi2 = 0, one smoothed position -- and which converges when refitted in the GUI.
///
///   Yassid's: THE GOOD FIT IS BELOW 5 %.  prod_dp_findscan.sh runs MINPCT=5 and the GUI's own
///   ladder also stops at 5 %, but in the GUI a box can be dragged to select any sub-range, so a
///   manual selection reaches places FIND cannot.  3/2/1 % of 984 clusters is 30/20/10 clusters,
///   all above the 8-cluster guard: legal fits the production never tries.
///
///   Mine: FITTER STATE.  The production runs ~16 fits back to back and then refits the winner as
///   the 17th, and AtGenfitter.cxx:648-650 never checks that refit.  If it is state, the same
///   prefix fits in a fresh fitter and collapses in a reused one.
///
/// Running both ways separates them.  Written as ONE function with no namespace and no helper
/// returning a unique_ptr: that structure made the cling parser assert outright ("backtrack pos
/// points inside the annotated tokens") and the macro never ran at all.
///
///   root -l -b -q 'find_repro_dp.C("run_0016_multifit",1808,0)'
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>

void find_repro_dp(TString runTag = "run_0016_multifit", Long64_t entry = 1808, int trackID = 0,
                   TString gfDir = "/mnt/f/a1975/gf_dp_find805/",
                   TString geoName = "ATTPC_D300torr_v2_geomanager.root",
                   TString eLossTable = "proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf = TFile::Open(dir + "/geometry/" + geoName); gf->Get("FAIRGeom"); }

   TFile *f = TFile::Open(gfDir + runTag + "_genfitter_p.root");
   if (!f || f->IsZombie()) { printf("cannot open the production file\n"); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   t->GetEntry(entry);
   AtTrackingEvent *ev = (AtTrackingEvent *)te->At(0);
   if (!ev) { printf("no event at entry %lld\n", entry); return; }
   AtTrack *src = nullptr;
   for (auto &tr : ev->GetTrackArray())
      if (tr.GetTrackID() == trackID) src = &tr;
   if (!src) { printf("no trackID %d\n", trackID); return; }
   AtTrack track0 = *src;
   const int ncl = (int)track0.GetHitClusterArray()->size();
   printf("\n=== %s entry %lld track %d : %d clusters, GeoTheta %.2f deg ===\n",
          runTag.Data(), entry, trackID, ncl, track0.GetGeoTheta() * TMath::RadToDeg());

   TString elossPath =
      eLossTable.Length() ? dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable : TString("");

   std::vector<double> keF(11, -999), ndfF(11, 0);
   const int rungs[11] = {100, 80, 60, 40, 20, 10, 5, 3, 2, 1, 0};
   const int nr = 10;   // the trailing 0 is a terminator, not a rung

   // ONE fitter, built once and reused, exactly as fit_one_dp.C does.  Constructing and DELETING
   // an AtGenfitter per rung segfaults -- tearing it down takes genfit/TGeo state with it -- so
   // the fresh-fitter arm of this test is not possible this way and is dropped.  That means this
   // macro tests ONE hypothesis: whether a good fit exists BELOW FIND's 5 % floor.
   EventFit::AtGenfitter *fit =
      new EventFit::AtGenfitter(-2.85, 2212, 1.00782503207, 1, std::string(elossPath.Data()));
   fit->SetCatimaMaterial(kTRUE, kTRUE);
   fit->SetCatimaELoss(kTRUE, kFALSE);
   if (elossPath.Length()) fit->SetELossHybrid(kTRUE, 6.61e-5);
   fit->SetZPadPlane(1000.0);
   fit->SetMeasSigma(4.0);
   fit->SetSeedFromSpyral(kFALSE);
   fit->SetMatEffectsFallback(kFALSE);
   fit->SetThetaWindow(10.0, 170.0);
   fit->SetBackwardSeedFix(kTRUE);
   fit->SetBackExtrapToAxis(kTRUE);
   fit->SetUseClusterOrder(kFALSE);
   fit->Init();

   printf("\n   %4s %6s %8s %10s %6s %9s %8s %8s\n", "pct", "nUsed", "ndf", "chi2/ndf", "pts",
          "KE", "theta", "vtx|xy|");
   for (int i = 0; i < nr; ++i) {
      const int pct = rungs[i];
      const int nUse = (int)std::lround(ncl * pct / 100.0);
      if (nUse < 8) { printf("   %4d %6d   (below the 8-cluster guard, not attempted)\n", pct, nUse); continue; }
      fit->SetTruncatePercent(pct >= 100 ? 0 : pct);      // 0 disables truncation
      AtTrack track = track0;
      AtPatternEvent pe;
      pe.AddTrack(track);
      AtTrackingEvent tev;
      fit->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
      if (tev.GetFittedTracks().empty()) { printf("   %4d %6d   -- no fitted track returned --\n", pct, nUse); continue; }
      AtFittedTrack *ft = tev.GetFittedTracks().front().get();
      if (!ft) { printf("   %4d %6d   -- null fitted track --\n", pct, nUse); continue; }
      double ndf = ft->GetTrackMetadata()->GetNdf();
      double chi2 = ft->GetTrackMetadata()->GetChi2();
      auto k = ft->GetKinematics();
      auto v = ft->GetVertex();
      char c2s[32];
      if (ndf > 0) snprintf(c2s, sizeof(c2s), "%.2f", chi2 / ndf);
      else snprintf(c2s, sizeof(c2s), "COLLAPSED");
      printf("   %4d %6d %8.0f %10s %6d %9.3f %8.2f %8.1f\n", pct, nUse, ndf, c2s,
             (int)ft->GetSmoothedPositions().size(), k.kineticEnergy,
             k.theta * TMath::RadToDeg(), std::hypot(v.X(), v.Y()));
      keF[i] = k.kineticEnergy; ndfF[i] = ndf;
   }

   printf("\n   FIND's production floor is 5 %%.  A good fit at 3, 2 or 1 %% is one the ladder\n"
          "   cannot reach, and is what a hand-drawn selection in the GUI can still find.\n");
}
