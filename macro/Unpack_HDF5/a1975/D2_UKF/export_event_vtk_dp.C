/// @file export_event_vtk_dp.C
/// @brief Dump ONE (d,p) event, every reconstruction stage, to JSON for the VTK animation
/// (vtk_event_anim.py): all hits (AtEventH), the pattern-recognition tracks and noise
/// (AtPatternEvent), the hit clusters of every track, and the genfit result of the chosen track.
///
/// Reco and genfit files are entry-aligned (the genfitter reads the reco tree sequentially);
/// this is CHECKED, not assumed: the chosen track's hit count must agree between the two files,
/// otherwise the macro stops.
///
/// Frames: hits and clusters are in the TPC frame. Genfit smoothed positions and the vertex are in
/// the FITTER frame with z mirrored -- z = 1000 - z_fit (reference_genfit_smoothed_positions_frame).
///
///   root -l -b -q 'export_event_vtk_dp.C("0020", 18249, 0, "ev.json")'
#include <fstream>
void export_event_vtk_dp(TString run = "0020", Long64_t entry = 18249, int tid = 0, TString out = "event.json",
                         TString recoDir = "/mnt/f/a1975/reco_d2_dv1104/",
                         TString gfDir = "/mnt/f/a1975/gf_dp_find805/",
                         TString gfFullDir = "/mnt/f/a1975/gf_dp_cateloss/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   const double kZPadFit = 1000.0;

   auto *fr = TFile::Open(recoDir + "run_" + run + "_multifit_reco.root");
   auto *fg = TFile::Open(gfDir + "run_" + run + "_multifit_genfitter_p.root");
   if (!fr || fr->IsZombie() || !fg || fg->IsZombie()) { printf("  cannot open reco/genfit\n"); return; }
   auto *tr = (TTree*)fr->Get("cbmsim"); auto *tg = (TTree*)fg->Get("cbmsim");
   TClonesArray *aEv = nullptr, *aPat = nullptr, *aTrk = nullptr, *aPid = nullptr;
   tr->SetBranchAddress("AtEventH", &aEv); tr->SetBranchAddress("AtPatternEvent", &aPat);
   tg->SetBranchAddress("AtTrackingEvent", &aTrk); tg->SetBranchAddress("AtPIDEvent", &aPid);
   tr->GetEntry(entry); tg->GetEntry(entry);
   auto *ev = (AtEvent*)aEv->At(0); auto *pat = (AtPatternEvent*)aPat->At(0);
   auto *tev = (AtTrackingEvent*)aTrk->At(0);
   if (!ev || !pat || !tev) { printf("  missing event objects at entry %lld\n", entry); return; }

   auto trkArr = tev->GetTrackArray();   // BY VALUE: keep the copy alive
   AtTrack *gtrk = nullptr;
   for (auto &t : trkArr) if (t.GetTrackID() == tid) gtrk = &t;
   AtFittedTrack *ft = nullptr;
   for (auto &f : tev->GetFittedTracks()) if (f && f->GetTrackID() == tid) ft = f.get();
   const AtTrack *ptrk = nullptr;
   for (auto &t : pat->GetTrackCand()) if (t.GetTrackID() == tid) ptrk = &t;
   if (!gtrk || !ft || !ptrk) { printf("  track %d not found in both files\n", tid); return; }
   if (gtrk->GetHitArray().size() != ptrk->GetHitArray().size()) {
      printf("  ALIGNMENT FAILED: genfit track has %zu hits, reco track %zu\n", gtrk->GetHitArray().size(),
             ptrk->GetHitArray().size());
      return;
   }

   std::ofstream o(out.Data());
   auto P = [&](const ROOT::Math::XYZPoint &p, double q) {
      return std::string(TString::Format("[%.1f,%.1f,%.1f,%.0f]", p.X(), p.Y(), p.Z(), q).Data()); };
   // Spyral rigidity of the same track (T m): an independent, fit-free-ish momentum referee
   double brho = -1;
   if (aPid && aPid->GetEntries() > 0)
      for (auto &sr : ((AtPIDEvent*)aPid->At(0))->GetSpyral()) if (sr.valid && sr.trackID == tid) brho = sr.brho;
   o << "{\"brho\":" << brho << ",\"run\":\"" << run << "\",\"entry\":" << entry << ",\"tid\":" << tid << ",\"hits\":[";
   bool first = true;
   for (auto &h : ev->GetHits()) { o << (first ? "" : ",") << P(h->GetPosition(), h->GetCharge()); first = false; }
   o << "],\"noise\":[";
   first = true;
   for (auto &h : pat->GetNoiseHits()) { o << (first ? "" : ",") << P(h->GetPosition(), h->GetCharge()); first = false; }
   o << "],\"tracks\":[";
   bool ft1 = true;
   for (auto &t : trkArr) {
      o << (ft1 ? "" : ",") << "{\"id\":" << t.GetTrackID() << ",\"hits\":[";
      first = true;
      for (auto &h : t.GetHitArray()) { o << (first ? "" : ",") << P(h->GetPosition(), h->GetCharge()); first = false; }
      o << "],\"clusters\":[";
      first = true;
      for (auto &c : *t.GetHitClusterArray()) { o << (first ? "" : ",") << P(c.GetPosition(), c.GetCharge()); first = false; }
      o << "]}";
      ft1 = false;
   }
   o << "],";
   auto fitJson = [&](AtFittedTrack *f, const char *key) {
      auto &k = f->GetKinematicsXtr();
      double ndf = f->GetTrackMetadata()->GetNdf();
      auto v = f->GetVertex();
      o << "\"" << key << "\":{\"ke\":" << k.kineticEnergy << ",\"theta\":" << k.theta * TMath::RadToDeg()
        << ",\"phi\":" << k.phi * TMath::RadToDeg() << ",\"chi2ndf\":"
        << (ndf > 0 ? f->GetTrackMetadata()->GetChi2() / ndf : -1) << ",\"rvtx\":" << std::hypot(v.X(), v.Y())
        << ",\"vertex\":"
        << P(ROOT::Math::XYZPoint(v.X(), v.Y(), kZPadFit - v.Z()), 0) << ",\"smoothed\":[";
      bool fs = true;
      for (auto &p : f->GetSmoothedPositions()) {
         o << (fs ? "" : ",") << P(ROOT::Math::XYZPoint(p.X(), p.Y(), kZPadFit - p.Z()), 0); fs = false; }
      o << "]}";
   };
   fitJson(ft, "fit");
   // the full-length production fit of the same track, when that production has it
   auto *ff = TFile::Open(gfFullDir + "run_" + run + "_multifit_genfitter_p.root");
   if (ff && !ff->IsZombie()) {
      auto *tf = (TTree*)ff->Get("cbmsim"); TClonesArray *af = nullptr;
      tf->SetBranchAddress("AtTrackingEvent", &af); tf->GetEntry(entry);
      auto *fev = (AtTrackingEvent*)af->At(0);
      if (fev) for (auto &f : fev->GetFittedTracks()) if (f && f->GetTrackID() == tid) { o << ","; fitJson(f.get(), "fitFull"); }
   }
   o << "}";
   o.close();
   printf("  wrote %s : %d hits, %zu noise, %zu tracks, chosen track %zu clusters, fit %zu points\n", out.Data(),
          ev->GetNumHits(), pat->GetNoiseHits().size(), trkArr.size(), gtrk->GetHitClusterArray()->size(),
          ft->GetSmoothedPositions().size());
}
