// Dual-fitter event display for a1975 16C+p: overlays BOTH Kalman filters'
// fitted trajectories on the reconstructed hits.
//
//   main 3D view  : AtEventCorrected hits (AtTabMain)
//   red  polyline : genfit  (AtTrackingEventGenfit, from AtGenfitter)
//   blue polyline : UKF     (AtTrackingEventUKF,    from AtFitterUKF)
//
// The hits + AtPatternEvent come from <run>_reco.root; the two fit branches come
// from <run>_both.root (produced by fitBoth_a1975.C), added as a TTree friend so
// both live on the same entry. The two files must be entry-aligned, i.e. run
// fitBoth_a1975.C over the SAME run with nEvents covering the events you browse.
//
// Both fitters store smoothed positions in the lab frame (z_lab = ZPadPlane - z_digi);
// SetZPadPlane(1000) maps them back to the digi frame so they overlay the hits.
//
// Usage:
//   root -l 'run_eve_both_a1975.C("run_0106", "/mnt/f/a1975/reco/", "/tmp/")'
//     arg1 = run name, arg2 = dir with _reco.root, arg3 = dir with _both.root
void run_eve_both_a1975(TString runName = "run_0106", TString recoDir = "/mnt/f/a1975/reco/",
                        TString fitDir = "/tmp/")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtEventDisplay.so");

   TString recoFile = recoDir + runName + "_reco.root"; // AtEventCorrected + AtPatternEvent
   TString fitFile = fitDir + runName + "_both.root";   // AtTrackingEventGenfit + AtTrackingEventUKF
   if (gSystem->AccessPathName(recoFile)) {
      std::cout << "\033[1;31mMissing " << recoFile << "\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(fitFile)) {
      std::cout << "\033[1;31mMissing " << fitFile << " — run fitBoth_a1975.C first.\033[0m\n";
      return;
   }
   std::cout << "hits/pattern : " << recoFile << "\nfits (friend): " << fitFile << "\n";

   TString dir = getenv("VMCWORKDIR");
   TString geoFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   TString mapDir = dir + "/scripts/ANL2023.xml";

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(recoFile);
   source->AddFriend(fitFile); // merge the two fit branches onto the same entries
   fRun->SetSource(source);
   fRun->SetSink(new FairRootFileSink(fitDir + runName + "_both.reco_display.root"));
   fRun->SetGeomFile(geoFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   rtdb->setFirstInput(new FairParRootFileIo());

   auto fMap = std::make_shared<AtTpcMap>();
   fMap->ParseXMLMap(mapDir.Data());
   auto eveMan = new AtViewerManager(fMap);

   // 3D hits (select the "AtEventCorrected" branch in the sidebar if not default)
   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100);

   // genfit fitted trajectory — red, solid
   auto tabGenfit = std::make_unique<AtTabFitted>("Genfit", "AtTrackingEventGenfit");
   tabGenfit->SetDrawSmoothed(true);
   tabGenfit->SetZPadPlane(1000.0);
   tabGenfit->SetTrackColor(kRed + 1);
   tabGenfit->SetLineWidth(2);

   // UKF fitted trajectory — blue, dashed (so overlapping fits stay distinguishable)
   auto tabUKF = std::make_unique<AtTabFitted>("UKF", "AtTrackingEventUKF");
   tabUKF->SetDrawSmoothed(true);
   tabUKF->SetZPadPlane(1000.0);
   tabUKF->SetTrackColor(kAzure + 1);
   tabUKF->SetLineStyle(2);
   tabUKF->SetLineWidth(2);

   eveMan->AddTab(std::move(tabMain));
   eveMan->AddTab(std::move(tabGenfit));
   eveMan->AddTab(std::move(tabUKF));

   eveMan->Init();
   std::cout << "Finished init — red = genfit, blue dashed = UKF, over the hits.\n";
}
