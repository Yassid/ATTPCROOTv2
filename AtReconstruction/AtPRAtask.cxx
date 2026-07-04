#include "AtPRAtask.h"

#include "AtDigiPar.h"       // for AtDigiPar
#include "AtEvent.h"         // for AtEvent
#include "AtPRA.h"           // for AtPRA
#include "AtPatternEvent.h"  // for AtPatternEvent
#include "AtTrackFinderRiemann.h" // for AtTrackFinderRiemann
#include "AtTrackFinderTC.h"      // for AtTrackFinderHC

#include <FairLogger.h>      // for LOG, FairLogger
#include <FairRootManager.h> // for FairRootManager
#include <FairRun.h>         // for FairRun
#include <FairRuntimeDb.h>   // for FairRuntimeDb

#include <TObject.h> // for TObject

#include <iostream>  // for operator<<, basic_ostream, cout, ostream
#include <memory>    // for unique_ptr<>::element_type, unique_ptr
#include <stdexcept> // for runtime_error
#include <utility>   // for move
#include <vector>    // for allocator, vector

AtPRAtask::AtPRAtask()
   : fInputBranchName("AtEventH"), fOutputBranchName("AtPatternEvent"), FairTask("AtPRAtask"),
     fPatternEventArray("AtPatternEvent", 1)
{

   LOG(debug) << "Default Constructor of AtPRAtask";
   fPar = nullptr;
   fPRAlgorithm = 0;
   kIsPersistence = kFALSE;
   fMinNumHits = 10;
   fMaxNumHits = 5000;

   fHCs = 0.3;
   fHCk = 19;
   fHCn = 2;
   fHCm = 15;
   fHCr = 2.0;
   fHCa = 0.03;
   fHCt = 4.0;
   fHCpadding = 0.0;

   kSetPrunning = kFALSE;
   fKNN = 5;
   fStdDevMulkNN = 0.0;
   fkNNDist = 10.0;
}

AtPRAtask::AtPRAtask(std::unique_ptr<AtPATTERN::AtPRA> pra) : AtPRAtask()
{
   fPRA = std::move(pra);
   fInjected = true;
}

AtPRAtask::~AtPRAtask()
{
   LOG(debug) << "Destructor of AtPRAtask";
}

void AtPRAtask::SetPersistence(Bool_t value)
{
   kIsPersistence = value;
}
void AtPRAtask::SetPRAlgorithm(Int_t value)
{
   fPRAlgorithm = value;
}

void AtPRAtask::SetParContainers()
{
   LOG(debug) << "SetParContainers of AtPRAtask";

   FairRun *run = FairRun::Instance();
   if (!run)
      LOG(fatal) << "No analysis run!";

   FairRuntimeDb *db = run->GetRuntimeDb(); // NOLINT
   if (!db)
      LOG(fatal) << "No runtime database!";

   fPar = (AtDigiPar *)db->getContainer("AtDigiPar"); // NOLINT
   if (!fPar)
      LOG(fatal) << "AtDigiPar not found!!";
}

InitStatus AtPRAtask::Init()
{
   LOG(debug) << "Initilization of AtPRAtask";

   // Legacy path: construct the algorithm from fPRAlgorithm. Skipped when an algorithm
   // was injected via the unique_ptr constructor (the caller configured it already);
   // only the runtime AtDigiPar diffusion bridge below then applies to it.
   if (!fInjected) {
   if (fPRAlgorithm == 0) {
      LOG(info) << "Using Track Finder TriplClust algorithm";

      fPRA.reset(new AtPATTERN::AtTrackFinderTC());
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetTcluster(fHCt);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetScluster(fHCs);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetKtriplet(fHCk);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetNtriplet(fHCn);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetMcluster(fHCm);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetRsmooth(fHCr);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetAtriplet(fHCa);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetPadding(fHCpadding);

      std::cout << " Track Finder TriplClust parameters (see Dalitz et al.) "
                << "\n";
      std::cout << " T Cluster : " << fHCt << "\n";
      std::cout << " S Cluster : " << fHCs << "\n";
      std::cout << " K Triplet : " << fHCk << "\n";
      std::cout << " N Triplet : " << fHCn << "\n";
      std::cout << " M Cluster : " << fHCm << "\n";
      std::cout << " R Smooth  : " << fHCr << "\n";
      std::cout << " A Triplet : " << fHCa << "\n";

      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetClusterRadius(fClusterRadius);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetClusterDistance(fClusterDistance);
      dynamic_cast<AtPATTERN::AtTrackFinderTC *>(fPRA.get())->SetUseSelectAndMerge(fTCUseSelectAndMerge);
      fPRA->SetChargeFromCenter(fChargeFromCenter); // base setter (robust charge sign on shallow arcs)

      std::cout << " Track finder - Parameters for clusterization "
                << "\n";
      std::cout << " Cluster radius " << fClusterRadius << "\n";
      std::cout << " Cluster distance " << fClusterDistance << "\n";

   } else if (fPRAlgorithm == 1) {
      LOG(info) << "Using RANSAC algorithm";

   } else if (fPRAlgorithm == 2) {
      LOG(info) << "Using Hough transform algorithm";
      // fPSA = new AtPSAProto();
   } else if (fPRAlgorithm == 3) {
      LOG(info) << "Using Riemann-paraboloid track finder";

      auto *riemann = new AtPATTERN::AtTrackFinderRiemann();
      riemann->SetInlierDist(fRiemannInlierDist);
      riemann->SetMinHitsPerTrack(fRiemannMinHits);
      riemann->SetMaxTracks(fRiemannMaxTracks);
      riemann->SetMaxIterations(fRiemannMaxIter);
      riemann->SetClusterRadius(fClusterRadius);
      riemann->SetClusterDistance(fClusterDistance);
      riemann->SetUseSelectAndMerge(fRiemannUseSelectAndMerge);
      riemann->SetUseArcWalkExtend(fRiemannUseArcWalkExtend);
      riemann->SetArcWalkWindow(fRiemannArcWalkWindow);
      riemann->SetArcWalkMaxMiss(fRiemannArcWalkMaxMiss);
      fPRA.reset(riemann);

      std::cout << " Riemann track finder parameters\n";
      std::cout << "   Inlier distance : " << fRiemannInlierDist << " mm\n";
      std::cout << "   Min hits/track  : " << fRiemannMinHits << "\n";
      std::cout << "   Max tracks      : " << fRiemannMaxTracks << "\n";
      std::cout << "   RANSAC iters    : " << fRiemannMaxIter << "\n";
      std::cout << "   SelectAndMerge  : " << (fRiemannUseSelectAndMerge ? "on" : "off") << "\n";
      std::cout << "   ArcWalk extend  : " << (fRiemannUseArcWalkExtend ? "on" : "off")
                << "  win=" << fRiemannArcWalkWindow
                << "  maxMiss=" << fRiemannArcWalkMaxMiss << "\n";
   }
   } // end legacy construction (skipped when an algorithm was injected)

   // Arc-walk clustering options
   if (fUseArcWalk && fPRA) {
      fPRA->SetUseArcWalk(true);
      fPRA->SetTargetClusters(fTargetClusters);
      fPRA->SetMinHitsPerCluster(fMinHitsPerCluster);
      fPRA->SetArcWalkKNN(fArcWalkKNN);
      std::cout << " Arc-walk clustering enabled: target " << fTargetClusters << " clusters/track, min "
                << fMinHitsPerCluster << " hits/cluster, kNN=" << fArcWalkKNN << "\n";
   }

   // Prunning options
   std::cout << " Track prunning : " << kSetPrunning << "\n";
   if (kSetPrunning) {
      fPRA->SetPrunning();
      std::cout << " Number of k-nearest neighbors (kNN) : " << fKNN << "\n";
      fPRA->SetkNN(fKNN);
      std::cout << " Std deviation multiplier : " << fStdDevMulkNN << "\n";
      fPRA->SetStdDevMulkNN(fStdDevMulkNN);
      std::cout << " kNN Distance threshold : " << fkNNDist << "\n";
      fPRA->SetkNNDist(fkNNDist);
   }

   // Pass diffusion parameters from AtDigiPar to the track transformer for covariance calculation
   if (fPar && fPRA) {
      double tbTime = fPar->GetTBTime() * 1e-3; // ns → us
      fPRA->SetDiffusionParams(fPar->GetCoefDiffusionTrans(), fPar->GetCoefDiffusionLong(),
                               fPar->GetDriftVelocity(), tbTime, fPadResXYOverride);
      LOG(info) << "AtPRAtask: diffusion params from AtDigiPar — CoefT=" << fPar->GetCoefDiffusionTrans()
                << " CoefL=" << fPar->GetCoefDiffusionLong() << " DriftVel=" << fPar->GetDriftVelocity()
                << " TBTime=" << tbTime << " us"
                << (fPadResXYOverride > 0 ? Form(" PadResXY=%.3g mm (override)", fPadResXYOverride) : "");
   }
   if (fUseDriftAwareWeights && fPRA) {
      fPRA->SetUseDriftAwareWeights(true);
      LOG(info) << "AtPRAtask: drift-aware σ_xy² weighting enabled for circle radius fit";
   }
   if (fPreclusterRadiusFit && fPRA) {
      fPRA->SetPreclusterRadiusFit(true);
      fPRA->SetPreclusterBin(fPreclusterBin_mm);
      LOG(info) << "AtPRAtask: PRA radius-fit pre-clustering on, bin = " << fPreclusterBin_mm << " mm";
   }

   // Get a handle from the IO manager
   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "Cannot find RootManager!";
      return kERROR;
   }

   fEventHArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranchName));
   if (fEventHArray == nullptr) {
      LOG(error) << "Cannot find AtEvent array!";
      return kERROR;
   }

   ioMan->Register(fOutputBranchName, "AtTPC", &fPatternEventArray, kIsPersistence);

   return kSUCCESS;
}

void AtPRAtask::Exec(Option_t *option)
{
   LOG(debug) << "Exec of AtPRAtask";

   fPatternEventArray.Delete();

   if (fEventHArray->GetEntriesFast() == 0)
      return;

   AtEvent &event = *(dynamic_cast<AtEvent *>(fEventHArray->At(0))); // TODO: Make sure we are not copying
   auto &hitArray = event.GetHits();

   std::cout << "  -I- AtPRAtask -  Event Number :  " << event.GetEventID() << "\n";

   try {

      if (hitArray.size() > fMinNumHits && hitArray.size() < fMaxNumHits) {
         auto patternEvent = fPRA->FindTracks(event);
         new (fPatternEventArray[0]) AtPatternEvent(std::move(*patternEvent));
      }

   } catch (std::runtime_error e) {
      std::cout << "Analyzation failed! Error: " << e.what() << std::endl;
   }

   // fEvent  = (AtEvent *) fEventHArray -> At(0);
}

void AtPRAtask::Finish()
{
   LOG(debug) << "Finish of AtPRAtask";
}
