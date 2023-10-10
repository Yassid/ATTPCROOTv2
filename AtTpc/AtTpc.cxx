/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
#include "AtTpc.h"

#include "AtDetectorList.h"
#include "AtMCPoint.h"
#include "AtStack.h"
#include "AtVertexPropagator.h"

#include <FairDetector.h>
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRun.h>
#include <FairRuntimeDb.h>
#include <FairVolume.h>

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TVirtualMC.h>
#include <TVirtualMCStack.h>

#include <iostream>

using std::cout;
using std::endl;

constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";
constexpr auto cBLUE = "\033[1;34m";

AtTpc::AtTpc() : AtTpc("AtTpc", true) {}

AtTpc::AtTpc(const char *name, Bool_t active)
   : FairDetector(name, active, kAtTpc), fTrackID(-1), fVolumeID(-1), fPos(), fMom(), fTime(-1.), fLength(-1.),
     fELoss(-1), fPosIndex(-1), fAtTpcPointCollection(new TClonesArray("AtMCPoint")), fELossAcc(-1)
{
}

AtTpc::~AtTpc()
{
   if (fAtTpcPointCollection) {
      fAtTpcPointCollection->Delete();
      delete fAtTpcPointCollection;
   }
}

void AtTpc::Initialize()
{
   FairDetector::Initialize();
   FairRuntimeDb *rtdb = FairRun::Instance()->GetRuntimeDb();
   rtdb->getContainer("AtTpcGeoPar");
}

void AtTpc::trackEnteringVolume()
{
   auto AZ = DecodePdG(gMC->TrackPid());
   fELoss = 0.;
   fELossAcc = 0.;
   fTime = gMC->TrackTime() * 1.0e09;
   fLength = gMC->TrackLength();
   gMC->TrackPosition(fPosIn);
   gMC->TrackMomentum(fMomIn);
   fTrackID = gMC->GetStack()->GetCurrentTrackNumber();

   // Position of the first hit of the beam in the TPC volume ( For tracking purposes in the TPC)
   if (AtVertexPropagator::Instance()->GetBeamEvtCnt() % 2 != 0 && fTrackID == 0 &&
       (fVolName == "drift_volume" || fVolName == "cell"))
      InPos = fPosIn;

   Int_t VolumeID = 0;

   if (AtVertexPropagator::Instance()->GetBeamEvtCnt() % 2 != 0)
      LOG(debug) << cGREEN << " AtTPC: Beam Event ";
   else if (AtVertexPropagator::Instance()->GetDecayEvtCnt() % 2 == 0)
      LOG(debug) << cBLUE << " AtTPC: Reaction/Decay Event ";

   LOG(debug) << " AtTPC: First hit in Volume " << fVolName;
   LOG(debug) << " Particle : " << gMC->ParticleName(gMC->TrackPid());
   LOG(debug) << " PID PdG : " << gMC->TrackPid();
   LOG(debug) << " Atomic Mass : " << AZ.first;
   LOG(debug) << " Atomic Number : " << AZ.second;
   LOG(debug) << " Volume ID " << gMC->CurrentVolID(VolumeID);
   LOG(debug) << " Track ID : " << fTrackID;
   LOG(debug) << " Position : " << fPosIn.X() << " " << fPosIn.Y() << "  " << fPosIn.Z();
   LOG(debug) << " Momentum : " << fMomIn.X() << " " << fMomIn.Y() << "  " << fMomIn.Z();
   LOG(debug) << " Total relativistic energy " << gMC->Etot();
   LOG(debug) << " Mass of the Beam particle (gAVTP) : " << AtVertexPropagator::Instance()->GetBeamMass();
   LOG(debug) << " Mass of the Tracked particle (gMC) : " << gMC->TrackMass(); // NB: with electrons
   LOG(debug) << " Initial energy of the beam particle in this volume : "
              << ((gMC->Etot() - AtVertexPropagator::Instance()->GetBeamMass() * 0.93149401) *
                  1000.); // Relativistic Mass
   LOG(debug) << " Total energy of the current track (gMC) : "
              << ((gMC->Etot() - gMC->TrackMass()) * 1000.); // Relativistic Mass
   LOG(debug) << " ==================================================== " << cNORMAL;
}

void AtTpc::getTrackParametersFromMC()
{
   fELoss = gMC->Edep();
   fELossAcc += fELoss;
   fTime = gMC->TrackTime() * 1.0e09;
   fLength = gMC->TrackLength();
   gMC->TrackPosition(fPosIn);
   gMC->TrackMomentum(fMomIn);
   fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
}

void AtTpc::getTrackParametersWhileExiting()
{
   fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
   gMC->TrackPosition(fPosOut);
   gMC->TrackMomentum(fMomOut);

   // Correct fPosOut
   if (gMC->IsTrackExiting()) {
      correctPosOut();
      if ((fVolName.Contains("drift_volume") || fVolName.Contains("cell")) &&
          AtVertexPropagator::Instance()->GetBeamEvtCnt() % 2 != 0 && fTrackID == 0)
         resetVertex();
   }
}

void AtTpc::resetVertex()
{
   AtVertexPropagator::Instance()->ResetVertex();
   LOG(INFO) << cRED << " - AtTpc Warning : Beam punched through the AtTPC. Reseting Vertex! " << cNORMAL << std::endl;
}

void AtTpc::correctPosOut()
{
   const Double_t *oldpos = nullptr;
   const Double_t *olddirection = nullptr;
   Double_t newpos[3];
   Double_t newdirection[3];
   Double_t safety = 0;

   gGeoManager->FindNode(fPosOut.X(), fPosOut.Y(), fPosOut.Z());
   oldpos = gGeoManager->GetCurrentPoint();
   olddirection = gGeoManager->GetCurrentDirection();

   for (Int_t i = 0; i < 3; i++) {
      newdirection[i] = -1 * olddirection[i];
   }

   gGeoManager->SetCurrentDirection(newdirection);
   safety = gGeoManager->GetSafeDistance(); // Get distance to boundry?
   gGeoManager->SetCurrentDirection(-newdirection[0], -newdirection[1], -newdirection[2]);

   for (Int_t i = 0; i < 3; i++) {
      newpos[i] = oldpos[i] - (3 * safety * olddirection[i]);
   }

   fPosOut.SetX(newpos[0]);
   fPosOut.SetY(newpos[1]);
   fPosOut.SetZ(newpos[2]);
}

bool AtTpc::reactionOccursHere()
{
   bool atEnergyLoss = fELossAcc * 1000 > AtVertexPropagator::Instance()->GetRndELoss();
   bool isPrimaryBeam = AtVertexPropagator::Instance()->GetBeamEvtCnt() % 2 != 0 && fTrackID == 0;
   bool isInRightVolume = fVolName.Contains("drift_volume") || fVolName.Contains("cell");
   return atEnergyLoss && isPrimaryBeam && isInRightVolume;
}
Bool_t AtTpc::ProcessHits(FairVolume *vol)
{
   /** This method is called from the MC stepping */

   auto *stack = dynamic_cast<AtStack *>(gMC->GetStack());
   fVolName = gMC->CurrentVolName();
   fVolumeID = vol->getMCid();
   fDetCopyID = vol->getCopyNo();

   if (gMC->IsTrackEntering())
      trackEnteringVolume();

   getTrackParametersFromMC();

   if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
      getTrackParametersWhileExiting();

   addHit();

   // Reaction Occurs here
   if (reactionOccursHere())
      startReactionEvent();

   // Increment number of AtTpc det points in TParticle
   stack->AddPoint(kAtTpc);
   return kTRUE;
}

void AtTpc::startReactionEvent()
{

   gMC->StopTrack();
   AtVertexPropagator::Instance()->ResetVertex();

   TLorentzVector StopPos;
   TLorentzVector StopMom;
   gMC->TrackPosition(StopPos);
   gMC->TrackMomentum(StopMom);
   Double_t StopEnergy = ((gMC->Etot() - AtVertexPropagator::Instance()->GetBeamMass() * 0.93149401) * 1000.);

   LOG(debug) << cYELLOW << " Beam energy loss before reaction : " << fELossAcc * 1000;
   LOG(debug) << " Mass of the Tracked particle : " << gMC->TrackMass();
   LOG(debug) << " Mass of the Beam particle (gAVTP)  : " << AtVertexPropagator::Instance()->GetBeamMass();
   LOG(debug) << " Total energy of the Beam particle before reaction : " << StopEnergy << cNORMAL; // Relativistic Mass

   AtVertexPropagator::Instance()->SetVertex(StopPos.X(), StopPos.Y(), StopPos.Z(), InPos.X(), InPos.Y(), InPos.Z(),
                                             StopMom.Px(), StopMom.Py(), StopMom.Pz(), StopEnergy);
}

void AtTpc::addHit()
{
   auto AZ = DecodePdG(gMC->TrackPid());

   Double_t EIni = 0;
   Double_t AIni = 0;

   // We assume that the beam-like particle is fTrackID==0 since it is the first one added in the
   // primary generator
   if (AtVertexPropagator::Instance()->GetBeamEvtCnt() % 2 != 0 && fTrackID == 0) {
      EIni = 0;
      AIni = 0;
   } else if (AtVertexPropagator::Instance()->GetDecayEvtCnt() % 2 == 0) {
      EIni = AtVertexPropagator::Instance()->GetTrackEnergy(fTrackID);
      AIni = AtVertexPropagator::Instance()->GetTrackAngle(fTrackID);
   }

   AddHit(fTrackID, fVolumeID, fVolName, fDetCopyID, TVector3(fPosIn.X(), fPosIn.Y(), fPosIn.Z()),
          TVector3(fMomIn.Px(), fMomIn.Py(), fMomIn.Pz()), fTime, fLength, fELoss, EIni, AIni, AZ.first, AZ.second);
}

void AtTpc::EndOfEvent()
{

   fAtTpcPointCollection->Clear();
}

void AtTpc::Register()
{
   FairRootManager::Instance()->Register("AtTpcPoint", "AtTpc", fAtTpcPointCollection, kTRUE);
}

TClonesArray *AtTpc::GetCollection(Int_t iColl) const
{
   if (iColl == 0) {
      return fAtTpcPointCollection;
   } else {
      return nullptr;
   }
}

void AtTpc::Reset()
{
   fAtTpcPointCollection->Clear();
}

void AtTpc::Print(Option_t *option) const
{
   Int_t nHits = fAtTpcPointCollection->GetEntriesFast();
   LOG(INFO) << "AtTPC: " << nHits << " points registered in this event";
}

void AtTpc::ConstructGeometry()
{
   TString fileName = GetGeometryFileName();
   if (fileName.EndsWith(".geo")) {
      LOG(INFO) << "Constructing AtTPC geometry from ASCII file " << fileName;
      // ConstructASCIIGeometry();
   } else if (fileName.EndsWith(".root")) {
      LOG(INFO) << "Constructing AtTPC geometry from ROOT file " << fileName;
      ConstructRootGeometry();

   } else {
      std::cout << "Geometry format not supported." << std::endl;
   }
}

void AtTpc::ConstructOpGeometry()
{

   TString fileName = GetGeometryFileName();
   if (!fileName.Contains("OTPC"))
      return;

   auto cf4_id = gMC->MediumId("CF4_250mbar");
   std::cout << " Material id  " << cf4_id << "\n";

   const Int_t iNbEntries = 300;
   Double_t CF4PhotonMomentum[iNbEntries] = {
      6.2e-09,         6.138613861e-09, 6.078431373e-09, 6.019417476e-09, 5.961538462e-09, 5.904761905e-09,
      5.849056604e-09, 5.794392523e-09, 5.740740741e-09, 5.688073394e-09, 5.636363636e-09, 5.585585586e-09,
      5.535714286e-09, 5.486725664e-09, 5.438596491e-09, 5.391304348e-09, 5.344827586e-09, 5.299145299e-09,
      5.254237288e-09, 5.210084034e-09, 5.166666667e-09, 5.123966942e-09, 5.081967213e-09, 5.040650407e-09,
      5e-09,           4.96e-09,        4.920634921e-09, 4.881889764e-09, 4.84375e-09,     4.80620155e-09,
      4.769230769e-09, 4.732824427e-09, 4.696969697e-09, 4.661654135e-09, 4.626865672e-09, 4.592592593e-09,
      4.558823529e-09, 4.525547445e-09, 4.492753623e-09, 4.460431655e-09, 4.428571429e-09, 4.397163121e-09,
      4.366197183e-09, 4.335664336e-09, 4.305555556e-09, 4.275862069e-09, 4.246575342e-09, 4.217687075e-09,
      4.189189189e-09, 4.161073826e-09, 4.133333333e-09, 4.105960265e-09, 4.078947368e-09, 4.052287582e-09,
      4.025974026e-09, 4e-09,           3.974358974e-09, 3.949044586e-09, 3.924050633e-09, 3.899371069e-09,
      3.875e-09,       3.850931677e-09, 3.827160494e-09, 3.803680982e-09, 3.780487805e-09, 3.757575758e-09,
      3.734939759e-09, 3.71257485e-09,  3.69047619e-09,  3.668639053e-09, 3.647058824e-09, 3.625730994e-09,
      3.604651163e-09, 3.583815029e-09, 3.563218391e-09, 3.542857143e-09, 3.522727273e-09, 3.502824859e-09,
      3.483146067e-09, 3.463687151e-09, 3.444444444e-09, 3.425414365e-09, 3.406593407e-09, 3.387978142e-09,
      3.369565217e-09, 3.351351351e-09, 3.333333333e-09, 3.315508021e-09, 3.29787234e-09,  3.28042328e-09,
      3.263157895e-09, 3.246073298e-09, 3.229166667e-09, 3.212435233e-09, 3.195876289e-09, 3.179487179e-09,
      3.163265306e-09, 3.147208122e-09, 3.131313131e-09, 3.115577889e-09, 3.1e-09,         3.084577114e-09,
      3.069306931e-09, 3.054187192e-09, 3.039215686e-09, 3.024390244e-09, 3.009708738e-09, 2.995169082e-09,
      2.980769231e-09, 2.966507177e-09, 2.952380952e-09, 2.938388626e-09, 2.924528302e-09, 2.910798122e-09,
      2.897196262e-09, 2.88372093e-09,  2.87037037e-09,  2.857142857e-09, 2.844036697e-09, 2.831050228e-09,
      2.818181818e-09, 2.805429864e-09, 2.792792793e-09, 2.780269058e-09, 2.767857143e-09, 2.755555556e-09,
      2.743362832e-09, 2.731277533e-09, 2.719298246e-09, 2.707423581e-09, 2.695652174e-09, 2.683982684e-09,
      2.672413793e-09, 2.660944206e-09, 2.64957265e-09,  2.638297872e-09, 2.627118644e-09, 2.616033755e-09,
      2.605042017e-09, 2.594142259e-09, 2.583333333e-09, 2.572614108e-09, 2.561983471e-09, 2.551440329e-09,
      2.540983607e-09, 2.530612245e-09, 2.520325203e-09, 2.510121457e-09, 2.5e-09,         2.489959839e-09,
      2.48e-09,        2.470119522e-09, 2.46031746e-09,  2.450592885e-09, 2.440944882e-09, 2.431372549e-09,
      2.421875e-09,    2.412451362e-09, 2.403100775e-09, 2.393822394e-09, 2.384615385e-09, 2.375478927e-09,
      2.366412214e-09, 2.357414449e-09, 2.348484848e-09, 2.339622642e-09, 2.330827068e-09, 2.322097378e-09,
      2.313432836e-09, 2.304832714e-09, 2.296296296e-09, 2.287822878e-09, 2.279411765e-09, 2.271062271e-09,
      2.262773723e-09, 2.254545455e-09, 2.246376812e-09, 2.238267148e-09, 2.230215827e-09, 2.222222222e-09,
      2.214285714e-09, 2.206405694e-09, 2.19858156e-09,  2.190812721e-09, 2.183098592e-09, 2.175438596e-09,
      2.167832168e-09, 2.160278746e-09, 2.152777778e-09, 2.14532872e-09,  2.137931034e-09, 2.130584192e-09,
      2.123287671e-09, 2.116040956e-09, 2.108843537e-09, 2.101694915e-09, 2.094594595e-09, 2.087542088e-09,
      2.080536913e-09, 2.073578595e-09, 2.066666667e-09, 2.059800664e-09, 2.052980132e-09, 2.04620462e-09,
      2.039473684e-09, 2.032786885e-09, 2.026143791e-09, 2.019543974e-09, 2.012987013e-09, 2.006472492e-09,
      2e-09,           1.993569132e-09, 1.987179487e-09, 1.980830671e-09, 1.974522293e-09, 1.968253968e-09,
      1.962025316e-09, 1.955835962e-09, 1.949685535e-09, 1.943573668e-09, 1.9375e-09,      1.931464174e-09,
      1.925465839e-09, 1.919504644e-09, 1.913580247e-09, 1.907692308e-09, 1.901840491e-09, 1.896024465e-09,
      1.890243902e-09, 1.88449848e-09,  1.878787879e-09, 1.873111782e-09, 1.86746988e-09,  1.861861862e-09,
      1.856287425e-09, 1.850746269e-09, 1.845238095e-09, 1.839762611e-09, 1.834319527e-09, 1.828908555e-09,
      1.823529412e-09, 1.818181818e-09, 1.812865497e-09, 1.807580175e-09, 1.802325581e-09, 1.797101449e-09,
      1.791907514e-09, 1.786743516e-09, 1.781609195e-09, 1.776504298e-09, 1.771428571e-09, 1.766381766e-09,
      1.761363636e-09, 1.756373938e-09, 1.751412429e-09, 1.746478873e-09, 1.741573034e-09, 1.736694678e-09,
      1.731843575e-09, 1.727019499e-09, 1.722222222e-09, 1.717451524e-09, 1.712707182e-09, 1.707988981e-09,
      1.703296703e-09, 1.698630137e-09, 1.693989071e-09, 1.689373297e-09, 1.684782609e-09, 1.680216802e-09,
      1.675675676e-09, 1.67115903e-09,  1.666666667e-09, 1.662198391e-09, 1.657754011e-09, 1.653333333e-09,
      1.64893617e-09,  1.644562334e-09, 1.64021164e-09,  1.635883905e-09, 1.631578947e-09, 1.627296588e-09,
      1.623036649e-09, 1.618798956e-09, 1.614583333e-09, 1.61038961e-09,  1.606217617e-09, 1.602067183e-09,
      1.597938144e-09, 1.593830334e-09, 1.58974359e-09,  1.585677749e-09, 1.581632653e-09, 1.577608142e-09,
      1.573604061e-09, 1.569620253e-09, 1.565656566e-09, 1.561712846e-09, 1.557788945e-09, 1.553884712e-09};
   Double_t CF4Scintillation_Fast[iNbEntries] = {
      0.0029, 0.0029, 0.0017, 0.0024, 0.0018, 0.0011, 0.0027, 0.0009, 0.0003, 0.0019, 0.0030, 0.0024, 0.0023, 0.0036,
      0.0039, 0.0056, 0.0049, 0.0061, 0.0053, 0.0052, 0.0056, 0.0064, 0.0072, 0.0064, 0.0080, 0.0071, 0.0056, 0.0069,
      0.0053, 0.0070, 0.0060, 0.0057, 0.0071, 0.0066, 0.0066, 0.0055, 0.0082, 0.0076, 0.0093, 0.0089, 0.0106, 0.0109,
      0.0105, 0.0102, 0.0120, 0.0121, 0.0102, 0.0097, 0.0120, 0.0126, 0.0097, 0.0103, 0.0097, 0.0084, 0.0119, 0.0112,
      0.0096, 0.0171, 0.0235, 0.0078, 0.0089, 0.0071, 0.0065, 0.0074, 0.0073, 0.0074, 0.0074, 0.0080, 0.0143, 0.0522,
      0.0069, 0.0076, 0.0042, 0.0059, 0.0039, 0.0053, 0.0054, 0.0185, 0.0077, 0.0599, 0.0048, 0.0034, 0.0041, 0.0041,
      0.0047, 0.0059, 0.0046, 0.0065, 0.0128, 0.0037, 0.0167, 0.0053, 0.0038, 0.0042, 0.0046, 0.0032, 0.0037, 0.0073,
      0.0049, 0.0067, 0.0116, 0.0054, 0.0077, 0.0111, 0.0042, 0.0043, 0.0037, 0.0046, 0.0041, 0.0028, 0.0055, 0.0031,
      0.0048, 0.0057, 0.0056, 0.0035, 0.0039, 0.0068, 0.0051, 0.0037, 0.0054, 0.0048, 0.0061, 0.0033, 0.0050, 0.0052,
      0.0047, 0.0014, 0.0043, 0.0041, 0.0023, 0.0062, 0.0036, 0.0038, 0.0039, 0.0043, 0.0049, 0.0049, 0.0036, 0.0048,
      0.0039, 0.0023, 0.0035, 0.0025, 0.0036, 0.0010, 0.0044, 0.0013, 0.0041, 0.0021, 0.0016, 0.0046, 0.0040, 0.0034,
      0.0027, 0.0026, 0.0034, 0.0004, 0.0037, 0.0004, 0.0036, 0.0029, 0.0029, 0.0036, 0.0055, 0.0034, 0.0034, 0.0025,
      0.0028, 0.0055, 0.0064, 0.0037, 0.0029, 0.0047, 0.0058, 0.0040, 0.0062, 0.0055, 0.0029, 0.0067, 0.0070, 0.0080,
      0.0060, 0.0094, 0.0082, 0.0072, 0.0089, 0.0117, 0.0102, 0.0134, 0.0131, 0.0131, 0.0120, 0.0135, 0.0096, 0.0107,
      0.0179, 0.0210, 0.0172, 0.0165, 0.0167, 0.0176, 0.0137, 0.0196, 0.0217, 0.0175, 0.0223, 0.0192, 0.0222, 0.0188,
      0.0184, 0.0183, 0.0156, 0.0098, 0.0198, 0.0268, 0.0188, 0.0236, 0.0208, 0.0171, 0.0229, 0.0228, 0.0227, 0.0204,
      0.0184, 0.0190, 0.0185, 0.0145, 0.0138, 0.0122, 0.0180, 0.0132, 0.0146, 0.0087, 0.0039, 0.0147, 0.0000, 0.0000,
      0.0137, 0.0084, 0.0094, 0.0114, 0.0078, 0.0100, 0.0069, 0.0055, 0.0164, 0.0113, 0.0148, 0.0053, 0.0054, 0.0065,
      0.0092, 0.0000, 0.0047, 0.0000, 0.0071, 0.0000, 0.0057, 0.0063, 0.0064, 0.0050, 0.0077, 0.0034, 0.0025, 0.0000,
      0.0041, 0.0025, 0.0019, 0.0042, 0.0030, 0.0000, 0.0030, 0.0000, 0.0000, 0.0000, 0.0027, 0.0000, 0.0000, 0.0000,
      0.0000, 0.0006, 0.0051, 0.0083, 0.0000, 0.0000, 0.0064, 0.0003, 0.0002, 0.0074, 0.0038, 0.0000, 0.0000, 0.0000,
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};
   Double_t CF4Scintillation_Slow[iNbEntries] = {
      0.0029, 0.0029, 0.0017, 0.0024, 0.0018, 0.0011, 0.0027, 0.0009, 0.0003, 0.0019, 0.0030, 0.0024, 0.0023, 0.0036,
      0.0039, 0.0056, 0.0049, 0.0061, 0.0053, 0.0052, 0.0056, 0.0064, 0.0072, 0.0064, 0.0080, 0.0071, 0.0056, 0.0069,
      0.0053, 0.0070, 0.0060, 0.0057, 0.0071, 0.0066, 0.0066, 0.0055, 0.0082, 0.0076, 0.0093, 0.0089, 0.0106, 0.0109,
      0.0105, 0.0102, 0.0120, 0.0121, 0.0102, 0.0097, 0.0120, 0.0126, 0.0097, 0.0103, 0.0097, 0.0084, 0.0119, 0.0112,
      0.0096, 0.0171, 0.0235, 0.0078, 0.0089, 0.0071, 0.0065, 0.0074, 0.0073, 0.0074, 0.0074, 0.0080, 0.0143, 0.0522,
      0.0069, 0.0076, 0.0042, 0.0059, 0.0039, 0.0053, 0.0054, 0.0185, 0.0077, 0.0599, 0.0048, 0.0034, 0.0041, 0.0041,
      0.0047, 0.0059, 0.0046, 0.0065, 0.0128, 0.0037, 0.0167, 0.0053, 0.0038, 0.0042, 0.0046, 0.0032, 0.0037, 0.0073,
      0.0049, 0.0067, 0.0116, 0.0054, 0.0077, 0.0111, 0.0042, 0.0043, 0.0037, 0.0046, 0.0041, 0.0028, 0.0055, 0.0031,
      0.0048, 0.0057, 0.0056, 0.0035, 0.0039, 0.0068, 0.0051, 0.0037, 0.0054, 0.0048, 0.0061, 0.0033, 0.0050, 0.0052,
      0.0047, 0.0014, 0.0043, 0.0041, 0.0023, 0.0062, 0.0036, 0.0038, 0.0039, 0.0043, 0.0049, 0.0049, 0.0036, 0.0048,
      0.0039, 0.0023, 0.0035, 0.0025, 0.0036, 0.0010, 0.0044, 0.0013, 0.0041, 0.0021, 0.0016, 0.0046, 0.0040, 0.0034,
      0.0027, 0.0026, 0.0034, 0.0004, 0.0037, 0.0004, 0.0036, 0.0029, 0.0029, 0.0036, 0.0055, 0.0034, 0.0034, 0.0025,
      0.0028, 0.0055, 0.0064, 0.0037, 0.0029, 0.0047, 0.0058, 0.0040, 0.0062, 0.0055, 0.0029, 0.0067, 0.0070, 0.0080,
      0.0060, 0.0094, 0.0082, 0.0072, 0.0089, 0.0117, 0.0102, 0.0134, 0.0131, 0.0131, 0.0120, 0.0135, 0.0096, 0.0107,
      0.0179, 0.0210, 0.0172, 0.0165, 0.0167, 0.0176, 0.0137, 0.0196, 0.0217, 0.0175, 0.0223, 0.0192, 0.0222, 0.0188,
      0.0184, 0.0183, 0.0156, 0.0098, 0.0198, 0.0268, 0.0188, 0.0236, 0.0208, 0.0171, 0.0229, 0.0228, 0.0227, 0.0204,
      0.0184, 0.0190, 0.0185, 0.0145, 0.0138, 0.0122, 0.0180, 0.0132, 0.0146, 0.0087, 0.0039, 0.0147, 0.0000, 0.0000,
      0.0137, 0.0084, 0.0094, 0.0114, 0.0078, 0.0100, 0.0069, 0.0055, 0.0164, 0.0113, 0.0148, 0.0053, 0.0054, 0.0065,
      0.0092, 0.0000, 0.0047, 0.0000, 0.0071, 0.0000, 0.0057, 0.0063, 0.0064, 0.0050, 0.0077, 0.0034, 0.0025, 0.0000,
      0.0041, 0.0025, 0.0019, 0.0042, 0.0030, 0.0000, 0.0030, 0.0000, 0.0000, 0.0000, 0.0027, 0.0000, 0.0000, 0.0000,
      0.0000, 0.0006, 0.0051, 0.0083, 0.0000, 0.0000, 0.0064, 0.0003, 0.0002, 0.0074, 0.0038, 0.0000, 0.0000, 0.0000,
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000};

   const Int_t iNbEntries_1 = 3;
   Bool_t spline = true;
   Bool_t noSpline = false;
   Double_t CF4PhotonMomentum_1[iNbEntries_1] = {6.2e-09, 6.138613861e-09, 6.078431373e-09};
   Double_t CF4RefractiveIndex[iNbEntries_1] = {1.004, 1.004, 1.004};
   Double_t CF4AbsorbtionLength[iNbEntries_1] = {10., 10., 10.};
   Double_t CF4ScatteringLength[iNbEntries_1] = {3., 3., 3.};

   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONYIELD", 1200); //
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONCOMPONENT1", iNbEntries, CF4PhotonMomentum, CF4Scintillation_Fast,
                            false, spline);
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONCOMPONENT2", iNbEntries, CF4PhotonMomentum, CF4Scintillation_Slow,
                            false, spline);
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONTIMECONSTANT1", 3.0e-09);  // 3.*ns
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONTIMECONSTANT2", 10.0e-09); // 10.*ns
   gMC->SetMaterialProperty(cf4_id, "RESOLUTIONSCALE", 1.0);
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONYIELD1", 1.0);
   gMC->SetMaterialProperty(cf4_id, "SCINTILLATIONYIELD2", 1.0);
   gMC->SetMaterialProperty(cf4_id, "RINDEX", iNbEntries_1, CF4PhotonMomentum_1, CF4ScatteringLength);

   /*G4MaterialPropertiesTable *CF4PropertiesTable = new G4MaterialPropertiesTable();
   CF4PropertiesTable->AddProperty("FASTCOMPONENT", CF4PhotonMomentum, CF4Scintillation_Fast, iNbEntries);
   CF4PropertiesTable->AddProperty("SLOWCOMPONENT", CF4PhotonMomentum, CF4Scintillation_Slow, iNbEntries);
   CF4PropertiesTable->AddProperty("RINDEX", CF4PhotonMomentum_1, CF4RefractiveIndex, iNbEntries_1);
   CF4PropertiesTable->AddProperty("ABSLENGTH", CF4PhotonMomentum_1, CF4AbsorbtionLength, iNbEntries_1);
   CF4PropertiesTable->AddProperty("RAYLEIGH", CF4PhotonMomentum_1, CF4ScatteringLength, iNbEntries_1);
   CF4PropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 2500./keV);  // for electron recoil
   CF4PropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
   CF4PropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
   CF4PropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
   CF4PropertiesTable->AddConstProperty("YIELDRATIO", 1.0);
   CF4->SetMaterialPropertiesTable(CF4PropertiesTable); */
}

Bool_t AtTpc::CheckIfSensitive(std::string name)
{

   TString tsname = name;
   if (tsname.Contains("drift_volume") || tsname.Contains("window") || tsname.Contains("cell")) {
      LOG(INFO) << " AtTPC geometry: Sensitive volume found: " << tsname;
      return kTRUE;
   }
   return kFALSE;
}

AtMCPoint *
AtTpc::AddHit(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length, Double_t eLoss)
{
   TClonesArray &clref = *fAtTpcPointCollection;
   Int_t size = clref.GetEntriesFast();
   return new (clref[size]) AtMCPoint(trackID, detID, pos, mom, time, length, eLoss);
}

// -----   Private method AddHit   --------------------------------------------
AtMCPoint *AtTpc::AddHit(Int_t trackID, Int_t detID, TString VolName, Int_t detCopyID, TVector3 pos, TVector3 mom,
                         Double_t time, Double_t length, Double_t eLoss, Double_t EIni, Double_t AIni, Int_t A, Int_t Z)
{
   TClonesArray &clref = *fAtTpcPointCollection;
   Int_t size = clref.GetEntriesFast();
   if (fVerboseLevel > 1)
      LOG(INFO) << "AtTPC: Adding Point at (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ") cm,  detector "
                << detID << ", track " << trackID << ", energy loss " << eLoss * 1e06 << " keV";

   return new (clref[size])
      AtMCPoint(trackID, detID, pos, mom, time, length, eLoss, VolName, detCopyID, EIni, AIni, A, Z);
}

std::pair<Int_t, Int_t> AtTpc::DecodePdG(Int_t PdG_Code)
{
   Int_t A = PdG_Code / 10 % 1000;
   Int_t Z = PdG_Code / 10000 % 1000;

   std::pair<Int_t, Int_t> nucleus;

   if (PdG_Code == 2212) {
      nucleus.first = 1;
      nucleus.second = 1;
   } else if (PdG_Code == 2112) {
      nucleus.first = 1;
      nucleus.second = 0;
   } else {
      nucleus.first = A;
      nucleus.second = Z;
   }

   return nucleus;
}

ClassImp(AtTpc)
