#include "AtPSA.h"

// FairRoot classes
#include "FairRuntimeDb.h"
#include "FairRun.h"

// ROOT classes
#include "TClonesArray.h"
#include "TSpectrum.h"
#include "TVector3.h"
#include "TMath.h"

// STL
#include <algorithm>

// AtTPCROOT classes
#include "AtRawEvent.h"
#include "AtEvent.h"
#include "AtDigiPar.h"
#include "AtLangevin.h"
#include "AtCalibration.h"
#include "AtHit.h"
#include "AtTpcPoint.h"

using std::distance;
using std::max_element;
using std::min_element;

AtPSA::AtPSA()
{
   std::cout << "Calling AtPSA Constructor" << std::endl;

   // TODO:Move to class that needs them
   fIniTB = 0;
   fEndTB = 512;

   fThreshold = -1;
   fThresholdlow = -1;
   fUsingLowThreshold = kFALSE;

   fIsGainCalibrated = kFALSE;
   fIsJitterCalibrated = kFALSE;
}

AtPSA::~AtPSA()
{

   delete fCalibration;
}

void AtPSA::Init()
{
   fCalibration = new AtCalibration();

   FairRun *run = FairRun::Instance();
   if (!run)
      LOG(FATAL) << "No analysis run!";

   FairRuntimeDb *db = run->GetRuntimeDb();
   if (!db)
      LOG(FATAL) << "No runtime database!";

   fPar = (AtDigiPar *)db->getContainer("AtDigiPar");
   if (!fPar)
      LOG(FATAL) << "AtDigiPar not found!!";

   fPadPlaneX = fPar->GetPadPlaneX();
   fPadSizeX = fPar->GetPadSizeX();
   fPadSizeZ = fPar->GetPadSizeZ();
   fPadRows = fPar->GetPadRows();
   fPadLayers = fPar->GetPadLayers();
   fNumTbs = fPar->GetNumTbs();
   fTBTime = fPar->GetTBTime();
   fDriftVelocity = fPar->GetDriftVelocity();
   fMaxDriftLength = fPar->GetDriftLength();
   fBField = fPar->GetBField();
   fEField = fPar->GetEField();
   fTiltAng = fPar->GetTiltAngle();
   fTB0 = fPar->GetTB0();
   fZk = fPar->GetZPadPlane();
   fEntTB = (Int_t)fPar->GetTBEntrance();

   // Was hardcoded to -103.0 deg, which contradicts both AtDigiPar (ThetaPad) and every
   // other user of the parameter (e.g. AtMCQMinimization, which rotates by +GetThetaPad()).
   fThetaPad = fPar->GetThetaPad() * TMath::Pi() / 180.0;
   fTiltAzim = fPar->GetThetaRot();

   std::cout << " ==== Parameters for Pulse Shape Analysis Task ==== " << std::endl;
   std::cout << " ==== Magnetic Field : " << fBField << " T " << std::endl;
   std::cout << " ==== Electric Field : " << fEField << " V/cm " << std::endl;
   std::cout << " ==== Sampling Rate : " << fTBTime << " ns " << std::endl;
   std::cout << " ==== Tilting Angle : " << fTiltAng << " deg " << std::endl;
   std::cout << " ==== Drift Velocity : " << fDriftVelocity << " cm/us " << std::endl;
   std::cout << " ==== TB0 : " << fTB0 << std::endl;
   std::cout << " ==== NumTbs : " << fNumTbs << std::endl;
}

void AtPSA::SetSimulatedEvent(TClonesArray *MCSimPointArray)
{
   fMCSimPointArray = MCSimPointArray;
}

void AtPSA::SetThreshold(Int_t threshold)
{
   fThreshold = threshold;
   if (!fUsingLowThreshold)
      fThresholdlow = threshold;
}

void AtPSA::SetThresholdLow(Int_t thresholdlow)
{
   fThresholdlow = thresholdlow;
   fUsingLowThreshold = kTRUE;
}

Double_t AtPSA::CalculateX(Double_t row)
{
   return (row + 0.5) * fPadSizeX - fPadPlaneX / 2.;
}

Double_t AtPSA::CalculateZ(Double_t peakIdx)
{
   // DEPRECAtED
   return (fNumTbs - peakIdx) * fTBTime * fDriftVelocity / 100.;
}

Double_t AtPSA::CalculateZGeo(Double_t peakIdx)
{

   // This function must be consistent with the re-calibrations done before.
   return fZk - (fEntTB - peakIdx) * fTBTime * fDriftVelocity / 100.;
}

Double_t AtPSA::CalculateY(Double_t layer)
{
   return (layer + 0.5) * fPadSizeZ;
}

// The drift time must be referenced to the same origin as CalculateZGeo, otherwise the
// transverse and longitudinal corrections disagree about where z=0 is. CalculateZGeo puts
// z=0 at tb = fEntTB - fZk*100/(fTBTime*fDriftVelocity); fTB0 (98, flagged DEPRECATED in
// the parameter file) is not that value. Note the drift-time origin only translates the
// event -- it cannot rotate it -- so it moves the vertex, never a track direction.
Double_t AtPSA::DriftTimeUs(Int_t tb) const
{
   Double_t tbZero = fEntTB - fZk * 100.0 / (fTBTime * fDriftVelocity);
   return (tb - tbZero) * fTBTime * 1E-3;
}

Double_t AtPSA::CalculateXCorr(Double_t xvalue, Int_t Tbx)
{
   return xvalue - fLorentzVector.X() * DriftTimeUs(Tbx) * 10.0; // cm/us * us -> cm -> mm
}

Double_t AtPSA::CalculateYCorr(Double_t yvalue, Int_t Tby)
{
   return yvalue - fLorentzVector.Y() * DriftTimeUs(Tby) * 10.0;
}

Double_t AtPSA::CalculateZCorr(Double_t zvalue, Int_t Tbz)
{
   Double_t zcorr = fLorentzVector.Z() * (Tbz - fTB0) * fTBTime * 1E-2;
   return zcorr;
}

void AtPSA::CalcLorentzVector()
{
   // Single source of truth, shared with AtClusterizeTask (the simulation's forward
   // model). See AtParameter/AtLangevin.h for the physics and for why the pad-frame
   // rotation and the general azimuth are both necessary.
   auto v = AtTools::LangevinDrift(fDriftVelocity, fBField, fEField, fTiltAng, fTiltAzim,
                                   fThetaPad * 180.0 / TMath::Pi());
   fLorentzVector.SetXYZ(v.x, v.y, v.z);

   static bool reported = false;
   if (!reported) {
      reported = true;
      std::cout << " ==== Lorentz drift vector (pad frame) : (" << v.x << ", " << v.y << ", " << v.z
                << ") cm/us,  omega*tau = " << v.omegaTau << std::endl;
   }
}

TVector3 AtPSA::RotateDetector(Double_t x, Double_t y, Double_t z, Int_t tb)
{

   // DEPRECAtED because of timebucket calibration (-271.0)
   TVector3 posRot;
   TVector3 posDet;

   posRot.SetX(x * TMath::Cos(fThetaPad) - y * TMath::Sin(fThetaPad));
   posRot.SetY(x * TMath::Sin(fThetaPad) + y * TMath::Cos(fThetaPad));
   posRot.SetZ((-271.0 + tb) * fTBTime * fDriftVelocity / 100. + fZk);

   Double_t TiltAng = -fTiltAng * TMath::Pi() / 180.0;

   posDet.SetX(posRot.X());
   posDet.SetY(-(fZk - posRot.Z()) * TMath::Sin(TiltAng) + posRot.Y() * TMath::Cos(TiltAng));
   posDet.SetZ(posRot.Z() * TMath::Cos(TiltAng) - posRot.Y() * TMath::Sin(TiltAng));

   return posDet;
}

void AtPSA::SetGainCalibration(TString gainFile)
{
   fCalibration->SetGainFile(gainFile);
}

void AtPSA::SetJitterCalibration(TString jitterFile)
{
   fCalibration->SetJitterFile(jitterFile);
}

void AtPSA::SetTBLimits(std::pair<Int_t, Int_t> limits)
{
   if (limits.first >= limits.second) {
      std::cout << " Warning AtPSA::SetTBLimits -  Wrong Time Bucket limits. Setting default limits (0,512) ... "
                << "\n";
      fIniTB = 0;
      fEndTB = 512;

   } else {
      fIniTB = limits.first;
      fEndTB = limits.second;
   }
}

void AtPSA::TrackMCPoints(std::multimap<Int_t, std::size_t> &map, AtHit *hit)
{
   typedef std::multimap<Int_t, std::size_t>::iterator MCMapIterator;

   // Find every simulated point ID for each valid pad
   std::pair<MCMapIterator, MCMapIterator> result = map.equal_range(hit->GetHitPadNum());
   for (MCMapIterator it = result.first; it != result.second; it++) {
      // std::cout<<fMCSimPointArray->GetEntries()<<"\n";

      int count = std::distance(result.first, result.second);

      // if(count>1){
      //  std::cout<<" Count "<<count<<"\n";

      if (fMCSimPointArray != 0) {
         AtTpcPoint *MCPoint = (AtTpcPoint *)fMCSimPointArray->At(it->second);
         AtHit::MCSimPoint mcpoint(it->second, MCPoint->GetTrackID(), MCPoint->GetEIni(), MCPoint->GetEnergyLoss(),
                                   MCPoint->GetAIni(), MCPoint->GetMassNum(), MCPoint->GetAtomicNum());
         hit->SetMCSimPoint(mcpoint);
         // std::cout << " Pad Num : "<<hit->GetHitPadNum()<<" MC Point ID : "<<it->second << std::endl;
         // std::cout << " Track ID : "<<MCPoint->GetTrackID()<<" Energy (MeV) : "<<MCPoint->GetEIni()<<" Angle (deg) :
         // "<<MCPoint->GetAIni()<<"\n"; std::cout << " Mass Number : "<<MCPoint->GetMassNum()<<" Atomic Number
         // "<<MCPoint->GetAtomicNum()<<"\n";
      }
      //  }
   }
}

ClassImp(AtPSA)
