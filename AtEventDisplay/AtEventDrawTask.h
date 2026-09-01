/**
 * @brief Event display task
 * Used to draw AT-TPC events
 * @author JungWoo Lee (Korea Univ.)
 *         Adapted for AtTPCROOT by Yassid Ayyad (NSCL)
 */
#ifndef ATEVENTDRAWTASK_H
#define ATEVENTDRAWTASK_H

// ATTPCROOT classes
// NOTE: the RANSAC / Hough / pattern-recognition / analysis classes are built only when
// PCL is available (see AtReconstruction/CMakeLists.txt). This display is deliberately
// PCL-free, so it does not reference them at all.
class AtEvent;
class AtEventManager;
class AtHit;
class AtRawEvent;
class AtMap;
class AtTpcMap;
class AtPatternEvent;

// FairRoot classes
#include "FairTask.h"

// ROOT classes
class TCanvas;
class TClonesArray;
class TEvePointSet;
class TEveGeoShape;
class TEveBoxSet;
class TEveLine;
class TGraph;
class TH1I;
class TH1D;
class TH1F;
class TH2D;
class TH2F;
class TH2Poly;
class TH3F;
class TPaletteAxis;
class TVector3;

#include "TEveRGBAPalette.h"

#include <Rtypes.h>
#include <fstream>
#include "TF1.h"

#ifndef __CINT__ // Boost
#include <boost/multi_array.hpp>
#endif //__CINT__

enum DetectorId { kAtTpc, kGADGETII, kSpecMAT };

class AtEventDrawTask : public FairTask {
public:
   AtEventDrawTask();
   AtEventDrawTask(TString modes);

   ~AtEventDrawTask();

   InitStatus Init();
   void Exec(Option_t *option);
   void Reset();

   // void Set2DPlotRange(Int_t uaIdx);
   void SetThreshold(Int_t val) { fThreshold = val; }
   void UnpackHoughSpace() { fUnpackHough = kTRUE; }
   void SetHitAttributes(Color_t, Size_t, Style_t);
   void Set3DHitStyleBar();
   void Set3DHitStyleBox();
   void SetSaveTextData();
   void SetLine(double t, std::vector<Double_t> p, double &x, double &y, double &z);
   void SetLine6(double t, std::vector<Double_t> p, double &x, double &y, double &z);
   void SetRawEventBranch(TString branchName);
   void SetEventBranch(TString branchName);

   static void SelectPad(const char *rawevt);
   /// Tilt for the beam-frame overlay, in degrees. Set this from the macro: the display
   /// cannot open a parameter file (FairParAsciiFileIo is an incomplete type in cling once
   /// the Eve dictionaries are loaded), so AtDigiPar is only used when it happens to be
   /// available. azimDeg is the beam azimuth in the PAD plane -- see AtLangevin.h.
   void SetTilt(Double_t tiltDeg, Double_t azimDeg)
   {
      fTiltAng = tiltDeg;
      fThetaRot = azimDeg;
      fHasTilt = (tiltDeg != 0.);
   }

   void DrawWave(Int_t PadNum);
   void SetMultiHit(Int_t hitMax);
   void SetAlgorithm(Int_t val) { fRANSACAlg = val; };
   inline void SelectDetectorId(DetectorId val) { fDetectorId = val; }

private:
   void DrawPadPlane();
   void DrawPadWave();
   void DrawPadAll();
   void DrawQEvent();
   void DrawRhoVariance();
   void DrawHoughSpace();
   void DrawPhiReco();
   void DrawMesh();
   void Draw3DHist();
   void DrawRad();
   void DrawTheta();
   void DrawThetaxPhi();
   void DrawMC();
   void DrawAux();

   AtMap *fAtMapPtr;
   void UpdateCvsPadPlane();
   void UpdateCvsPadWave();
   void UpdateCvsPadAll();
   void UpdateCvsQEvent();
   void UpdateCvsRhoVariance();
   void UpdateCvsHoughSpace();
   void UpdateCvsPhi();
   void UpdateCvsMesh();
   void UpdateCvs3DHist();
   void UpdateCvsRad();
   void UpdateCvsTheta();
   void UpdateCvsThetaxPhi();
   void UpdateCvsQuadrants();
   void UpdateCvsMC();
   void UpdateCvsAux();

   void ResetPadAll();
   void ResetPhiDistr();

   void DrawHitPoints();
   void DrawHSpace();
   void DrawMeshSpace();
   // void DrawHitClusterPoints();
   // void DrawRiemannHits();

   EColor GetTrackColor(int i);

   Bool_t fIs2DPlotRange;
   Bool_t fUnpackHough;
   Bool_t fIsCircularHough;
   Bool_t fIsLinearHough;
   static const Int_t fNumPads = 1000; // Maximum number of pads to draw for DrawAllPads option

   TString fRawEventBranchName;
   TString fEventBranchName;

   TClonesArray *fHitArray;
   TClonesArray *fRawEventArray;
   TClonesArray *fPatternEventArray;

   AtEventManager *fEventManager;
   AtRawEvent *fRawevent;

   AtMap *fDetmap;

   Int_t fThreshold;
   TString fMap;

   TEvePointSet *fHitSet;
   TEvePointSet *fHitSetMin;

   // Beam-frame overlay. The display used to TOGGLE between raw and corrected hits, so the
   // two could never be compared in one picture. These are drawn together instead:
   //   fHitSetCorr  Lorentz/tilt-corrected hits (GetPositionCorr)
   //   fHitSetRot   the same hits rotated onto the detector axis -- the "as if the detector
   //                had never been tilted" view that validates the correction on data
   //   fBeamAxis    the beam direction (== B), tilted by TiltAng from the pad-plane normal
   //   fDetAxis     the detector axis itself, for reference
   TEvePointSet *fHitSetCorr = nullptr;
   TEvePointSet *fHitSetRot = nullptr;
   TEveLine *fBeamAxis = nullptr;
   TEveLine *fDetAxis = nullptr;

   Double_t fTiltAng = 0.;  ///< deg, from AtDigiPar
   Double_t fThetaRot = 0.; ///< deg, azimuth of the beam in the PAD plane (see AtLangevin.h)
   Bool_t fHasTilt = kFALSE;
   TVector3 RotateToAxis(const TVector3 &v) const;

   TEvePointSet *fHitSetMC[5];    // For MC results
   TEvePointSet *fHitSetTFHC[20]; // for TrackFinderHC
   TEveBoxSet *fHitClusterSet[20]; // Track clusterization

   // TEveGeoShape* x;
   // std::vector<TEveGeoShape*> hitSphereArray;

   TEveBoxSet *fhitBoxSet;

   TPaletteAxis *fPadPlanePal;

   Color_t fHitColor;
   Size_t fHitSize;
   Style_t fHitStyle;

   TCanvas *fCvsPadPlane;
   TH2Poly *fPadPlane;
   TCanvas *fCvsPadWave;
   TH1I *fPadWave;
   TCanvas *fCvsPadAll;
   TH1I *fPadAll[fNumPads];
   TCanvas *fCvsQEvent;
   TH1D *fQEventHist;
   TH1D *fQEventHist_H;
   TCanvas *fCvsHoughSpace;
   TH2F *fHoughSpace;
   TCanvas *fCvsRhoVariance;
   TH1D *fRhoVariance;
   TCanvas *fCvsPhi;
   TH1D *fPhiDistr[5];
   TCanvas *fCvsMesh;
   TH1F *fMesh;
   TCanvas *fCvs3DHist;
   TH3F *f3DHist;
   TCanvas *fCvsRad;
   TH2F *fRadVSTb;
   TCanvas *fCvsTheta;
   TH2F *fTheta;
   TCanvas *fCvsThetaxPhi;
   TH2F *fThetaxPhi;
   TCanvas *fCvsQuadrant1;
   TH2F *fQuadrant1;
   TCanvas *fCvsQuadrant2;
   TH2F *fQuadrant2;
   TCanvas *fCvsQuadrant3;
   TH2F *fQuadrant3;
   TCanvas *fCvsQuadrant4;
   TH2F *fQuadrant4;
   TH1F *fAuxChannels[9];
   TCanvas *fCvsAux;

   TH2F *fThetaxPhi_Ini;
   TH2F *fThetaxPhi_Ini_RANSAC;

   TCanvas *fCvsMC_XY;
   TGraph *fMC_XY;
   TGraph *fMC_XY_exp;
   TGraph *fMC_XY_int;
   TGraph *fMC_XY_back;
   TCanvas *fCvsMC_Z;
   TGraph *fMC_ZX;
   TGraph *fMC_ZX_int;
   TGraph *fMC_ZX_back;
   TGraph *fMC_ZY;
   TGraph *fMC_ZY_int;
   TGraph *fMC_ZY_back;

   Int_t fNQuads;

   Int_t fMinZ;
   Int_t fMaxZ;
   Int_t fMinX;
   Int_t fMaxX;

   Int_t f3DHitStyle;
   Int_t fMultiHit;
   Bool_t fSaveTextData;
   Float_t f3DThreshold;
   Bool_t fIsRawData;
   Int_t fRANSACAlg;
   DetectorId fDetectorId;
   Int_t fDetNumPads;
   TF1 *fHoughLinearFit;
   TF1 *fRansacLinearFit;
   AtHit const *fIniHit;
   AtHit const *fIniHitRansac;

   // std::vector<TEveLine*> fLineArray;
   TEveLine *fLineArray[5];
   TEvePointSet *fVertex = nullptr;
   Int_t fLineNum;
   Int_t fTrackNum;
   // TEveLine* fLine;

   TEveRGBAPalette *fRGBAPalette;

   ClassDef(AtEventDrawTask, 1);
};

#endif
