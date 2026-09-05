/// @file Ar46_telescope.C
/// @brief Forward dE-E telescope for 46Ar(3He,d)47K: two position-sensitive DSSDs behind the
///        cathode of the REVERSED AT-TPC.
///
///   root -l geometry/Ar46_telescope.C
///
/// In the reversed detector the beam enters through the pad plane and leaves through the cathode,
/// so the exit aperture can be large without removing pads -- which is what makes this telescope
/// possible. It sits just past the end of the drift volume and catches the heavy residual.
///
/// THICKNESSES: 20 um dE + 200 um E, SET BY YASSID (2026-09-05). Si_forward_telescope.C's
/// 500 um dE is far too thick for this beam whichever stopping model you believe.
///
/// !! THE STOPPING POWER HERE IS UNRESOLVED -- DO NOT TREAT THE NUMBERS BELOW AS SETTLED. !!
/// Two calculations disagree by a factor ~4.5 for a 461 MeV 47K in silicon:
///
///     AtELossCATIMA        60 MeV in 100 um   (range ~570 um)
///     GEANT4, as measured in this simulation      279 MeV in 100 um   (stops in ~150-200 um)
///
/// A scaling of this simulation's own gas energy loss (the beam loses 96 MeV over 8.3e-3 g/cm2,
/// i.e. ~11600 MeV/(g/cm2), scaled to Si by Z/A) lands near the GEANT4 value -- but Yassid has
/// flagged the argument as not right and it is to be revisited. So the thicknesses here are a
/// DECISION, not a derivation, and the CATIMA/GEANT4 discrepancy is an open item that also
/// affects AtELossCATIMA wherever else this analysis uses it.
///
/// The CsI of the original is dropped -- nothing gets that far under either model.
///
/// WHY A dE-E TELESCOPE AT ALL: Z TAGGING. The residual (Z = 19) and the unreacted beam (Z = 18)
/// arrive at the same place with similar energies, and dE goes as Z^2, so they separate by
/// (19/18)^2 - 1 = 11.4 % in dE. That ratio is INDEPENDENT of the dE thickness, so the thickness
/// is chosen purely on signal size and punch-through, as above.
///
/// SIZE AND STANDOFF. The residual stays within 3.33 deg of the axis over the proposal's
/// theta_cm 15-80 deg window, but the vertex is spread over the whole metre of gas, so the lever
/// arm runs from ~100 cm (reaction at the pad plane) down to a few cm (reaction at the cathode)
/// and the spot radius with it. For a 10 x 10 cm DSSD:
///
///     gap beyond the cathode:   3 cm    5 cm   10 cm   15 cm   20 cm
///     acceptance:              96.3 %  95.6 % 93.5 %  91.3 %  88.9 %
///
/// and what is clipped is always the large-theta_cm corner, where the DWBA is weakest -- so the
/// real loss is smaller than those flat numbers. 5 cm is the default here as a compromise with
/// whatever the vessel actually needs; move it with kZFront.
///
/// POSITION. Both layers are position sensitive. Strip pitch is NOT in this geometry -- Geant4
/// records the true hit position in AtSiPoint and the strips are applied in analysis, which is
/// also how the AtSiArray/HELIOS setup works. Do not go finer than ~1 mm when you do: the
/// residual crosses up to a metre of gas first, giving 0.084 deg of multiple scattering, and with
/// the recoil Jacobian of 0.04 deg(lab)/deg(cm) that is already ~2 deg in theta_cm. A 1 mm strip
/// at a 70 cm lever is also ~2 deg. Finer strips would measure the scattering, not the reaction.
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"

#include <iostream>

const TString geoVersion = "Ar46_telescope_v1.0";
const TString FileName = geoVersion + ".root";
const TString FileName1 = geoVersion + "_geomanager.root";

const TString MediumSi = "silicon";
const TString MediumVacuum = "vacuum4";

TGeoManager *gGeoMan = new TGeoManager("Ar46Tel", "Ar46Tel");

// --- the numbers, all in cm, all justified in the header --------------------------------------
const Double_t kXSize = 10.0;  ///< active width  [cm]
const Double_t kYSize = 10.0;  ///< active height [cm]
const Double_t kDriftEnd = 100.0; ///< z of the cathode = end of the drift volume [cm]
const Double_t kGap = 5.0;     ///< clearance between the cathode and the first DSSD [cm]
const Double_t kZFront = kDriftEnd + kGap;
const Double_t kDEThick = 0.0020; ///< dE DSSD,  20 um  (set by Yassid; see the header)
const Double_t kEThick = 0.0200;  ///< E  DSSD, 200 um  (set by Yassid; see the header)
const Double_t kSep = 1.0;        ///< gap between the two DSSDs [cm]

void create_materials_from_media_file();
TGeoVolume *create_detector();

void Ar46_telescope()
{
   create_materials_from_media_file();

   gGeoMan = (TGeoManager *)gROOT->FindObject("FAIRGeom");
   gGeoMan->SetVisLevel(7);

   TGeoVolume *top = new TGeoVolumeAssembly("TOP");
   gGeoMan->SetTopVolume(top);

   TGeoMedium *vac = gGeoMan->GetMedium(MediumVacuum);
   TGeoVolume *topvac = new TGeoVolumeAssembly(geoVersion);
   topvac->SetMedium(vac);
   top->AddNode(topvac, 1);

   create_detector();

   std::cout << "Voxelizing." << std::endl;
   top->Voxelize("");
   gGeoMan->CloseGeometry();

   // Overlaps are checked at a 10 um tolerance because the dE layer is only 100 um thick: the
   // default 1 mm would be larger than the object being placed.
   gGeoMan->CheckOverlaps(0.001);
   gGeoMan->PrintOverlaps();
   gGeoMan->Test();

   TFile *outfile = new TFile(FileName, "RECREATE");
   top->Write();
   outfile->Close();

   TFile *outfile1 = new TFile(FileName1, "RECREATE");
   gGeoMan->Write();
   outfile1->Close();

   std::cout << "\n  wrote " << FileName << " and " << FileName1 << "\n"
             << "  dE  DSSD: " << kDEThick * 1e4 << " um at z = " << kZFront << " cm\n"
             << "  E   DSSD: " << kEThick * 1e4 << " um at z = " << kZFront + kSep << " cm\n"
             << "  active  : " << kXSize << " x " << kYSize << " cm, "
             << kGap << " cm beyond the cathode at z = " << kDriftEnd << " cm\n\n";

   // Opens the OpenGL viewer when run interactively (root -l). Skipped under -b so the batch
   // rebuild used by the campaign scripts stays headless.
   if (!gROOT->IsBatch())
      top->Draw("ogl");
}

void create_materials_from_media_file()
{
   FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   TString geoPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = geoPath + "/geometry/media.geo";
   geoFace->setMediaFile(geoFile);
   geoFace->readMedia();

   FairGeoMedia *geoMedia = geoFace->getMedia();
   FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

   // silicon and vacuum4 already exist in media.geo -- NOTHING is appended to that file. Adding
   // media past Ar90CF4_250mbar has hung FairGeoMedia and taken the machine out of memory before.
   FairGeoMedium *silicon = geoMedia->getMedium("silicon");
   FairGeoMedium *vacuum4 = geoMedia->getMedium("vacuum4");
   if (!silicon || !vacuum4) {
      std::cerr << "MISSING MEDIUM: silicon or vacuum4 not found in media.geo\n";
      return;
   }
   geoBuild->createMedium(silicon);
   geoBuild->createMedium(vacuum4);
}

TGeoVolume *create_detector()
{
   TGeoMedium *silicon = gGeoMan->GetMedium(MediumSi);

   // !! THE VOLUME NAME MUST CONTAIN THE SUBSTRING "silicon" OR THE DETECTOR RECORDS NOTHING. !!
   // AtSiArray::CheckIfSensitive (AtSiArray.cxx) is literally
   //     if (tsname.Contains("silicon")) return kTRUE;
   // and FairDetector only calls ProcessHits for volumes that test declares sensitive. Name them
   // anything else and the telescope is still built, still drawn, still reported as "Constructing
   // Si Array geometry from ROOT file", and produces ZERO hits with no error whatsoever -- which
   // is exactly what the first version of this file did with "Ar46_dE"/"Ar46_E". Note that
   // geometry/Si_forward_telescope.C names its layers "dESi"/"ESi" and so has the same latent
   // defect; it appears never to have been used.
   //
   // Beyond that, the name is the contract with the analysis: AtSiArray::ProcessHits stores the
   // volume name on every AtSiPoint, and that is how dE and E are told apart downstream.
   TGeoVolume *dE = gGeoManager->MakeBox("silicon_Ar46_dE", silicon, kXSize / 2, kYSize / 2, kDEThick / 2);
   gGeoMan->GetVolume(geoVersion)->AddNode(dE, 0, new TGeoTranslation(0.0, 0.0, kZFront + kDEThick / 2));
   dE->SetLineColor(kGreen);

   TGeoVolume *E = gGeoManager->MakeBox("silicon_Ar46_E", silicon, kXSize / 2, kYSize / 2, kEThick / 2);
   gGeoMan->GetVolume(geoVersion)->AddNode(E, 0, new TGeoTranslation(0.0, 0.0, kZFront + kSep + kEThick / 2));
   E->SetLineColor(kAzure + 2);

   return dE;
}
