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
/// THE STACK IS THE REAL HARDWARE (Yassid, 2026-09-05):
///
///     dE    500 um silicon DSSD
///     E    1000 um silicon DSSD
///     CsI  array of 18 x 18 mm elements, 25 mm deep
///
/// which is also what geometry/Si_forward_telescope.C had for the silicon -- my earlier 20/200
/// suggestion came from my own stopping-power estimate and was wrong. The CsI is there precisely
/// because ions punch through 1.5 mm of silicon, so the CsI is not optional decoration; it is
/// where the remaining energy is measured.
///
/// !! AN OPEN DISCREPANCY, RECORDED SO IT IS NOT REDISCOVERED. !! GEANT4 in this simulation has
/// a 461 MeV 47K depositing 279 MeV in only 100 um of silicon, i.e. stopping in ~150-200 um,
/// which cannot be reconciled with a stack that deliberately puts CsI behind 1.5 mm of silicon.
/// AtELossCATIMA disagrees with GEANT4 by a factor ~4.5 on the same case. Something in that chain
/// is wrong and it has NOT been resolved -- so read the deposits this geometry produces as a
/// measurement to be checked against the real detector, not as a prediction to trust. The same
/// question hangs over AtELossCATIMA wherever else this analysis uses it.
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

TString geoVersion = "Ar46_telescope_v1.0";
TString FileName = geoVersion + ".root";
TString FileName1 = geoVersion + "_geomanager.root";

const TString MediumSi = "silicon";
const TString MediumVacuum = "vacuum4";
const TString MediumCsI = "CsI";

TGeoManager *gGeoMan = new TGeoManager("Ar46Tel", "Ar46Tel");

// --- the numbers, all in cm, all justified in the header --------------------------------------
const Double_t kXSize = 10.0;  ///< active width  [cm]
const Double_t kYSize = 10.0;  ///< active height [cm]
const Double_t kDriftEnd = 100.0; ///< z of the cathode = end of the drift volume [cm]
const Double_t kGap = 5.0;     ///< clearance between the cathode and the first DSSD [cm]
const Double_t kZFront = kDriftEnd + kGap;
Double_t kDEThick = 0.0500; ///< dE DSSD,  500 um  (overridden by the argument)
Double_t kEThick = 0.1000;  ///< E  DSSD, 1000 um  (overridden by the argument)
// CsI array: 18 x 18 mm entrance face, 25 mm deep. 5 x 5 covers 90 x 90 mm, i.e. the 10 x 10 cm
// silicon in front of it -- the element size is the hardware, the 5 x 5 count is an assumption
// and is the first thing to change if the real array differs.
const Int_t kNCsI = 5;            ///< elements per side
const Double_t kCsIFace = 1.8;    ///< entrance face [cm]
const Double_t kCsIDepth = 2.5;   ///< depth [cm]
const Double_t kCsIGap = 0.5;     ///< gap between the E DSSD and the CsI front face [cm]
const Double_t kSep = 1.0;        ///< gap between the two DSSDs [cm]

void create_materials_from_media_file();
TGeoVolume *create_detector();

/// @param dEum   dE silicon thickness in um. Default 500 = the hardware.
/// @param Eum    E  silicon thickness in um. Default 1000 = the hardware.
/// @param tag    appended to the geometry name, so alternative stacks get their OWN files instead
///               of overwriting the hardware one. Empty keeps Ar46_telescope_v1.0.
void Ar46_telescope(Double_t dEum = 500., Double_t Eum = 1000., TString tag = "")
{
   kDEThick = dEum * 1e-4;
   kEThick = Eum * 1e-4;
   if (tag.Length()) {
      geoVersion = "Ar46_telescope_" + tag;
      FileName = geoVersion + ".root";
      FileName1 = geoVersion + "_geomanager.root";
   }
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
             << "  CsI     : " << kNCsI << " x " << kNCsI << " of " << kCsIFace * 10 << " x "
             << kCsIFace * 10 << " mm, " << kCsIDepth * 10 << " mm deep, front face at z = "
             << kZFront + kSep + kEThick + kCsIGap << " cm\n"
             << "  Si active: " << kXSize << " x " << kYSize << " cm, "
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
   FairGeoMedium *csi = geoMedia->getMedium("CsI");
   if (!silicon || !vacuum4 || !csi) {
      std::cerr << "MISSING MEDIUM: silicon, vacuum4 or CsI not found in media.geo\n";
      return;
   }
   geoBuild->createMedium(silicon);
   geoBuild->createMedium(vacuum4);
   geoBuild->createMedium(csi);
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

   // --- CsI array ------------------------------------------------------------------------------
   // Each element is its own volume and its own node, so AtSiPoint carries which element fired and
   // the array is position sensitive at the element level without any strip decoding.
   //
   // The name must contain "CsI": AtSiArray::CheckIfSensitive tests substrings of the volume name,
   // and anything it does not recognise is silently non-sensitive.
   TGeoMedium *csi = gGeoMan->GetMedium(MediumCsI);
   if (!csi) {
      std::cerr << "MISSING MEDIUM: CsI not found -- the array will not be built\n";
      return dE;
   }
   const Double_t zCsI = kZFront + kSep + kEThick + kCsIGap;
   for (Int_t ix = 0; ix < kNCsI; ++ix) {
      for (Int_t iy = 0; iy < kNCsI; ++iy) {
         const Double_t x = (ix - (kNCsI - 1) / 2.0) * kCsIFace;
         const Double_t y = (iy - (kNCsI - 1) / 2.0) * kCsIFace;
         TGeoVolume *el = gGeoManager->MakeBox(Form("CsI_Ar46_%d_%d", ix, iy), csi, kCsIFace / 2,
                                               kCsIFace / 2, kCsIDepth / 2);
         el->SetLineColor(kYellow - 9);
         gGeoMan->GetVolume(geoVersion)->AddNode(el, ix * kNCsI + iy,
                                                 new TGeoTranslation(x, y, zCsI + kCsIDepth / 2));
      }
   }

   return dE;
}
