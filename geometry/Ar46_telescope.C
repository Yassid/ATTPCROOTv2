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
/// !! THE 500/1000 um STACK CANNOT WORK FOR THIS RESIDUAL, AND THE OLD "FACTOR 4.5" NOTE WAS
/// WRONG -- BUT THE CODES DO NOT AGREE EITHER. !! (measured 2026-09-09)
///
/// RANGE. The 47K residual arrives at 9.5-10.9 MeV/u (arrival KE 477 +- 14 MeV, measured) and its
/// range in silicon is 119-172 um. catima called directly gives 133.5 um at 461 MeV, and
/// AtELossCATIMA::GetRange gives 133.7 um -- that path is a direct catima::range call and is not
/// step-dependent. A Bethe-Barkas hand check agrees: a 10 MeV proton in Si loses 35.8 MeV cm2/g
/// (PSTAR says 35.5), and Z_eff = 17.5 for K at this velocity gives Z_eff^2 = 307, hence
/// 2561 MeV/mm -- catima's own entrance value to 0.1 %. The old note here claiming a 573 um range
/// is simply wrong.
///
/// dE/dx. GEANT4 in this simulation is LOW against all of that, by 15-30 %:
///
///     layer      arrival      GEANT4 dE      catima dE     catima/GEANT4
///     100 um     461 MeV       278.6 MeV      315.2 MeV        1.13
///      40 um     476 MeV        81.5 MeV      107.3 MeV        1.32
///
/// (the ratio is not constant because at 100 um the ion is already into its Bragg peak, and the
/// two codes put that peak in different places.) catima is the one that matches the hand check,
/// so GEANT4's heavy-ion stopping here is the suspect -- most likely the physics list's ion
/// parameterisation. NOT RESOLVED, but it is a 25 % question, not a factor of 4.5.
///
/// !! THE TRAP THAT MANUFACTURED THE OLD DISAGREEMENT !! AtELossCATIMA integrates in steps of
/// fRangeStepSize, DEFAULT 0.1 mm -- one single step for a 100 um layer. Its deposit for the case
/// above runs 256.0 / 279.1 / 306.8 / 315.2 MeV at steps of 0.1 / 0.05 / 0.01 / 0.001 mm. The
/// default underestimates by 19 %, and 256 MeV happens to sit close to GEANT4's 279, which is what
/// made the two look like they agreed. CALL SetRangeStepSize() BEFORE USING GetEnergy OR
/// GetEnergyLoss ON ANYTHING THINNER THAN A FEW MILLIMETRES. GetRange(E, 0) is unaffected.
///
/// WHAT THIS MEANS FOR THE STACK -- and it is the same answer under either code. The hardware
/// stack is a TOTAL-ENERGY detector for 47K and carries no dE-E information: measured on 12000
/// events per state, 5680 of 5688 residuals stop in the 500 um dE, NOTHING reaches the E layer,
/// and the CsI -- behind 1.5 mm of silicon -- never fires. The CsI is presumably there for light
/// ejectiles, but in this reaction the deuteron goes BACKWARD (theta_lab 60-134 deg) and never
/// reaches the telescope at all.
///
/// A WORKING STACK NEEDS A THIN dE. Build one with the thickness arguments:
///   root -l 'geometry/Ar46_telescope.C(40,200,"40_200")'
/// Measured on 12000 events: 5655 of 5665 residuals deposit in the 40 um dE and 5650 punch through
/// into the E layer -- 99.8 % -- with a mean dE of 81.4 MeV, and the CsI still empty. 200 um of E
/// stops the rest with margin. The species separate by ~9 % in dE per unit of Z (46Ar / 47K /
/// 48Ca), nearly independent of the dE thickness, while the energy-loss straggling is 0.05-0.3 %,
/// so the separation is limited by DETECTOR resolution, not physics: ~1 % silicon leaves it
/// ~9 sigma. That 9 % is a Z_eff^2 scaling both codes share, so the GEANT4/catima gap above does
/// not put it in doubt.
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
/// MEASURED, 2026-09-09, and it BEATS that table: 99.9 % of the 47K residuals deposit in the dE
/// (5680 of 5688 on 12000 events, at 2.85 T with a 5 cm gap). The analytic numbers above are
/// straight lines and are therefore pessimistic -- the solenoid FOCUSES a residual emitted from a
/// point on its own axis, and the median radius at the dE comes out 2.15 cm against a
/// straight-line 2.43 cm. Do not read the field as a nuisance here; it helps.
///
/// That 99.9 % only became visible once the residual stopped being truncated in flight: see
/// AtGenerators/AtTPCIonGenerator.cxx. Before that fix this telescope measured 46 %, and the
/// deficit was NOT geometric.
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

#include <cmath>
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
// CsI array: 18 x 18 mm entrance face, 25 mm deep -- the element size and depth are the hardware.
// THE COUNT IS DERIVED, NOT CHOSEN, so the array always COVERS THE DSSD: a fixed 5 x 5 spans only
// 90 x 90 mm against a 100 x 100 mm silicon and leaves the corners -- and therefore the largest
// theta_cm, which is where the residual lands -- unwatched. ceil() rounds up, so the array
// slightly overhangs the silicon rather than falling short, and it follows kXSize/kYSize
// automatically if the DSSD size is ever changed.
const Int_t kNCsI = (Int_t)std::ceil(std::max(kXSize, kYSize) / 1.8); ///< elements per side
const Double_t kCsIFace = 1.8;    ///< entrance face [cm] -- hardware
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

   // WRITE INTO $VMCWORKDIR/geometry, NOT THE CURRENT DIRECTORY. FairModule::SetGeometryFileName
   // searches only the standard geometry path, and a file it cannot find there is reported as
   // "[FATAL] ... not found in standard path" -- which does NOT stop the job. A 12000-event run
   // was produced that way with a telescope present in the module list, no geometry behind it and
   // ZERO AtSiPoints, and it looked like a physics result (0 of 5673 residuals in the dE) rather
   // than a missing file.
   TString outDir = gSystem->Getenv("VMCWORKDIR");
   if (outDir.IsNull())
      outDir = "..";
   outDir += "/geometry/";
   FileName = outDir + FileName;
   FileName1 = outDir + FileName1;

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
