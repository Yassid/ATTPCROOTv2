/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *         GNU Lesser General Public Licence version 3 (LGPL) version 3,        *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// in root all sizes are given in cm

#include "TFile.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"
#include "TGeoVolume.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"

#include <iostream>

// Name of geometry version and output file
const TString geoVersion = "nOTPC_CF4_250mbar";
const TString FileName = geoVersion + ".root";
const TString FileName1 = geoVersion + "_geomanager.root";

// Names of the different used materials which are used to build the modules
// The materials are defined in the global media.geo file
const TString MediumGas = "CF4_250mbar"; //"Ar90CF4_250mbar";
const TString MediumVacuum = "vacuum4";
// const TString MediumWindow = "aramid";
const TString CylinderVolumeMedium = "steel";
const TString FlangeMaterial = "Aluminum5083";
const TString FieldCageMaterial = "Epoxy";
const TString ScintillatorWindowMaterial = "Quartz";
const TString ScintillatorCrystalMaterial = "CeBr3";
const TString ScintillatorReflectorMaterial = "TiO2";
const TString ModeratorMaterial = "polypropylene";

// Parameters for SpecMAT detector geometry
/*const int TotalCrystNb = 45; // Number of scintillators in the scintillator array
TGeoVolumeAssembly *sci_seg[15];
const TString geoSegments = "sci_seg_";*/

// some global variables
TGeoManager *gGeoMan = new TGeoManager("ATTPC", "ATTPC");
;                     // Pointer to TGeoManager instance
TGeoVolume *gModules; // Global storage for module types

// Forward declarations
void create_materials_from_media_file();
TGeoVolume *create_detector();
void position_detector();
void add_alignable_volumes();

void nOTPC_CF4_250mbar()
{
   // Load the necessary FairRoot libraries
   // gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
   // basiclibs();
   gSystem->Load("libGeoBase");
   gSystem->Load("libParBase");
   gSystem->Load("libBase");

   // Load needed material definition from media.geo file
   create_materials_from_media_file();

   // Get the GeoManager for later usage
   gGeoMan = (TGeoManager *)gROOT->FindObject("FAIRGeom");
   gGeoMan->SetVisLevel(7);

   // Create the top volume

   TGeoVolume *top = new TGeoVolumeAssembly("TOP");
   gGeoMan->SetTopVolume(top);

   TGeoMedium *gas = gGeoMan->GetMedium(MediumVacuum);
   TGeoVolume *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(gas);
   top->AddNode(tpcvac, 1);

   // Build the detector
   gModules = create_detector();

   // position_detector();

   cout << "Voxelizing." << endl;
   top->Voxelize("");
   gGeoMan->CloseGeometry();

   // add_alignable_volumes();

   gGeoMan->CheckOverlaps(0.001);
   gGeoMan->PrintOverlaps();
   gGeoMan->Test();

   TFile *outfile = new TFile(FileName, "RECREATE");
   top->Write();
   outfile->Close();

   TFile *outfile1 = new TFile(FileName1, "RECREATE");
   gGeoMan->Write();
   outfile1->Close();

   top->Draw("ogl");
   // top->Raytrace();
}

void create_materials_from_media_file()
{
   // Use the FairRoot geometry interface to load the media which are already defined
   FairGeoLoader *geoLoad = new FairGeoLoader("TGeo", "FairGeoLoader");
   FairGeoInterface *geoFace = geoLoad->getGeoInterface();
   TString geoPath = gSystem->Getenv("VMCWORKDIR");
   TString geoFile = geoPath + "/geometry/media.geo";
   geoFace->setMediaFile(geoFile);
   geoFace->readMedia();

   const TString MediumGas = "CF4_250mbar";
   const TString CylinderVolumeMedium = "steel";
   const TString MediumVacuum = "vacuum4";
   const TString FieldCageCyl = "Aluminum";
   const TString ModeratorMedium = "polypropylene";
   //   const TString fc_rings = "Aluminum";

   // Read the required media and create them in the GeoManager
   FairGeoMedia *geoMedia = geoFace->getMedia();
   FairGeoBuilder *geoBuild = geoLoad->getGeoBuilder();

   FairGeoMedium *ATTPCGas = geoMedia->getMedium(MediumGas);
   FairGeoMedium *steel = geoMedia->getMedium("steel");
   FairGeoMedium *vacuum4 = geoMedia->getMedium("vacuum4");
   // FairGeoMedium *aramid = geoMedia->getMedium("aramid");
   FairGeoMedium *Aluminum5083 = geoMedia->getMedium(FlangeMaterial);
   FairGeoMedium *Epoxy = geoMedia->getMedium(FieldCageMaterial);
   FairGeoMedium *TiO2 = geoMedia->getMedium(ScintillatorReflectorMaterial);
   FairGeoMedium *CeBr3 = geoMedia->getMedium(ScintillatorCrystalMaterial);
   FairGeoMedium *Quartz = geoMedia->getMedium(ScintillatorWindowMaterial);
   FairGeoMedium *heco2 = geoMedia->getMedium("heco2");
   FairGeoMedium *G10 = geoMedia->getMedium("G10");
   FairGeoMedium *Aluminum = geoMedia->getMedium("Aluminum");
   FairGeoMedium *PP = geoMedia->getMedium("polypropylene");

   /***** Unused IC media *****
   FairGeoMedium* ICpolypropylene   = geoMedia->getMedium("ICpolypropylene");
   FairGeoMedium* ICIsoButane       = geoMedia->getMedium("ICIsoButane");
   FairGeoMedium* Aluminium         = geoMedia->getMedium("Aluminum");
   geoBuild->createMedium(ICpolypropylene);
   geoBuild->createMedium(ICIsoButane);
   // geoBuild->createMedium(aramid);
   */

   // include check if all media are found
   geoBuild->createMedium(ATTPCGas);
   geoBuild->createMedium(steel);
   geoBuild->createMedium(vacuum4);
   geoBuild->createMedium(Aluminum5083);
   geoBuild->createMedium(Aluminum);
   geoBuild->createMedium(Epoxy);
   geoBuild->createMedium(TiO2);
   geoBuild->createMedium(CeBr3);
   geoBuild->createMedium(Quartz);
   geoBuild->createMedium(G10);
   geoBuild->createMedium(PP);
}

TGeoVolume *create_detector()
{

   // needed materials
   TGeoMedium *OuterCylinder = gGeoMan->GetMedium(CylinderVolumeMedium);
   TGeoMedium *gas = gGeoMan->GetMedium(MediumGas);
   // TGeoMedium *windowmat = gGeoMan->GetMedium(MediumWindow);
   TGeoMedium *Vacuummat = gGeoMan->GetMedium("vacuum4");
   TGeoMedium *Aluminum5083mat = gGeoMan->GetMedium(FlangeMaterial);
   TGeoMedium *Epoxymat = gGeoMan->GetMedium(FieldCageMaterial);
   TGeoMedium *TiO2mat = gGeoMan->GetMedium(ScintillatorReflectorMaterial);
   TGeoMedium *CeBr3mat = gGeoMan->GetMedium(ScintillatorCrystalMaterial);
   TGeoMedium *Quartzmat = gGeoMan->GetMedium(ScintillatorWindowMaterial);
   TGeoMedium *ModMat = gGeoMan->GetMedium(ModeratorMaterial);
   // TGeoMedium *fc = gGeoMan->GetMedium(FieldCageCyl);

   /*** Unused IC media ***
   TGeoMedium* ICgas         = gGeoMan->GetMedium(ICMediumGas);
   TGeoMedium* ICwindowmat   = gGeoMan->GetMedium(ICWindowMedium);
   TGeoMedium* ICAlwindowmat = gGeoMan->GetMedium(ICAlWindowMedium);
   */

   //  Main drift volume

   // double tpc_rot = TMath::Pi();
   double tpc_rot = 0;
   double active_volume_radius = 25.0;
   double active_volume_halflength = 25.0 / 2.0;

   TGeoVolume *drift_volume =
      gGeoManager->MakeTube("drift_volume", gas, 0, active_volume_radius, active_volume_halflength);
   // TGeoVolume *drift_volume = gGeoManager->MakeBox("drift_volume", gas,  100./2, 100./2, 100./2);
   // gGeoMan->GetVolume(geoVersion)->AddNode(drift_volume,1, new TGeoTranslation(0,0,drift_length/2));
   drift_volume->SetLineColor(kCyan);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(drift_volume, 1, new TGeoCombiTrans(0, 0, 0.0, new TGeoRotation("drift_volume", 0, tpc_rot, 0)));
   drift_volume->SetTransparency(80);

   double fieldcage_innerradius = active_volume_radius;
   double fieldcage_outerradius = 25.65;
   double fieldcage_halflength = 25.2 / 2.0;

   TGeoVolume *field_cage =
      gGeoManager->MakeTube("field_cage", Epoxymat, active_volume_radius, fieldcage_outerradius, fieldcage_halflength);
   field_cage->SetLineColor(kOrange + 2);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(field_cage, 2, new TGeoCombiTrans(0.0, 0.0, 0.0, new TGeoRotation("field_cage", 0, tpc_rot, 0)));
   field_cage->SetTransparency(35);

   TGeoVolume *moderator = gGeoManager->MakeTube("drift_volume_moderator", ModMat, 0, active_volume_radius, 50.0);

   moderator->SetLineColor(kGray);
   moderator->SetTransparency(80);
   gGeoMan->GetVolume(geoVersion)
      ->AddNode(moderator, 3, new TGeoCombiTrans(0, 0, 100.0, new TGeoRotation("moderator", 0, tpc_rot, 0)));
   drift_volume->SetTransparency(80);

   return drift_volume;
}

void position_detector()
{

   /* TGeoTranslation* det_trans=NULL;

    Int_t numDets=0;
    for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {
      det_trans
        = new TGeoTranslation("", 0., 0., First_Z_Position+(numDets*Z_Distance));
      gGeoMan->GetVolume(geoVersion)->AddNode(gModules, numDets, det_trans);
      numDets++;

    }*/
}

void add_alignable_volumes()
{

   /* TString volPath;
    TString symName;
    TString detStr   = "Tutorial4/det";
    TString volStr   = "/TOP_1/tutorial4_1/tut4_det_";

    for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {

      volPath  = volStr;
      volPath += detectorPlanes;

      symName  = detStr;
      symName += Form("%02d",detectorPlanes);

      cout<<"Path: "<<volPath<<", "<<symName<<endl;
  //    gGeoMan->cd(volPath);

      gGeoMan->SetAlignableEntry(symName.Data(),volPath.Data());

    }
      cout<<"Nr of alignable objects: "<<gGeoMan->GetNAlignable()<<endl;*/
}
