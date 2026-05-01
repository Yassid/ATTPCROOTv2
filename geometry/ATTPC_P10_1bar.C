/********************************************************************************
 *    ATTPC geometry — P10 gas (90% Ar + 10% CH4 by volume) at 1 bar.
 *    Self-contained build: defines materials directly via TGeoMixture instead of
 *    going through FairGeoLoader (which hangs reading media.geo on this setup).
 *    Produces ATTPC_P10_1bar.root and ATTPC_P10_1bar_geomanager.root in $PWD.
 ********************************************************************************/
// All sizes in cm.

#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>

void ATTPC_P10_1bar()
{
   const TString geoVersion = "ATTPC_P10_1bar";
   const TString FileName = geoVersion + ".root";
   const TString FileName1 = geoVersion + "_geomanager.root";

   const Float_t tpc_diameter = 50.; // cm
   const Float_t drift_length = 100.; // cm

   auto *geo = new TGeoManager("FAIRGeom", "FAIR geometry");

   // ---- Materials ------------------------------------------------------
   // Vacuum
   auto *matVac = new TGeoMaterial("vacuum4", 1.008, 1, 1e-16);
   auto *medVac = new TGeoMedium("vacuum4", 1, matVac, nullptr);

   // Steel (Fe approx) for the vessel
   auto *matSteel = new TGeoMaterial("steel", 55.85, 26, 7.87);
   auto *medSteel = new TGeoMedium("steel", 2, matSteel, nullptr);

   // Aramid for the window (Kevlar approx: C14 H10 N2 O2 — use weight composition)
   auto *matAramid = new TGeoMixture("aramid", 4, 1.44);
   matAramid->AddElement(12.011, 6, 0.7058); // C
   matAramid->AddElement(1.008, 1, 0.0423);  // H
   matAramid->AddElement(14.007, 7, 0.1176); // N
   matAramid->AddElement(15.999, 8, 0.1343); // O
   auto *medAramid = new TGeoMedium("aramid", 3, matAramid, nullptr);

   // P10 gas: 90% Ar + 10% CH4 by volume → mass fractions Ar 0.9573, C 0.0320, H 0.0107
   // Density at 1 bar, 273 K: 1.654e-3 g/cm^3
   auto *matP10 = new TGeoMixture("P10_1bar", 3, 1.654e-3);
   matP10->AddElement(39.948, 18, 0.9573); // Ar
   matP10->AddElement(12.011, 6, 0.0320);  // C
   matP10->AddElement(1.008, 1, 0.0107);   // H

   // Medium params: id, ifield, fieldm, tmaxfd, stemax, deemax, epsil, stmin
   // Use a 1 mm step ceiling in the gas so MIP pions produce dense MC points.
   Double_t params[10] = {1., 1., 20., 0., 0.1, 0.001, 0.001, 0., 0., 0.};
   auto *medP10 = new TGeoMedium("P10_1bar", 4, matP10, params);

   // ---- Volumes --------------------------------------------------------
   auto *top = new TGeoVolumeAssembly("TOP");
   geo->SetTopVolume(top);

   auto *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(medVac);
   top->AddNode(tpcvac, 1);

   double tpc_rot = 0;

   TGeoVolume *drift_volume =
      geo->MakeTube("drift_volume", medP10, 0., tpc_diameter / 2., drift_length / 2.);
   tpcvac->AddNode(drift_volume, 1,
                   new TGeoCombiTrans(0.0, 0.0, drift_length / 2.0,
                                      new TGeoRotation("drift_volume", 0, tpc_rot, 0)));
   drift_volume->SetTransparency(80);

   TGeoVolume *tpc_window = geo->MakeTube("tpc_window", medAramid, 0., 0.5, 0.00018);
   tpc_window->SetLineColor(kBlue);
   tpcvac->AddNode(tpc_window, 1,
                   new TGeoCombiTrans(0.0, 0.0, 0.0, new TGeoRotation("tpc_window", 0, tpc_rot, 0)));
   tpc_window->SetTransparency(50);

   TGeoVolume *vessel_volume =
      geo->MakeTube("vessel_volume", medSteel, tpc_diameter / 2., (tpc_diameter + 2.0) / 2., drift_length / 2.);
   tpcvac->AddNode(vessel_volume, 1,
                   new TGeoCombiTrans(0.0, 0.0, drift_length / 2.0,
                                      new TGeoRotation("vessel_volume", 0, tpc_rot, 0)));
   vessel_volume->SetTransparency(90);

   std::cout << "Voxelizing." << std::endl;
   top->Voxelize("");
   geo->CloseGeometry();

   TFile *outfile = new TFile(FileName, "RECREATE");
   top->Write();
   outfile->Close();

   TFile *outfile1 = new TFile(FileName1, "RECREATE");
   geo->Write();
   outfile1->Close();

   std::cout << "Wrote " << FileName << " and " << FileName1 << std::endl;
}
