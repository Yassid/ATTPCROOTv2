/********************************************************************************
 *    HYDRA Prototype (GLAD-TPC) geometry — P10 gas at 1 bar.
 *
 *    Port of R3BRootGroup/glad-tpc HYDRAprototype_FileSetup.par into
 *    ATTPCROOT's drift-along-z convention.
 *
 *    HYDRA Prototype params:
 *        TPCL          19.7 × 31.4 × 47.8 cm  (outer envelope, x × y × z)
 *        ActiveRegion  8.8 × 29.4 × 25.6 cm   (HYDRA x × y × z)
 *        PadSize       5 mm
 *        Gas           P10 (90% Ar + 10% CH4)
 *        Target        at HYDRA-(x=-10.825, y=0, z=-91.85) cm — outside
 *                      the chamber, ~80 cm upstream of the entrance face.
 *
 *    Axis mapping (HYDRA → ATTPCROOT, keeping drift along z):
 *        HYDRA-z (beam, 256 mm)        → ATTPCROOT-x  (pad-plane long axis)
 *        HYDRA-x (transverse, 88 mm)   → ATTPCROOT-y  (pad-plane short axis)
 *        HYDRA-y (drift, 294 mm)       → ATTPCROOT-z  (drift)
 *
 *    In this file all sizes are cm.
 ********************************************************************************/

#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TROOT.h"
#include "TString.h"

#include <iostream>

void ATTPC_HYDRA_Prototype_P10_1bar()
{
   const TString geoVersion = "ATTPC_HYDRA_Prototype_P10_1bar";
   const TString FileName = geoVersion + ".root";
   const TString FileName1 = geoVersion + "_geomanager.root";

   // Active drift volume (ATTPCROOT-convention sizes, cm)
   const Double_t active_x = 25.6;  // pad-plane long axis (HYDRA-z = beam)
   const Double_t active_y = 8.8;   // pad-plane transverse (HYDRA-x)
   const Double_t active_z = 29.4;  // drift (HYDRA-y)

   // Outer envelope
   const Double_t outer_x = 47.8;
   const Double_t outer_y = 19.7;
   const Double_t outer_z = 31.4;

   auto *geo = new TGeoManager("FAIRGeom", "FAIR geometry");

   // ---- Materials --------------------------------------------------------
   auto *matVac = new TGeoMaterial("vacuum4", 1.008, 1, 1e-16);
   auto *medVac = new TGeoMedium("vacuum4", 1, matVac, nullptr);

   auto *matSteel = new TGeoMaterial("steel", 55.85, 26, 7.87);
   auto *medSteel = new TGeoMedium("steel", 2, matSteel, nullptr);

   auto *matMylar = new TGeoMixture("mylar", 3, 1.397);
   matMylar->AddElement(12.011, 6, 0.6250);
   matMylar->AddElement(1.008, 1, 0.0420);
   matMylar->AddElement(15.999, 8, 0.3330);
   auto *medMylar = new TGeoMedium("mylar", 3, matMylar, nullptr);

   auto *matP10 = new TGeoMixture("P10_1bar", 3, 1.654e-3);
   matP10->AddElement(39.948, 18, 0.9573);
   matP10->AddElement(12.011, 6, 0.0320);
   matP10->AddElement(1.008, 1, 0.0107);
   Double_t params[10] = {1., 1., 20., 0., 0.1, 0.001, 0.001, 0., 0., 0.};
   auto *medP10 = new TGeoMedium("P10_1bar", 4, matP10, params);

   // ---- Volumes ----------------------------------------------------------
   auto *top = new TGeoVolumeAssembly("TOP");
   geo->SetTopVolume(top);

   auto *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(medVac);
   top->AddNode(tpcvac, 1);

   // Active region lower-left at (x, y, z) = (0, 0, 0). Drift_volume box
   // center is therefore at (active_x/2, active_y/2, active_z/2). Pad
   // plane is at z=0 (lower face), cathode at z=active_z.
   TGeoVolume *drift_volume =
      geo->MakeBox("drift_volume", medP10, active_x / 2., active_y / 2., active_z / 2.);
   tpcvac->AddNode(drift_volume, 1,
                   new TGeoCombiTrans(active_x / 2., active_y / 2., active_z / 2., new TGeoRotation()));
   drift_volume->SetTransparency(80);

   // Hollow steel vessel (composite subtraction so it doesn't overlap gas).
   // Same center as drift_volume so the vessel surrounds the active region.
   const Double_t wall = 1.0;
   auto *outerBox = new TGeoBBox("vessel_outer_solid", outer_x / 2., outer_y / 2., outer_z / 2.);
   auto *innerBox = new TGeoBBox("vessel_inner_void", outer_x / 2. - wall, outer_y / 2. - wall,
                                 outer_z / 2. - wall);
   auto *vesselShape = new TGeoCompositeShape("vessel_shell", "vessel_outer_solid - vessel_inner_void");
   TGeoVolume *vessel_outer = new TGeoVolume("vessel_outer", vesselShape, medSteel);
   tpcvac->AddNode(vessel_outer, 1,
                   new TGeoCombiTrans(active_x / 2., active_y / 2., active_z / 2., new TGeoRotation()));
   vessel_outer->SetTransparency(95);

   // Mylar entrance window on the upstream (-x) face of the vessel.
   // Beam now travels along +x; window face is normal to x.
   TGeoVolume *tpc_window =
      geo->MakeBox("tpc_window", medMylar, 0.005 / 2., active_y / 2. + 1.0, active_z / 2. + 1.0);
   tpc_window->SetLineColor(kBlue);
   tpcvac->AddNode(tpc_window, 1,
                   new TGeoCombiTrans(active_x / 2. - outer_x / 2., active_y / 2.,
                                      active_z / 2., new TGeoRotation()));
   tpc_window->SetTransparency(40);

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
   std::cout << "Active drift volume: " << active_x << " x " << active_y << " x "
             << active_z << " cm (x, y, z)\n";
   std::cout << "Pad plane: " << active_x * 10 << " x " << active_y * 10
             << " mm = 128 x 44 pads at 2 mm pitch (X=beam, Y=transverse)\n";
   std::cout << "Drift length: " << active_z * 10 << " mm along z\n";
}
