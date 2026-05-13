/********************************************************************************
 *    HYDRA (GLAD-TPC) geometry — P10 gas at 1 bar.
 *
 *    Port of the R3BRoot/glad-tpc detector (HYDRAFullBeamIn) into ATTPCROOT's
 *    drift-along-z convention. Parameters from
 *    https://github.com/R3BRootGroup/glad-tpc, params/HYDRAFullBeamIn_FileSetup.par:
 *        ActiveRegion 7.2 × 15 × 100 cm (HYDRA x × y × z, drift along y)
 *        Pad pitch    5 mm
 *        Gas          P10 (90% Ar + 10% CH4)
 *        Outer box    36 × 20 × 120 cm
 *
 *    Axis mapping (HYDRA → ATTPCROOT):
 *        HYDRA-x (transverse, 72 mm)   → ATTPCROOT-x
 *        HYDRA-z (along beam, 1000 mm) → ATTPCROOT-y
 *        HYDRA-y (drift, 150 mm)       → ATTPCROOT-z
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

void ATTPC_HYDRA_P10_1bar()
{
   const TString geoVersion = "ATTPC_HYDRA_P10_1bar";
   const TString FileName = geoVersion + ".root";
   const TString FileName1 = geoVersion + "_geomanager.root";

   // ---- Active drift volume (ATTPCROOT-convention sizes) -----------------
   const Double_t active_x = 7.2;   // cm — pad-plane horizontal
   const Double_t active_y = 100.0; // cm — pad-plane "long axis" (along beam)
   const Double_t active_z = 15.0;  // cm — drift direction

   // Outer gas envelope (steel vessel inscribes this rectangular volume).
   const Double_t outer_x = 36.0;
   const Double_t outer_y = 120.0;
   const Double_t outer_z = 20.0;

   auto *geo = new TGeoManager("FAIRGeom", "FAIR geometry");

   // ---- Materials --------------------------------------------------------
   auto *matVac = new TGeoMaterial("vacuum4", 1.008, 1, 1e-16);
   auto *medVac = new TGeoMedium("vacuum4", 1, matVac, nullptr);

   auto *matSteel = new TGeoMaterial("steel", 55.85, 26, 7.87);
   auto *medSteel = new TGeoMedium("steel", 2, matSteel, nullptr);

   // Mylar window: C10 H8 O4 — mass fractions per glad-tpc media_tpc.geo.
   auto *matMylar = new TGeoMixture("mylar", 3, 1.397);
   matMylar->AddElement(12.011, 6, 0.6250); // C
   matMylar->AddElement(1.008, 1, 0.0420);  // H
   matMylar->AddElement(15.999, 8, 0.3330); // O
   auto *medMylar = new TGeoMedium("mylar", 3, matMylar, nullptr);

   // P10 at 1 bar, 273 K → 1.654e-3 g/cm³. Mass fractions 0.9573 Ar / 0.0320 C / 0.0107 H.
   auto *matP10 = new TGeoMixture("P10_1bar", 3, 1.654e-3);
   matP10->AddElement(39.948, 18, 0.9573);
   matP10->AddElement(12.011, 6, 0.0320);
   matP10->AddElement(1.008, 1, 0.0107);

   // Medium params: 1 mm step ceiling in gas for dense MC points.
   Double_t params[10] = {1., 1., 20., 0., 0.1, 0.001, 0.001, 0., 0., 0.};
   auto *medP10 = new TGeoMedium("P10_1bar", 4, matP10, params);

   // ---- Volumes ----------------------------------------------------------
   auto *top = new TGeoVolumeAssembly("TOP");
   geo->SetTopVolume(top);

   auto *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(medVac);
   top->AddNode(tpcvac, 1);

   // Drift volume: rectangular box. Centered at (0, 0, active_z/2) so the
   // pad plane sits at z=0 and the cathode at z = active_z. Same convention
   // as ATTPC_P10_1bar.
   TGeoVolume *drift_volume =
      geo->MakeBox("drift_volume", medP10, active_x / 2., active_y / 2., active_z / 2.);
   tpcvac->AddNode(drift_volume, 1,
                   new TGeoCombiTrans(0., 0., active_z / 2., new TGeoRotation()));
   drift_volume->SetTransparency(80);

   // Vessel: hollow steel box around the active volume. Built as a composite
   // subtraction (outer minus inner) so it does NOT overlap the drift_volume
   // gas — Geant4/TGeo's "last placed wins" would otherwise turn the gas
   // into steel for any volume that's added after drift_volume.
   const Double_t wall = 1.0; // cm shell thickness
   auto *outerBox = new TGeoBBox("vessel_outer_solid", outer_x / 2., outer_y / 2., outer_z / 2.);
   auto *innerBox = new TGeoBBox("vessel_inner_void", outer_x / 2. - wall, outer_y / 2. - wall,
                                 outer_z / 2. - wall);
   auto *vesselShape = new TGeoCompositeShape("vessel_shell", "vessel_outer_solid - vessel_inner_void");
   TGeoVolume *vessel_outer = new TGeoVolume("vessel_outer", vesselShape, medSteel);
   tpcvac->AddNode(vessel_outer, 1,
                   new TGeoCombiTrans(0., 0., active_z / 2., new TGeoRotation()));
   vessel_outer->SetTransparency(95);

   // Mylar entrance window (upstream end of the beam axis = -y face). HYDRA
   // window dimensions: 15 × 5 × 0.005 cm in HYDRA-(x, y, z); under the
   // axis mapping the window face is normal to ATTPCROOT-y and has size
   // (x, z) = (15 cm, 5 cm). Thickness 50 µm.
   TGeoVolume *tpc_window =
      geo->MakeBox("tpc_window", medMylar, 15.0 / 2., 0.005 / 2., 5.0 / 2.);
   tpc_window->SetLineColor(kBlue);
   tpcvac->AddNode(tpc_window, 1,
                   new TGeoCombiTrans(0., -outer_y / 2., active_z / 2., new TGeoRotation()));
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
             << " mm = 14 x 200 pads at 5 mm pitch (use AtTpcSquareMap)\n";
   std::cout << "Drift length: " << active_z * 10 << " mm along z\n";
}
