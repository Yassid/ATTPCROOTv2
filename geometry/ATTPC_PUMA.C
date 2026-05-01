/********************************************************************************
 *  ATTPCROOT geometry for PUMA-style hollow annular TPC.
 *
 *  Mirrors the relevant volumes of the puma-tpc-simulation
 *  (ExN03DetectorConstruction): central Cu antiproton trap, Cu/Al cryostat
 *  cylinders, and an annular gas drift volume that carries the sensitive
 *  region. Built directly via TGeoMixture / TGeoMedium so this works on
 *  setups where FairGeoLoader hangs reading media.geo.
 *
 *  Dimensions (mm) follow ExN03Setup::ExN03Setup() defaults:
 *    Trap                R = 10               (Cu, length 150)
 *    TrapAl wall         10 < R < 14          (Cu)
 *    Cu2 cylinder        25.25 < R < ~29.25   (Cu)
 *    CryoAl1             ~28.0 < R < 29       (Al)
 *    CryoAl2             ~34.0 < R < 35       (Al)
 *    Gas drift_volume    58.5 < R < 125.1     (P10@1bar, length 300)
 *    World vacuum tube   R = 125.1, length 600
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

void ATTPC_PUMA()
{
   const TString geoVersion = "ATTPC_PUMA";
   const TString FileName = geoVersion + ".root";
   const TString FileName1 = geoVersion + "_geomanager.root";

   // Half-lengths in cm.
   const double trap_hZ = 7.5;     // 150 mm full
   const double drift_hZ = 15.0;   // 300 mm full
   const double world_hZ = 30.0;   // generous

   // Radii in cm.
   const double r_trap = 1.0;
   const double r_trapAl_in = 1.0, r_trapAl_out = 1.4;
   const double r_cu2_in = 2.525, r_cu2_out = 2.925; // 4 mm Cu sleeve
   const double r_cry1_in = 2.8, r_cry1_out = 2.9;
   const double r_cry2_in = 3.4, r_cry2_out = 3.5;
   const double r_drift_in = 5.85, r_drift_out = 12.51;
   const double r_world = 12.51;

   auto *geo = new TGeoManager("FAIRGeom", "PUMA TPC geometry");

   // ---- Materials ------------------------------------------------------
   auto *matVac = new TGeoMaterial("vacuum4", 1.008, 1, 1e-16);
   auto *medVac = new TGeoMedium("vacuum4", 1, matVac, nullptr);

   auto *matCu = new TGeoMaterial("copper", 63.546, 29, 8.96);
   auto *medCu = new TGeoMedium("copper", 2, matCu, nullptr);

   auto *matAl = new TGeoMaterial("aluminum", 26.982, 13, 2.70);
   auto *medAl = new TGeoMedium("aluminum", 3, matAl, nullptr);

   // P10 gas: 90% Ar + 10% CH4 (vol), rho = 1.654e-3 g/cm^3 at 1 bar / 273 K.
   auto *matP10 = new TGeoMixture("P10_1bar", 3, 1.654e-3);
   matP10->AddElement(39.948, 18, 0.9573);
   matP10->AddElement(12.011, 6, 0.0320);
   matP10->AddElement(1.008, 1, 0.0107);
   double medParams[10] = {1., 1., 20., 0., 0.1, 0.001, 0.001, 0., 0., 0.};
   auto *medP10 = new TGeoMedium("P10_1bar", 4, matP10, medParams);

   // ---- Volumes --------------------------------------------------------
   auto *top = new TGeoVolumeAssembly("TOP");
   geo->SetTopVolume(top);

   auto *tpcvac = new TGeoVolumeAssembly(geoVersion);
   tpcvac->SetMedium(medVac);
   top->AddNode(tpcvac, 1);

   // Place all volumes shifted so drift_volume spans z = [0, +2*drift_hZ]
   // (i.e. pad plane at z=0, electrons drift in +z direction toward it).
   // This matches the AT-TPC PSA convention CalculateZGeo expects:
   //   z_recon = ZPadPlane - drift_distance, with pad plane at ZPadPlane.
   const double driftCenterZ = drift_hZ; // shift so drift_volume z in [0, 2*drift_hZ]

   auto add = [&](const char *name, TGeoMedium *med, double rin, double rout, double hz, double zCenter, int color,
                  double vis = 80) {
      auto *vol = geo->MakeTube(name, med, rin, rout, hz);
      vol->SetLineColor(color);
      vol->SetTransparency(vis);
      tpcvac->AddNode(vol, 1, new TGeoTranslation(0., 0., zCenter));
      return vol;
   };

   add("Trap", medCu, 0., r_trap, trap_hZ, driftCenterZ, kOrange + 7, 50);
   add("TrapAl", medCu, r_trapAl_in, r_trapAl_out, trap_hZ, driftCenterZ, kOrange + 9, 60);
   add("Cu2", medCu, r_cu2_in, r_cu2_out, trap_hZ, driftCenterZ, kRed - 7, 70);
   add("CryoAl1", medAl, r_cry1_in, r_cry1_out, drift_hZ, driftCenterZ, kGray, 80);
   add("CryoAl2", medAl, r_cry2_in, r_cry2_out, drift_hZ, driftCenterZ, kGray + 1, 80);

   // Sensitive drift_volume (name matches AtTpc::CheckIfSensitive contains-check).
   TGeoVolume *drift = geo->MakeTube("drift_volume", medP10, r_drift_in, r_drift_out, drift_hZ);
   drift->SetLineColor(kAzure + 1);
   drift->SetTransparency(85);
   tpcvac->AddNode(drift, 1, new TGeoTranslation(0., 0., driftCenterZ));

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
