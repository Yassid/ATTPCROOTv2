// Closure test for MCMinimization::AtMCQMinimization (AtReconstruction/AtMinimization).
//
//   1. generate a track with the minimizer itself (a single trial, so the truth is known),
//   2. turn the simulated pads into AtHits ("data"),
//   3. fit those hits starting from a perturbed seed and compare with the truth.
//
// It is a standalone program rather than a unit test because it needs a parameter file and a pad
// plane map. Build and run it with (after sourcing build/config.sh):
//
//   O=$VMCWORKDIR
//   FAIRROOT=~/fair_install/FairRootInstall      # FAIRROOTPATH, config.sh does not export it
//   FAIRSOFT=~/fair_install/FairSoft/install     # SIMPATH (fairlogger, boost)
//   g++ -O2 -std=c++17 -o mcqClosure MCQMinimizationClosure.cxx \
//       $O/AtReconstruction/AtMinimization/AtMinimization.cxx \
//       $O/AtReconstruction/AtMinimization/AtMCQMinimization.cxx \
//       -I$O/AtReconstruction/AtMinimization -I$O/AtData -I$O/AtData/AtPattern -I$O/AtMap \
//       -I$O/AtParameter -I$O/AtTools \
//       -I$O/build/_deps/catima-src/include -I$O/build/_deps/catima-build/generated \
//       -isystem $FAIRROOT/include -isystem $FAIRSOFT/include -isystem $FAIRSOFT/include/vmc \
//       -L$O/build/lib -lAtData -lAtMap -lAtParameter -lAtTools -lcatima \
//       -L$FAIRROOT/lib -lBase -lParBase -lFairTools -lFairLogger \
//       $(root-config --cflags --libs) -lGenVector -lXMLParser
//   ./mcqClosure 50 1     # radius of the track [mm], use the position chi2 (0/1)
//
// Reference result of `./mcqClosure 50 1` (proton in D2 at 300 torr, 2.85 T): starting 4.55 deg
// and 6.2 mm away from the truth, the fit lands within 0.07 deg and 0.07 mm of it.

#include "AtDigiPar.h"
#include "AtELossCATIMA.h"
#include "AtEvent.h"
#include "AtHit.h"
#include "AtMCQMinimization.h"
#include "AtTpcMap.h"

#include <FairParAsciiFileIo.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <memory>

using MCMinimization::AtMCQMinimization;

namespace {
const char *kRepo = gSystem->Getenv("VMCWORKDIR");

std::shared_ptr<AtTools::AtELossCATIMA> MakeELoss()
{
   // D2 at 300 torr, 293 K -> 6.61e-5 g/cm3 (the "Density: 0.0661 mg/cm3" of the parameter file)
   auto eLoss = std::make_shared<AtTools::AtELossCATIMA>(6.61e-5);
   eLoss->SetProjectile(1, 1, 938.272 / 931.494);
   std::vector<std::tuple<int, int, int>> material;
   material.emplace_back(2, 1, 2); // D2
   eLoss->SetMaterial(material);
   return eLoss;
}

void Configure(AtMCQMinimization &min, const std::shared_ptr<AtTpcMap> &map,
               const std::shared_ptr<AtTools::AtELossCATIMA> &eLoss, bool usePosChi2)
{
   min.SetMap(map);
   min.SetParticle(1, 1);
   min.SetELossModel(eLoss);
   min.SetUsePosChi2(usePosChi2);
   min.SetBackwardPropagation(false);
}
} // namespace

int main(int argc, char **argv)
{
   const double radius = argc > 1 ? std::atof(argv[1]) : 50.;
   const bool usePos = argc > 2 && std::atoi(argv[2]) != 0;
   gRandom->SetSeed(20260817);

   auto *run = new FairRunAna();
   auto *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open(TString(kRepo) + "/parameters/ATTPC.a1975_deuterium.par", "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");
   rtdb->initContainers(0);

   auto map = std::make_shared<AtTpcMap>();
   map->ParseXMLMap(TString(kRepo) + "/scripts/Lookup20150611.xml");
   map->GeneratePadPlane();
   std::cout << "Map: " << map->GetNumPads() << " pads" << std::endl;

   auto eLoss = MakeELoss();
   std::cout << "CATIMA: dE/dx(3 MeV) = " << eLoss->GetdEdx(3.) << " MeV/mm, range(3 MeV) = " << eLoss->GetRange(3.)
             << " mm" << std::endl;

   /*** 1. Generate the track to fit ***/
   AtMCQMinimization generator;
   Configure(generator, map, eLoss, false);
   generator.SetNumIterations(1, 1); // a single trial around the seed
   generator.SetVerbose(false);

   AtMCQMinimization::TrackSeed truthSeed;
   truthSeed.fVertex = AtMCQMinimization::XYZPoint(0., 0., 0.);
   truthSeed.fVertexTB = 280;
   truthSeed.fTheta = 60. * TMath::DegToRad();
   truthSeed.fPhi = 30. * TMath::DegToRad();
   truthSeed.fRadius = radius; // mm

   AtEvent dummy;
   dummy.AddHit(0, AtHit::XYZPoint(0, 0, 0), 0.);
   if (!generator.Minimize(truthSeed, dummy)) {
      std::cout << "FAILED: could not generate the track" << std::endl;
      return 1;
   }

   const auto truth = generator.GetFitPar();
   const auto &pads = generator.GetSimPads();
   const auto &charges = generator.GetSimCharge();
   std::cout << "Truth: theta " << truth.fTheta * TMath::RadToDeg() << " deg, phi " << truth.fPhi * TMath::RadToDeg()
             << " deg, radius " << truth.fRadius << " mm, energy " << truth.fEnergy << " MeV/u, start ("
             << truth.fPos.X() << ", " << truth.fPos.Y() << ", " << truth.fPos.Z() << ") cm" << std::endl;
   std::cout << "Track: " << generator.GetSimTrack().size() << " points, " << pads.size() << " pads over threshold"
             << std::endl;

   if (pads.size() < 20) {
      std::cout << "FAILED: the generated track is too short to fit" << std::endl;
      return 1;
   }

   /*** 2. The generated pads become the experimental hits ***/
   AtEvent event;
   const double gain = generator.GetGain(); // GetSimCharge() is in primary electrons
   const double tbLength = 1.8176;          // drift velocity 1.136 cm/us x 160 ns, in mm
   const double zStart = generator.GetSimTrack().front().Z();
   for (size_t i = 0; i < pads.size(); i++) {
      const int padNum = map->GetPadNum(ROOT::Math::XYPoint(pads[i].X(), pads[i].Y()));
      if (padNum < 0)
         continue;
      // Time bucket of the hit, with the same convention the minimizer uses to sort the hits
      const int tb = truthSeed.fVertexTB - static_cast<int>(std::round((zStart - pads[i].Z()) / tbLength));
      auto &hit = event.AddHit(padNum, AtHit::XYZPoint(pads[i].X(), pads[i].Y(), pads[i].Z()), charges[i] * gain);
      hit.SetTimeStamp(tb);
   }
   std::cout << "Data: " << event.GetNumHits() << " hits" << std::endl;

   /*** 3. Fit them starting away from the truth ***/
   AtMCQMinimization fitter;
   Configure(fitter, map, eLoss, usePos);
   fitter.SetNumIterations(6, 60);

   AtMCQMinimization::TrackSeed seed = truthSeed;
   seed.fTheta += 4. * TMath::DegToRad();
   seed.fPhi -= 4. * TMath::DegToRad();
   seed.fRadius *= 1.08;
   seed.fVertex = AtMCQMinimization::XYZPoint(4., -4., 0.);
   seed.fNumExpPoints = event.GetNumHits();

   if (!fitter.Minimize(seed, event)) {
      std::cout << "FAILED: the fit did not run" << std::endl;
      return 1;
   }

   const auto fit = fitter.GetFitPar();
   const double dTheta = (fit.fTheta - truth.fTheta) * TMath::RadToDeg();
   const double dPhi = (fit.fPhi - truth.fPhi) * TMath::RadToDeg();
   const double dRadius = fit.fRadius - truth.fRadius;
   const double dEnergy = fit.fEnergy - truth.fEnergy;
   const double dPos =
      std::sqrt(std::pow(fit.fPos.X() - truth.fPos.X(), 2) + std::pow(fit.fPos.Y() - truth.fPos.Y(), 2) +
                std::pow(fit.fPos.Z() - truth.fPos.Z(), 2));

   const double seedTheta = (seed.fTheta - truth.fTheta) * TMath::RadToDeg();

   std::cout << "\n--- closure ---" << std::endl;
   std::cout << "  seed  - truth: theta " << seedTheta << " deg, radius " << (seed.fRadius - truth.fRadius) << " mm"
             << std::endl;
   std::cout << "  fit   - truth: theta " << dTheta << " deg, phi " << dPhi << " deg, radius " << dRadius
             << " mm, energy " << dEnergy << " MeV/u, start " << dPos << " cm" << std::endl;
   std::cout << "  chi2: " << fit.fChi2 << " (charge " << fit.fChi2Q << ")" << std::endl;

   const bool ok = std::abs(dTheta) < std::abs(seedTheta) && std::abs(dRadius) < std::abs(seed.fRadius - truth.fRadius);
   std::cout << (ok ? "IMPROVED" : "NOT IMPROVED") << std::endl;
   return ok ? 0 : 1;
}
