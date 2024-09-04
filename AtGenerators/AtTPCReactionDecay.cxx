#include "AtTPCReactionDecay.h"

#include "AtFormat.h"
#include "AtVertexPropagator.h"

#include <FairIon.h>
#include <FairLogger.h>
#include <FairParticle.h>
#include <FairPrimaryGenerator.h>
#include <FairRunSim.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <memory>
#include <utility>

constexpr float amu = 931.494;

AtTPCReactionDecay::AtTPCReactionDecay() {}

AtTPCReactionDecay::AtTPCReactionDecay(DecayIon &scatter, DecayIon &recoil, Int_t zBeam, Int_t aBeam, Double_t massBeam,
                                       Double_t massTarget)
   : fScatter(scatter), fRecoil(recoil), fZBeam(zBeam), fABeam(aBeam), fBeamMass(massBeam), fTargetMass(massTarget),
     fEnergy(0), fIon(0), fParticle(0), fPType(0), fVx(0.), fVy(0.), fVz(0.), fNumIons(0), fNumTracks(0),
     fIsDecay(kFALSE)
{


   RegisterIons(fScatter);
   RegisterIons(fRecoil);
}

Bool_t AtTPCReactionDecay::GenerateReaction(FairPrimaryGenerator *primGen)
{

   LOG(info) << cYELLOW << " =======  AtTPCReactionDecay - Generating decay  ========= " << cNORMAL << "\n";

   fNumTracks = 0;

   fScatter.exEnergy = AtVertexPropagator::Instance()->GetScatterEx() / 1000.0;
   fRecoil.exEnergy = AtVertexPropagator::Instance()->GetRecoilEx() / 1000.0;
   fScatter.momentum = AtVertexPropagator::Instance()->GetScatterP();
   fRecoil.momentum = AtVertexPropagator::Instance()->GetRecoilP();

   fVx = AtVertexPropagator::Instance()->GetVx();
   fVy = AtVertexPropagator::Instance()->GetVy();
   fVz = AtVertexPropagator::Instance()->GetVz();

   // auto scatterDecay = GetSimultaneousDecay(fScatter);
   // auto recoilDecay = GetSimultaneousDecay(fRecoil);

   std::vector<TLorentzVector> scatterDecay;
   std::vector<TLorentzVector> recoilDecay;

   if (fScatter.sequentialDecay) {
      scatterDecay = GetSequentialDecay(fScatter);
   } else {
      scatterDecay = GetSimultaneousDecay(fScatter);
   }

   if (fRecoil.sequentialDecay) {
      recoilDecay = GetSequentialDecay(fRecoil);
   } else {
      recoilDecay = GetSimultaneousDecay(fRecoil);
   }

   PropagateIons(scatterDecay, primGen);
   PropagateIons(recoilDecay, primGen);

   return kTRUE;
}

std::vector<TLorentzVector> AtTPCReactionDecay::GetSequentialDecay(DecayIon &ion)
{

   LOG(info) << cYELLOW << "           Sequential decay. track ID : " << ion.trackID << cNORMAL << "\n";

   // 2 step sequential decay
   Double_t beta;
   Double_t s = 0.0;
   Double_t mass[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   Double_t M_tot = 0;
   TLorentzVector fEnergyImpulsionLab_beam;
   TLorentzVector fEnergyImpulsionLab_Total;
   TLorentzVector fEnergyImpulsionLab_target;
   TLorentzVector fEnergyImpulsionFinal;

   TVector3 fImpulsionLab_beam;
   std::vector<TLorentzVector> phaseSpaceMomentum;
   TGenPhaseSpace eventParent;
   TGenPhaseSpace eventDaughter;

   TParticlePDG *thisPart0 = nullptr;
   thisPart0 = TDatabasePDG::Instance()->GetParticle(
      fIon.at(0)->GetName()); // NB: The first particle of the list must be the decaying ion
   int pdgType0 = thisPart0->PdgCode();

   LOG(info) << " Ejectile info : " << pdgType0 << " " << fBeamMass << " " << ion.momentum.X() << " "
             << ion.momentum.Y() << " " << ion.momentum.Z() << " " << fEnergy << " " << ion.exEnergy << "\n";

   fEnergy = AtVertexPropagator::Instance()->GetTrackEnergy(ion.trackID) / 1000.0; // 0 - scattered, 1-recoil
   fBeamMass = ion.parentMass;

   fImpulsionLab_beam = TVector3(ion.momentum.X(), ion.momentum.Y(), ion.momentum.Z());
   fEnergyImpulsionLab_beam = TLorentzVector(fImpulsionLab_beam, fBeamMass + fEnergy + ion.exEnergy);
   fEnergyImpulsionLab_target = TLorentzVector(TVector3(0, 0, 0), fTargetMass);

   fEnergyImpulsionLab_Total = fEnergyImpulsionLab_beam + fEnergyImpulsionLab_target;

   s = fEnergyImpulsionLab_Total.M2();

   for (auto i = 0; i < ion.parentMultiplicity; i++) {
      M_tot += ion.mass.at(i) * amu / 1000.0;
      mass[i] = ion.mass.at(i) * amu / 1000.0;
      LOG(info) << cYELLOW << "           - Parent ion :  " << i << ".   " << ion.mass.at(i) * amu / 1000.0 << " "
                << M_tot << " " << fBeamMass << " " << fBeamMass - M_tot << " " << ion.exEnergy << " " << sqrt(s) << " "
                << (sqrt(s) - M_tot) * 1000.0 << cNORMAL << std::endl;
   }

   if (eventParent.SetDecay(fEnergyImpulsionLab_Total, ion.parentMultiplicity, mass)) {
      eventParent.Generate();

      // NB Assume the first particle is the decaying ion daugther
      for (Int_t i = 1; i < ion.parentMultiplicity; i++) {
         phaseSpaceMomentum.push_back(
            *eventParent.GetDecay(i)); // Propagate particles that won sequentially decay. This can be controlled with
                                       // ion.decays.at(i) (Perhaps a better idea)
      }

      // Reset mass
      for (auto i = 0; i < 10; i++) {
         mass[i] = 0.0;
      }

      for (auto i = 0; i < ion.daughterMultiplicity; i++) {
         mass[i] = ion.mass.at(i + ion.parentMultiplicity) * amu / 1000.0;
         LOG(info) << cYELLOW << "           - Daughter ion :  " << i + ion.parentMultiplicity << ".   "
                   << ion.mass.at(i + ion.parentMultiplicity) * amu / 1000.0 << " " << cNORMAL << std::endl;
      }

      if (eventDaughter.SetDecay(*eventParent.GetDecay(0), ion.daughterMultiplicity, mass)) {
         eventDaughter.Generate();
         for (Int_t i = 0; i < ion.daughterMultiplicity; i++) {
            phaseSpaceMomentum.push_back(*eventDaughter.GetDecay(i)); // Same as before.
         }
      } else {
         LOG(info) << cYELLOW << "AtTPCIonDecay - Warning, kinematical conditions for daughter decay not fulfilled "
                   << cNORMAL << "\n";
      }

   } else {
      LOG(info) << cYELLOW << "AtTPCIonDecay - Warning, kinematical conditions for parent decay not fulfilled "
                << cNORMAL << "\n";
      LOG(info) << cYELLOW << " s = " << s << " - pow(M_tot,2) = " << pow(M_tot, 2) << cNORMAL << "\n";
   }

   return phaseSpaceMomentum;
}

std::vector<TLorentzVector> AtTPCReactionDecay::GetSimultaneousDecay(DecayIon &ion)
{

   Double_t beta;
   Double_t s = 0.0;
   Double_t mass[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   Double_t M_tot = 0;
   TLorentzVector fEnergyImpulsionLab_beam;
   TLorentzVector fEnergyImpulsionLab_Total;
   TLorentzVector fEnergyImpulsionLab_target;
   TLorentzVector fEnergyImpulsionFinal;

   TVector3 fImpulsionLab_beam;
   std::vector<TLorentzVector> phaseSpaceMomentum;
   TGenPhaseSpace event1;

   TParticlePDG *thisPart0 = nullptr;
   thisPart0 = TDatabasePDG::Instance()->GetParticle(fIon.at(0)->GetName());
   int pdgType0 = thisPart0->PdgCode();

   LOG(info) << " Ejectile info : " << pdgType0 << " " << fBeamMass << " " << ion.momentum.X() << " "
             << ion.momentum.Y() << " " << ion.momentum.Z() << " " << fEnergy << " " << ion.exEnergy << "\n";

   fEnergy = AtVertexPropagator::Instance()->GetTrackEnergy(ion.trackID) / 1000.0; // 0 - scattered, 1-recoil
   fBeamMass = ion.parentMass;

   fImpulsionLab_beam = TVector3(ion.momentum.X(), ion.momentum.Y(), ion.momentum.Z());
   fEnergyImpulsionLab_beam = TLorentzVector(fImpulsionLab_beam, fBeamMass + fEnergy + ion.exEnergy);
   fEnergyImpulsionLab_target =
      TLorentzVector(TVector3(0, 0, 0), fTargetMass); // NB: fTargetMass should be always 0 for decay

   fEnergyImpulsionLab_Total = fEnergyImpulsionLab_beam + fEnergyImpulsionLab_target;

   s = fEnergyImpulsionLab_Total.M2();
   // beta = fEnergyImpulsionLab_Total.Beta();

   for (auto i = 0; i < ion.multiplicity; i++) {
      M_tot += ion.mass.at(i) * amu / 1000.0;
      mass[i] = ion.mass.at(i) * amu / 1000.0;
      std::cout << ion.mass.at(i) * amu / 1000.0 << " " << M_tot << " " << fBeamMass << " " << fBeamMass - M_tot << " "
                << ion.exEnergy << " " << sqrt(s) << " " << (sqrt(s) - M_tot) * 1000.0 << std::endl;
   }

   if (s > pow(M_tot, 2)) {
      fIsDecay = kTRUE;
      event1.SetDecay(fEnergyImpulsionLab_Total, ion.multiplicity, mass);
      event1.Generate(); // to generate the random final state
      std::vector<Double_t> kineticEnergy;
      std::vector<Double_t> thetaLab;

      LOG(info) << " AtTPCReactionDecay -  Phase Space Information "
                << "\n";
      for (Int_t i = 0; i < ion.multiplicity; i++) {
         phaseSpaceMomentum.push_back(*event1.GetDecay(i));

         kineticEnergy.push_back((phaseSpaceMomentum.at(i).E() - mass[i]) * 1000.0);
         thetaLab.push_back(phaseSpaceMomentum.at(i).Theta() * 180. / TMath::Pi());

         LOG(info) << " Particle " << i << " - TKE (MeV) : " << kineticEnergy.at(i)
                   << " - Lab Angle (deg) : " << thetaLab.at(i) << "\n";
      }
      LOG(info) << cNORMAL;

   } else { // if kinematics condition

      LOG(info) << cYELLOW << "AtTPCIonDecay - Warning, kinematical conditions for decay not fulfilled " << cNORMAL
                << "\n";
      LOG(info) << cYELLOW << " s = " << s << " - pow(M_tot,2) = " << pow(M_tot, 2) << cNORMAL << "\n";
   }

   return phaseSpaceMomentum;
}

bool AtTPCReactionDecay::PropagateIons(std::vector<TLorentzVector> phaseSpaceMomentum, FairPrimaryGenerator *primGen)
{

   Int_t numTracks = 0;

   for (Int_t i = 0; i < phaseSpaceMomentum.size(); i++) {

      TParticlePDG *thisPart = nullptr;

      if (fPType.at(i + fNumTracks) == "Ion")
         thisPart = TDatabasePDG::Instance()->GetParticle(fIon.at(i + fNumTracks)->GetName());
      else if (fPType.at(i + fNumTracks) == "Proton")
         thisPart = TDatabasePDG::Instance()->GetParticle(fParticle.at(i + fNumTracks)->GetName());
      else if (fPType.at(i + fNumTracks) == "Neutron")
         thisPart = TDatabasePDG::Instance()->GetParticle(fParticle.at(i + fNumTracks)->GetName());

      if (!thisPart) {
         if (fPType.at(i + fNumTracks) == "Ion")
            std::cout << "-W- FairIonGenerator: Ion " << fIon.at(i + fNumTracks)->GetName() << " not found in database!"
                      << std::endl;
         else if (fPType.at(i + fNumTracks) == "Proton")
            std::cout << "-W- FairIonGenerator: Particle " << fParticle.at(i + fNumTracks)->GetName()
                      << " not found in database!" << std::endl;
         else if (fPType.at(i + fNumTracks) == "Neutron")
            std::cout << "-W- FairIonGenerator: Particle " << fParticle.at(i + fNumTracks)->GetName()
                      << " not found in database!" << std::endl;
         return kFALSE;
      }

      int pdgType = thisPart->PdgCode();

      auto momentum = phaseSpaceMomentum.at(i);

      LOG(info) << cGREEN << "       -I- FairIonGenerator: Generating "
                << " with mass " << thisPart->Mass() << " ions of type " << fIon.at(i + fNumTracks)->GetName()
                << " (PDG code " << pdgType << ")" << cNORMAL << std::endl;
      // std::cout << "    Momentum(" << fPx.at(i) << ", " << fPy.at(i) << ", " << fPz.at(i)<< ") Gev from vertex (" <<
      // fVx << ", " << fVy<< ", " << fVz << ") cm" << std::endl;

      primGen->AddTrack(pdgType, momentum.Px(), momentum.Py(), momentum.Pz(), fVx, fVy, fVz);
      ++numTracks;
   }

   fNumTracks += numTracks;

   return kTRUE;
}

void AtTPCReactionDecay::RegisterIons(DecayIon &ion)
{

   FairRunSim *run = FairRunSim::Instance();
   if (!run) {
      std::cout << "-E- FairIonGenerator: No FairRun instantised!" << std::endl;
      Fatal("FairIonGenerator", "No FairRun instantised!");
      return;
   }

   for (Int_t i = 0; i < ion.multiplicity; ++i) {

      if (ion.decays.at(i) == 1) {
         LOG(info) << cYELLOW << "Ion " << i
                   << " will decay and won't be registered for propagation. Mass: " << ion.mass.at(i) * amu / 1000.0
                   << cNORMAL << std::endl;
         continue;
      }

      std::unique_ptr<FairIon> IonBuff = nullptr;
      std::unique_ptr<FairParticle> ParticleBuff = nullptr;
      char buffer[40];
      sprintf(buffer, "Product_Ion_%d", i);
      if (ion.a.at(i) != 1) {
         IonBuff = std::make_unique<FairIon>(buffer, ion.z.at(i), ion.a.at(i), ion.q.at(i), 0.0,
                                             ion.mass.at(i) * amu / 1000.0);
         ParticleBuff = std::make_unique<FairParticle>("dummyPart", 1, 1, 1.0, 0, 0.0, 0.0);
         fPType.push_back("Ion");
         run->AddNewIon(IonBuff.get());

      } else if (ion.a.at(i) == 1 && ion.z.at(i) == 1) {
         IonBuff = std::make_unique<FairIon>(buffer, ion.z.at(i), ion.a.at(i), ion.q.at(i), 0.0,
                                             ion.mass.at(i) * amu / 1000.0);
         auto *kProton = new TParticle(); // NOLINT
         kProton->SetPdgCode(2212);
         ParticleBuff = std::make_unique<FairParticle>(2212, kProton);
         fPType.push_back("Proton");

      } else if (ion.a.at(i) == 1 && ion.z.at(i) == 0) {
         IonBuff = std::make_unique<FairIon>(buffer, ion.z.at(i), ion.a.at(i), ion.q.at(i), 0.0,
                                             ion.mass.at(i) * amu / 1000.0);
         auto *kNeutron = new TParticle(); // NOLINT
         kNeutron->SetPdgCode(2112);
         ParticleBuff = std::make_unique<FairParticle>(2112, kNeutron);
         fPType.push_back("Neutron");
      }

      // Filter out the ions that will not decay
      fIon.push_back(std::move(IonBuff));
      fParticle.push_back(std::move(ParticleBuff));
      fNumIons++;
   }

   std::cout << "Number of Ions registered: " << fNumIons << " " << fIon.size() << std::endl;
}

ClassImp(AtTPCReactionDecay)