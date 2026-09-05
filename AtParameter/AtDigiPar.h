
#ifndef AtDIGIPAR_H
#define AtDIGIPAR_H

#include <FairParGenericSet.h>

#include <Rtypes.h>

class FairParamList;
class TBuffer;
class TClass;
class TMemberInspector;

class AtDigiPar : public FairParGenericSet {
private:
   Bool_t fInitialized;

   Double_t fBField{};
   Double_t fEField{};

   // Detector geometry
   Int_t fTBEntrance{};
   Double_t fZPadPlane{};

   // Gas properties
   Double_t fEIonize{};       //< effective ionization energy [eV]
   Double_t fFano{};          //< Fano factor of gas
   Double_t fCoefL{};         //< longitudinal diffusion coefficient [cm^2/us]
   Double_t fCoefT{};         //< transversal diffusion coefficient [cm^2/us]
   Double_t fGasPressure{};   //< gas pressure [torr]
   Double_t fDensity{};       //< Gas density [kg/m^3]
   Double_t fDriftVelocity{}; //< Electron drift velocity [cm/us]
   Double_t fGain{};          //< gain factor from wire plane

   // Electronic info
   Int_t fSamplingRate{};
   Double_t fGETGain{};  //< Gain from get electronics in fC
   Int_t fPeakingTime{}; //< Peaking time of the electronics in ns

   /// APPENDED DELIBERATELY, AND ANY FUTURE MEMBER MUST BE TOO. The getters are inline, so a
   /// module compiled against an older header keeps the OLD member offsets; inserting a member
   /// mid-class shifts every member after it and such a module then reads the wrong field with no
   /// error. Appending leaves existing offsets alone, and it is also what ROOT schema evolution
   /// wants for a class with a streamer (AtDigiPar is `+` in AtParLinkDef.h). A full rebuild is
   /// still required after changing this header.
   ///
   /// REVERSED DETECTOR: the beam enters through the PAD PLANE and leaves through the cathode,
   /// so ionisation electrons drift *with* the beam instead of against it. 0 = normal (default).
   ///
   /// This is the ONLY thing that changes. The digi-frame hit z keeps its meaning in both modes
   /// (z_digi = ZPadPlane - z_beam, i.e. 0 at the pad plane), so the hit cloud, the pattern
   /// recognition, the fitters' z_lab = ZPadPlane - z_digi and every downstream convention are
   /// untouched. What differs is the DRIFT LENGTH of a given point: normal drift = ZPadPlane -
   /// z_beam, reversed drift = z_beam. That changes the diffusion each hit picks up and which
   /// time bucket it lands in, and nothing else.
   ///
   /// IT LIVES IN THE PAR, NOT IN A TASK SETTER, ON PURPOSE. Digitisation (AtClusterize) and
   /// reconstruction (AtPSA) must agree about which end the pad plane is on; a setter on each
   /// would let a job digitise reversed and reconstruct normal, which produces a plausible
   /// mirrored z rather than an error. Sharing one par makes that disagreement impossible.
   Int_t fReverseDrift{0};

public:
   // Constructors and Destructors
   AtDigiPar(const Char_t *name, const Char_t *title, const Char_t *context);
   ~AtDigiPar() = default;

   // Getters
   Double_t GetBField() const { return fBField; }
   Double_t GetEField() const { return fEField; }

   Int_t GetTBEntrance() const { return fTBEntrance; }
   Double_t GetZPadPlane() const { return fZPadPlane; }
   /// True when the beam enters through the pad plane (see fReverseDrift).
   Bool_t GetReverseDrift() const { return fReverseDrift != 0; }

   Double_t GetEIonize() const { return fEIonize; }
   Double_t GetFano() const { return fFano; }
   Double_t GetCoefDiffusionTrans() const { return fCoefT; }
   Double_t GetCoefDiffusionLong() const { return fCoefL; }
   Double_t GetGasPressure() const { return fGasPressure; }
   Double_t GetDensity() const { return fDensity; }
   Double_t GetDriftVelocity() const { return fDriftVelocity; }
   Double_t GetGain() const { return fGain; }

   Int_t GetTBTime() const;
   Double_t GetGETGain() const { return fGETGain; };
   Int_t GetPeakingTime() const { return fPeakingTime; };

   // Setters
   virtual void putParams(FairParamList *paramList) override;
   virtual Bool_t getParams(FairParamList *paramList) override;
   // Main methods

   ClassDefOverride(AtDigiPar, 4);
};

#endif
