// -------------------------------------------------------------------------
// -----                   AtTPCIonGenerator header file               -----
// -----                    Based on FairIonGenerator                  -----
// -----                 Created 30/01/15  by Y. Ayyad                 -----
// -------------------------------------------------------------------------

/** Include beam spot size **/

#ifndef AtTPCIONGENERAtOR_H
#define AtTPCIONGENERAtOR_H

#include <FairGenerator.h>

#include <Rtypes.h>

class FairPrimaryGenerator;
class FairIon;
class TBuffer;
class TClass;
class TMemberInspector;

class AtTPCIonGenerator : public FairGenerator {
protected:
   static Int_t fgNIon;    //! Number of the instance of this class
   Int_t fMult;            //< Multiplicity per event
   Double_t fPx, fPy, fPz; //< Momentum components [GeV] per nucleon

   Double32_t fR, fz, fOffsetX, fOffsetY; //< beam Spot radius [cm], z source, y source, x source
   Double_t fVx, fVy, fVz;                //< Vertex coordinates [cm]

   FairIon *fIon; //< Pointer to the FairIon to be generated
   Int_t fQ;      //< Electric charge [e]
   Int_t fNomEner{};
   Double_t fMaxEnLoss{}; //< Max energy loss before reation happens

   Bool_t fDoReact{true};

   /// Sets fVx, fVy, fVz depending on the type of ion generator.
   virtual void SetVertexCoordinates();

public:
   /** Default constructor **/
   AtTPCIonGenerator();

   /** Constructor with ion name
    ** For the generation of ions with pre-defined FairIon
    ** By default, the  excitation energy is zero. This can be changed with the
    ** respective modifiers.
    **@param ionName  Ion name
    **@param mult      Number of ions per event
    **@param px,py,pz  Momentum components [GeV] per nucleon!
    **@param vx,vy,vz  Vertex coordinates [cm]
    **/
   AtTPCIonGenerator(const Char_t *ionName, Int_t mult, Double_t px, Double_t py, Double_t pz);

   /** Default constructor
    ** For the generation of ions with atomic number z and mass number a.
    ** By default, the mass equals a times the proton mass and the
    ** excitation energy is zero. This can be changed with the
    ** respective modifiers.
    **@param z         Atomic number
    **@param a         Atomic mass
    **@param q         Electric charge [e]
    **@param mult      Number of ions per event
    **@param px,py,pz  Momentum components [GeV] per nucleon!
    **@param vx,vy,vz  Vertex coordinates [cm]
    **@param ener      Energy of the ion.
    **@param eLoss     Maximum energy loss before reaction happens. Defaults to ener.
    **/
   AtTPCIonGenerator(const char *name, Int_t z, Int_t a, Int_t q, Int_t mult, Double_t px, Double_t py, Double_t pz,
                     Double_t Ex, Double_t m, Double_t ener, Double_t eLoss = -1);

   /** Destructor **/
   virtual ~AtTPCIonGenerator() = default;

   /** Modifiers **/
   void SetCharge(Int_t charge) { fQ = charge; }
   void SetExcitationEnergy(Double_t eExc);
   void SetMass(Double_t mass);

   void SetSpotRadius(Double32_t r = 0, Double32_t z = 0, Double32_t offx = 0, Double32_t offy = 0);

   void SetDoReaction(Bool_t doReact) { fDoReact = doReact; }

   /** Method ReadEvent
   ** Generates <mult> of the specified ions and hands hem to the
   ** FairPrimaryGenerator.
   **/
   virtual Bool_t ReadEvent(FairPrimaryGenerator *primGen);

   ClassDef(AtTPCIonGenerator, 2)
};

#endif
