#ifndef AtFITTER_H
#define AtFITTER_H

#include "AtDigiPar.h"
#include "AtTrack.h"
#include "AtEvent.h"
#include "AtPatternEvent.h"

// FairRoot classes
#include "FairRootManager.h"
#include "FairLogger.h"

// GENFIT
#include "Track.h"
#include "TrackCand.h"
#include "RKTrackRep.h"
#include "Exception.h"

#define cRED "\033[1;31m"
#define cYELLOW "\033[1;33m"
#define cNORMAL "\033[0m"
#define cGREEN "\033[1;32m"

namespace AtFITTER {

class AtFitter : public TObject {

public:
   AtFitter();
   virtual ~AtFitter();
   // virtual std::vector<AtTrack> GetFittedTrack() = 0;
   virtual genfit::Track *FitTracks(AtTrack *track) = 0;
   virtual void Init() = 0;

   void MergeTracks(std::vector<AtTrack> *trackCandSource, std::vector<AtTrack> *trackJunkSource,
                    std::vector<AtTrack> *trackDest, bool fitDirection, bool simulationConv);
   Bool_t MergeTracks(std::vector<AtTrack *> *trackCandSource, std::vector<AtTrack> *trackDest,
                      Bool_t enableSingleVertexTrack);
   void ClusterizeSmooth3D(AtTrack &track, Float_t distance, Float_t radius);

protected:
   FairLogger *fLogger; ///< logger pointer
   AtDigiPar *fPar;     ///< parameter container
   std::tuple<Double_t, Double_t>
   GetMomFromBrho(Double_t A, Double_t Z,
                  Double_t brho); ///< Returns momentum (in GeV) from Brho assuming M (amu) and Z;
   Bool_t FindVertexTrack(AtTrack *trA, AtTrack *trB); ///< Lambda function to find track closer to vertex
   ClassDef(AtFitter, 1)
};

} // namespace AtFITTER

#endif
