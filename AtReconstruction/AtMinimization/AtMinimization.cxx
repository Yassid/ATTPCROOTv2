#include "AtMinimization.h"

#include "AtEvent.h"
#include "AtTrack.h"

namespace MCMinimization {

bool AtMinimization::Minimize(const TrackSeed &seed, const AtEvent &event)
{
   return Minimize(seed, event.GetHits());
}

bool AtMinimization::Minimize(const TrackSeed &seed, const AtTrack &track)
{
   return Minimize(seed, track.GetHitArray());
}

} // namespace MCMinimization
