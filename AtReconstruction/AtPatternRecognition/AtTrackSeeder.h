#ifndef ATTRACKSEEDER_H
#define ATTRACKSEEDER_H

class AtTrack;

namespace AtPATTERN {

class AtTrackSeeder {
public:
   void SetTrackInitialParameters(AtTrack &track, double radiusFitFraction, int minHitsRadius, int maxHitsRadius);
};

} // namespace AtPATTERN

#endif
