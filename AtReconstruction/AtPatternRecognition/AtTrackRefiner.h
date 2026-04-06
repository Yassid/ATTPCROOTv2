#ifndef ATTRACKREFINER_H
#define ATTRACKREFINER_H

#include <vector>

class AtTrack;

namespace AtTools {
class AtTrackTransformer;
}

namespace AtPATTERN {

class AtTrackRefiner {
public:
   void OrderClustersAlongTrack(AtTrack &track);
   void SelectAndMergeTracks(std::vector<AtTrack> &tracks, AtTools::AtTrackTransformer &trackTransformer,
                             double clusterRadius, double clusterDistance, double vertexRadiusXY = 80.0,
                             double mergeDist = 30.0, double minLabTheta = 10.0);
};

} // namespace AtPATTERN

#endif
