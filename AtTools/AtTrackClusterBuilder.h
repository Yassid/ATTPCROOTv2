#ifndef ATTRACKCLUSTERBUILDER_H
#define ATTRACKCLUSTERBUILDER_H

#include "AtHit.h"
#include "AtTrackTransformer.h"

#include <TMatrixDSymfwd.h>
#include <TMatrixTSym.h>

#include <memory>
#include <vector>

class AtHitCluster;

namespace AtTools {

struct AtTrackClusterBuilderConfig {
   double coefT{0.00009};
   double coefL{0.0000009};
   double driftVel{1.0};
   double samplingRate{0.320};
   double padResXY{2.3};
   double padResZ{3.45};
   AtTrackTransformer::CovarianceMode covarianceMode{AtTrackTransformer::CovarianceMode::TransformerDirect};
};

class AtTrackClusterBuilder {
public:
   explicit AtTrackClusterBuilder(AtTrackClusterBuilderConfig config);

   std::shared_ptr<AtHitCluster> BuildCluster(const std::vector<AtHit> &hits, int clusterID) const;

private:
   struct TransformerDirectClusterStats {
      double x{0};
      double y{0};
      double z{0};
      double sigmaX2{0};
      double sigmaY2{0};
      double sigmaZ2{0};
      double totalCharge{0};
      int timeStamp{0};
      int nHits{0};
      bool valid{false};
   };

   AtHit::XYZVector GetPerHitVariance(const AtHit &hit) const;
   TransformerDirectClusterStats BuildTransformerDirectClusterStats(const std::vector<AtHit> &hits) const;
   TMatrixDSym BuildTransformerDirectCovariance(const TransformerDirectClusterStats &stats) const;
   AtHitCluster BuildHitClusterOnlineCluster(const std::vector<AtHit> &hits) const;
   TMatrixDSym BuildDiagonalOnlyCovariance(const TMatrixDSym &src) const;

   AtTrackClusterBuilderConfig fConfig;
};

} // namespace AtTools

#endif
