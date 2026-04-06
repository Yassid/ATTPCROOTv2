#include "AtTrackClusterBuilder.h"

#include "AtHit.h"
#include "AtHitCluster.h"

#include <TMatrixTSym.h>

AtTools::AtTrackClusterBuilder::AtTrackClusterBuilder(AtTrackClusterBuilderConfig config) : fConfig(config) {}

AtHit::XYZVector AtTools::AtTrackClusterBuilder::GetPerHitVariance(const AtHit &hit) const
{
   auto pos = hit.GetPosition();
   double driftTime = pos.Z() / (10.0 * fConfig.driftVel);
   double varT = 100.0 * fConfig.coefT * 2.0 * driftTime;
   double varL = 100.0 * fConfig.coefL * 2.0 * driftTime;
   double tbResMM = fConfig.driftVel * fConfig.samplingRate * 10.0;
   double varTB = tbResMM * tbResMM / 12.0;
   return {fConfig.padResXY * fConfig.padResXY + varT, fConfig.padResXY * fConfig.padResXY + varT,
           fConfig.padResZ * fConfig.padResZ + varTB + varL};
}

AtTools::AtTrackClusterBuilder::TransformerDirectClusterStats
AtTools::AtTrackClusterBuilder::BuildTransformerDirectClusterStats(const std::vector<AtHit> &hits) const
{
   TransformerDirectClusterStats stats;
   double varX = 0;
   double varY = 0;
   double varZ = 0;

   for (const auto &hit : hits) {
      auto pos = hit.GetPosition();
      double q = hit.GetCharge();
      auto hitVar = GetPerHitVariance(hit);

      stats.x += pos.X() * q;
      stats.y += pos.Y() * q;
      stats.z += pos.Z();
      stats.totalCharge += q;
      stats.timeStamp += hit.GetTimeStamp();
      stats.nHits++;

      varX += q * q * hitVar.X();
      varY += q * q * hitVar.Y();
      varZ += q * q * hitVar.Z();
   }

   if (stats.nHits == 0 || stats.totalCharge <= 0)
      return stats;

   stats.x /= stats.totalCharge;
   stats.y /= stats.totalCharge;
   stats.z /= stats.nHits;
   stats.timeStamp /= stats.nHits;
   stats.sigmaX2 = varX / (stats.totalCharge * stats.totalCharge);
   stats.sigmaY2 = varY / (stats.totalCharge * stats.totalCharge);
   stats.sigmaZ2 = varZ / (stats.totalCharge * stats.totalCharge);
   stats.valid = true;
   return stats;
}

TMatrixDSym AtTools::AtTrackClusterBuilder::BuildTransformerDirectCovariance(const TransformerDirectClusterStats &stats) const
{
   TMatrixDSym cov(3);
   cov(0, 1) = 0;
   cov(1, 2) = 0;
   cov(2, 0) = 0;
   cov(0, 0) = stats.sigmaX2;
   cov(1, 1) = stats.sigmaY2;
   cov(2, 2) = stats.sigmaZ2;
   return cov;
}

AtHitCluster AtTools::AtTrackClusterBuilder::BuildHitClusterOnlineCluster(const std::vector<AtHit> &hits) const
{
   AtHitCluster cluster;
   for (const auto &hit : hits) {
      AtHit hitWithVar(hit);
      hitWithVar.SetPositionVariance(GetPerHitVariance(hitWithVar));
      cluster.AddHit(hitWithVar);
   }
   return cluster;
}

TMatrixDSym AtTools::AtTrackClusterBuilder::BuildDiagonalOnlyCovariance(const TMatrixDSym &src) const
{
   TMatrixDSym diag(src);
   for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
         if (row != col)
            diag(row, col) = 0.0;
      }
   }
   return diag;
}

std::shared_ptr<AtHitCluster> AtTools::AtTrackClusterBuilder::BuildCluster(const std::vector<AtHit> &hits,
                                                                           int clusterID) const
{
   auto stats = BuildTransformerDirectClusterStats(hits);
   if (!stats.valid)
      return nullptr;

   auto cluster = std::make_shared<AtHitCluster>();
   cluster->SetClusterID(clusterID);
   cluster->SetCharge(stats.totalCharge);
   cluster->SetPosition({stats.x, stats.y, stats.z});
   cluster->SetTimeStamp(stats.timeStamp);
   cluster->SetCovMatrix(BuildTransformerDirectCovariance(stats));

   if (fConfig.covarianceMode == AtTrackTransformer::CovarianceMode::TransformerDirect)
      return cluster;

   auto onlineCluster = BuildHitClusterOnlineCluster(hits);

   if (fConfig.covarianceMode == AtTrackTransformer::CovarianceMode::HitClusterOnline) {
      cluster->SetCovMatrix(onlineCluster.GetCovMatrix());
      return cluster;
   }

   if (fConfig.covarianceMode == AtTrackTransformer::CovarianceMode::HitClusterOnlineDiagOnly) {
      cluster->SetCovMatrix(BuildDiagonalOnlyCovariance(onlineCluster.GetCovMatrix()));
      return cluster;
   }

   if (fConfig.covarianceMode == AtTrackTransformer::CovarianceMode::HitClusterOnlineConsistent) {
      auto consistentCluster = std::make_shared<AtHitCluster>(onlineCluster);
      consistentCluster->SetClusterID(clusterID);
      consistentCluster->SetTimeStamp(stats.timeStamp);
      return consistentCluster;
   }

   return cluster;
}
