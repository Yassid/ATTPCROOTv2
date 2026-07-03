#ifndef ATDIRDEDXCLEANER_H
#define ATDIRDEDXCLEANER_H

#include "AtDataCleaner.h"

namespace AtTools {
namespace DataCleaning {

/**
 * @brief Direction + dE/dx continuity noise cleaner.
 *
 * A hit is kept (signal) if it has at least fMinDeg neighbours that are both DIRECTION-
 * continuous (the local tangents and the connecting segment are collinear) AND dE/dx-
 * continuous (smoothed local charge ratio above fQRatio). Isolated / off-track background
 * gets degree ~0 and is removed. On AT-TPC sim truth this keeps ~95% of signal (incl. the
 * sparse heavy recoil) while removing ~98.5% of background; it is a per-hit filter with no
 * clustering commitment, so it can front any downstream track finder.
 */
class AtDirDeDxCleaner : public AtDataCleaner {
protected:
   int fK;         //<! neighbours for tangent + candidate edges
   double fCosSeg; //<! min |segment . tangent| (both endpoints)
   double fCosTan; //<! min |tangent_i . tangent_j|
   double fRMax;   //<! max edge length [mm]
   double fQRatio; //<! min smoothed-charge ratio (dE/dx continuity)
   int fMinDeg;    //<! min continuous-neighbour count to keep a hit
   bool fSmoothQ;  //<! use median neighbourhood charge (rides Bragg profile)
   double fQKeep;  //<! charge-magnitude override: always keep a hit with charge > this (protects
                   //<! high-dE/dx signal, e.g. the 17C recoil, at ~zero noise cost since noise is
                   //<! baseline-low). 17C kept 65%->84% at fQKeep=150 with no change in noise removal.

public:
   AtDirDeDxCleaner(int k = 12, double cosSeg = 0.78, double cosTan = 0.72, double rMax = 32.0, double qRatio = 0.65,
                    int minDeg = 1, bool smoothQ = true, double qKeep = 150.0)
      : fK(k), fCosSeg(cosSeg), fCosTan(cosTan), fRMax(rMax), fQRatio(qRatio), fMinDeg(minDeg), fSmoothQ(smoothQ),
        fQKeep(qKeep)
   {
   }
   HitCloud CleanData(const HitCloud &hits) override;
};

} // namespace DataCleaning
} // namespace AtTools

#endif // ATDIRDEDXCLEANER_H
