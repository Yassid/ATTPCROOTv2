#ifndef AtPSA_H
#define AtPSA_H

#include "AtHit.h" // IWYU pragma: keep

#include <Rtypes.h>

#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <vector>

class TClonesArray;
class AtRawEvent;
class AtEvent;
class AtPad;
class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @brief Abstract base class for processing AtPads (traces) into AtHits.
 *
 * @defgroup PSA
 */
class AtPSA {
private:
   // Access in PSA methods through getThreshold()
   Int_t fThreshold{-1};    ///< threshold of ADC value
   Int_t fThresholdlow{-1}; ///< threshold for Central pads

protected:
   TClonesArray *fMCSimPointArray{};

   Bool_t fUsingLowThreshold{false};

   // Variables from parameter file
   Double_t fBField{};
   Double_t fEField{};
   Int_t fTB0{};

   // This was hard coded too many places to leave as a variable...
   Int_t fNumTbs{512};        //< the number of time buckets used in taking data
   Int_t fTBTime{};           //< time duration of a time bucket in ns
   Int_t fEntTB{};            //< Timebucket of the entrance window
   Double_t fDriftVelocity{}; //< drift velocity of electron in cm/us
   Double_t fZk{};            //< Relative position of micromegas-cathode

   // Optional Spyral-style two-point z calibration: z = (window-tb)/(window-mm)*length [mm].
   // fWindowTB<=0 (default) => use the par-file geometric calibration (CalculateZGeo).
   Double_t fWindowTB{0};       //< time bucket of the window plane (z=0 reference)
   Double_t fMicromegasTB{10};  //< time bucket of the micromegas/pad plane (z=length reference)
   Double_t fDetLength{1000};   //< detector active length in mm
   // Optional per-pad time-bucket offset (electronics-timing calibration; padNum -> offset in TB).
   std::map<Int_t, Double_t> fPadTimeOffset; //!< transient: loaded per run, not persisted
   /// Beam enters through the pad plane (AtDigiPar::ReverseDrift). Read from the par, never set
   /// by hand: digitisation and PSA must agree, and a mismatch mirrors z instead of erroring.
   ///
   /// APPENDED at the end of the members on purpose -- AtPSA has a streamer (`+` in
   /// AtReconstructionLinkDef.h) and inserting mid-class shifts every later member's offset for
   /// any translation unit still compiled against the old header.
   Bool_t fReverseDrift{kFALSE};

   using HitVector = std::vector<std::unique_ptr<AtHit>>;

public:
   AtPSA() = default;
   virtual ~AtPSA() = default;

   virtual void Init();

   void SetThreshold(Int_t threshold);
   void SetThresholdLow(Int_t thresholdlow);
   int GetThreshold() { return fThreshold; }
   int GetThresholdLow() { return fThresholdlow; }

   void SetSimulatedEvent(TClonesArray *MCSimPointArray);

   /// Enable Spyral-style two-point z calibration: z = (windowTB - peakTB)/(windowTB - mmTB) * length [mm].
   /// (a1975 D2: windowTB=560, mmTB=10, length=1000). Pass windowTB<=0 to revert to CalculateZGeo.
   void SetSpyralZ(Double_t windowTB, Double_t mmTB = 10, Double_t length = 1000)
   {
      fWindowTB = windowTB; fMicromegasTB = mmTB; fDetLength = length;
   }
   /// Load a per-pad time-bucket offset map from a CSV (one value per line, line index = pad number;
   /// an optional non-numeric header line is skipped). Same format as Spyral's pad_time_correction.csv.
   void LoadPadTimeOffsets(const char *csvFile);
   void SetPadTimeOffset(Int_t padNum, Double_t tb) { fPadTimeOffset[padNum] = tb; }
   void ClearPadTimeOffsets() { fPadTimeOffset.clear(); }

   AtEvent Analyze(AtRawEvent &rawEvent);
   virtual void Analyze(AtRawEvent *rawEvent, AtEvent *event);
   virtual HitVector AnalyzePad(AtPad *pad) = 0;

   // virtual HitVector AnalyzeTrace(const std::vector<double> &trace) = 0;
   virtual std::unique_ptr<AtPSA> Clone() = 0;

protected:
   // Protected functions
   void TrackMCPoints(std::multimap<Int_t, std::size_t> &map,
                      AtHit &hit); //< Assign MC Points kinematics to each hit.

   [[deprecated]] Double_t CalculateZ(Double_t peakIdx); ///< Calculate z position in mm using the peak index.

   Double_t CalculateZGeo(Double_t peakIdx);
   /// Unified z from a peak time bucket: applies the per-pad time offset (if loaded for padNum),
   /// then the Spyral two-point calibration if enabled, else the geometric CalculateZGeo.
   /// PSA implementations should call this (with the pad number) instead of CalculateZGeo directly.
   Double_t CalibrateZ(Double_t peakIdx, Int_t padNum = -1);
   Double_t getThreshold(int padSize = -1);

   virtual double getZhitVariance(double zLoc, double zLocVar) const;
   virtual std::pair<double, double> getXYhitVariance() const;
   ClassDef(AtPSA, 6)
};

#endif
