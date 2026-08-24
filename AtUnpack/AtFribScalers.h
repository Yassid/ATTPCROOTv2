#ifndef ATFRIBSCALERS_H
#define ATFRIBSCALERS_H

#include <Rtypes.h> // for ClassDef

#include <string>

class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @brief FRIBDAQ scaler counters for one run, or accumulated over several.
 *
 * Scalers are the DAQ's own counters. They are the only place in the data that records how much
 * BEAM went through the ion chamber, which makes them a route to the luminosity that owes nothing
 * to an elastic yield, an optical model or an acceptance.
 *
 * The unpackers in this directory never read them: AtFRIBHDFUnpacker walks `frib/evt` and stops,
 * and AtMergerHDFUnpacker documents `/scalers` as "not used here". Spyral does read them
 * (spyral/trace/frib_scalers.py); this class is that reader brought into the framework, and it is
 * validated bit-identical against Spyral's output on a1975 run 0016, all eleven counters.
 *
 * @warning THE COUNTERS ARE INCREMENTAL. Each scaler event holds the counts since the previous
 * one, so a run total is the SUM over every scaler event -- not the last value, which is what the
 * layout invites you to take. A run has a few thousand of them.
 *
 * @note Only the first eleven of the 32 stored words are named by FRIBDAQ; the rest are unassigned
 * and are not read.
 */
class AtFribScalers {
public:
   /// Index of each counter within a scaler dataset. Fixed by FRIBDAQ, not by us.
   enum EScaler {
      kClockFree = 0,  ///< time elapsed while the DAQ was running
      kClockLive = 1,  ///< time for which the DAQ could accept triggers
      kTriggerFree = 2,///< triggers received
      kTriggerLive = 3,///< triggers that actually produced events
      kIcSca = 4,      ///< ION CHAMBER counts -- the beam counter
      kMeshSca = 5,    ///< mesh signals
      kSi1Cfd = 6,     ///< Si detector 1
      kSi2Cfd = 7,     ///< Si detector 2
      kSiPM = 8,       ///< unclear; FRIBDAQ's own documentation says so
      kIcDs = 9,       ///< downscaled ion chamber
      kIcCfd = 10,     ///< unclear; equals kIcSca in every a1975 run checked
      kNScalers = 11
   };

   AtFribScalers() = default;

   /// @brief Sum every scaler event in one run file. Returns false if the file has no scaler group.
   /// Tries the legacy remerged layout (`frib/scaler`) first, then the raw merger layout
   /// (`scalers`), so it works on both without the caller having to know which it has.
   bool ReadRun(const std::string &fileName);

   /// @brief Accumulate another run into this one.
   void Add(const AtFribScalers &other);

   long long Get(EScaler which) const { return (which >= 0 && which < kNScalers) ? fCounts[which] : 0; }
   long GetScalerEvents() const { return fScalerEvents; }
   int GetRuns() const { return fRuns; }

   /// @brief Live fraction from the clock counters. Differs from the trigger version by ~2%, and
   /// which one belongs in a luminosity depends on what the yield was normalised to, so both exist.
   double GetLiveFractionClock() const;
   double GetLiveFractionTrigger() const;

   /// @brief Name of a counter, for printing.
   static const char *GetScalerName(EScaler which);

   void Print() const;

private:
   // A plain array, not std::array: the dictionary for std::array<long long,11> is not generated
   // here and loading the library segfaults before any user code runs.
   long long fCounts[kNScalers] = {};
   long fScalerEvents{0};
   int fRuns{0};

   ClassDef(AtFribScalers, 1);
};

#endif // ATFRIBSCALERS_H
