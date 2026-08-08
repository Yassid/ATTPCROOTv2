#ifndef ATPSAMULTIFIT_H
#define ATPSAMULTIFIT_H

#include "AtPSA.h" // for AtPSA, AtPSA::HitVector

#include "AtPad.h" // for AtPad::trace

#include <Rtypes.h> // for Double_t, Int_t, ClassDefOverride

#include <memory> // for unique_ptr
#include <vector>

class AtPad;
class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @brief Multi-peak PSA: search peaks, fit a SUM of GET response pulses to the raw trace.
 *
 * Finds candidate peaks (local maxima above threshold) in a pad trace, then fits a sum of the
 * electronics response function used in digitization
 * (ElectronicResponse::AtNominalResponse: exp(-3u) sin(u) u^3, u = (t-t0)/tau) -- one component
 * per peak -- directly to the trace. The fit simultaneously validates each pulse against the
 * true response and DECONVOLVES overlapping pulses (two tracks crossing the same pad at
 * different drift times), returning one AtHit per fitted component. This resolves the
 * single-peak limitation of AtPSAMax that loses hits on shared pads.
 *
 * @ingroup PSA
 */
class AtPSAMultiFit : public AtPSA {
private:
   Double_t fPeakingTime{0.720}; ///< electronics peaking time tau [us] (GET nominal ~0.72)
   Int_t fMaxPeaks{4};           ///< max simultaneous pulses to fit per pad
   Int_t fMinSep{8};             ///< min separation between seed peaks [TB]
   Double_t fSeedSigma{4.0};     ///< an ADDITIONAL seed peak must have prominence > fSeedSigma*noiseRMS
   Double_t fAbsProminence{0};   ///< if >0, use this ABSOLUTE prominence (ADC) for ALL peaks (Spyral-style)
                                 ///< instead of the noise-adaptive secondary-only fSeedSigma cut
   Double_t fMinSignif{3.0};     ///< FIT-QUALITY gate: keep a fitted pulse only if amplitude/error(amplitude)
                                 ///< >= fMinSignif (a real pulse has A >> sigma(A); noise has A ~ sigma(A))
   Bool_t fPromPrimary{false};   ///< prominence on PRIMARY peak. Default false: primary=threshold only (matches
                                 ///< AtPSAMax). true rejects broad baseline humps BUT also kills diffusion-broadened
                                 ///< far-drift pulses -> collapses track z. Only enable if noise pads are a problem.
                                 ///< (FWHM ~50-90 TB, low prominence) that a threshold-only primary would keep.
   Bool_t fFloatTau{false};      ///< let the peaking time float in the fit (else fixed)
   Bool_t fReportT0{false};      ///< report the fitted impulse arrival t0 instead of the response peak.
                                 ///< Default false keeps the historical peak-referenced time, which lags
                                 ///< the impulse by fRespPeakRedT*tau (~5 TB at tau=0.72 us, 6 MHz) and is
                                 ///< what every existing z calibration was tuned against. Turning it on
                                 ///< shifts every hit z by that lag -- correct, but recalibrate first.
   Double_t fFitChi2Cut{0};      ///< if >0, reject a peak whose LOCAL reduced chi2 (data vs fitted GET pulse)
                                 ///< exceeds this -> trace is not a real pulse SHAPE (noise hump). >0 also turns
                                 ///< on float-tau so a real BROAD pulse fits and survives. chi2 stored in
                                 ///< AtHit::ChargeVariance (computed whenever float-tau is on, for diagnostics).
   Double_t fFitRelErr{0.1};     ///< fractional model error: effective sigma^2 = noise^2 + (relErr*model)^2.
                                 ///< Keeps the chi2 amplitude-relative so a BRIGHT well-fit pulse (e.g. beam)
                                 ///< isn't penalized for the GET model's ~relErr shape approximation.
   // Spyral-style z calibration (SetSpyralZ) + per-pad time map are now provided by the base AtPSA.

   // response constants, precomputed in Init() (reduced-time units)
   double fRespPeakRedT{0};  ///< reduced time of the response maximum
   double fRespPeakVal{1};   ///< response value at its maximum (for amplitude init)
   double fRespIntegral{1};  ///< integral of the response over reduced time (for charge)

public:
   void Init() override;
   HitVector AnalyzePad(AtPad *pad) override;
   std::unique_ptr<AtPSA> Clone() override { return std::make_unique<AtPSAMultiFit>(*this); }

   void SetPeakingTime(Double_t tauUs) { fPeakingTime = tauUs; }
   void SetMaxPeaks(Int_t n) { fMaxPeaks = n; }
   void SetMinSeparation(Int_t tb) { fMinSep = tb; }
   void SetSeedSigma(Double_t s) { fSeedSigma = s; }
   void SetProminence(Double_t p) { fAbsProminence = p; } ///< absolute prominence (ADC), Spyral-style; 0=adaptive
   void SetMinSignificance(Double_t s) { fMinSignif = s; }
   void SetProminenceOnPrimary(Bool_t v) { fPromPrimary = v; } ///< false = old behavior (primary=threshold only)
   void SetFloatPeakingTime(Bool_t v) { fFloatTau = v; }
   void SetReportT0(Bool_t v) { fReportT0 = v; } ///< report impulse arrival instead of the response peak
   void SetFitChi2Cut(Double_t c) { fFitChi2Cut = c; } ///< local reduced-chi2 shape gate; 0 disables
   void SetFitRelErr(Double_t e) { fFitRelErr = e; }   ///< fractional model error in the chi2 effective sigma
   // SetSpyralZ / LoadPadTimeOffsets are inherited from AtPSA.

private:
   /// GET response to a unit charge arriving at t0 (all in TB), tau in TB. 0 for t<t0.
   static double Response(double t, double t0, double tauTB);
   /// Robust per-pad noise RMS from the trace (median + MAD; insensitive to the pulse).
   static double EstimateNoise(const AtPad::trace &adc);
   /// Prominent local maxima: smoothed local max, prominence > fSeedSigma*noise above local valleys,
   /// >= fMinSep apart, strongest first (up to fMaxPeaks). Rejects noise-wiggle seeds.
   std::vector<int> FindSeeds(const AtPad::trace &adc, double threshold, double noise, double tauTB) const;

   ClassDefOverride(AtPSAMultiFit, 1)
};

#endif
