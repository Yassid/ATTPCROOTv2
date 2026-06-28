#include "AtPSAMultiFit.h"

#include "AtHit.h"
#include "AtPad.h"

#include <FairLogger.h>

#include <Math/Point3D.h>    // for PositionVector3D
#include <Math/Point3Dfwd.h> // for XYZPoint
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>

#include <algorithm> // for max_element, sort, nth_element
#include <cmath>     // for exp, sin
#include <memory>    // for unique_ptr, make_unique
#include <numeric>   // for accumulate
#include <vector>

using XYZPoint = ROOT::Math::XYZPoint;

double AtPSAMultiFit::Response(double t, double t0, double tauTB)
{
   double u = (t - t0) / tauTB;
   if (u <= 0)
      return 0;
   return std::exp(-3.0 * u) * std::sin(u) * u * u * u; // GET nominal response (AtNominalResponse)
}

void AtPSAMultiFit::Init()
{
   AtPSA::Init();
   double du = 0.001, best = 0, bestU = 0, integ = 0;
   for (double u = du; u < M_PI; u += du) {
      double r = std::exp(-3.0 * u) * std::sin(u) * u * u * u;
      integ += r * du;
      if (r > best) {
         best = r;
         bestU = u;
      }
   }
   fRespPeakRedT = bestU;
   fRespPeakVal = best;
   fRespIntegral = integ;
   LOG(info) << "AtPSAMultiFit: tau=" << fPeakingTime << " us, response peak u=" << fRespPeakRedT << " (val "
             << fRespPeakVal << "), integral " << fRespIntegral << ", seedSigma=" << fSeedSigma;
}

double AtPSAMultiFit::EstimateNoise(const AtPad::trace &adc)
{
   // robust sigma = 1.4826 * MAD(adc) -- the pulse is a minority of samples so it barely moves the median.
   std::vector<double> v(adc.begin(), adc.end());
   size_t m = v.size() / 2;
   std::nth_element(v.begin(), v.begin() + m, v.end());
   double med = v[m];
   std::vector<double> dev(v.size());
   for (size_t i = 0; i < v.size(); ++i)
      dev[i] = std::abs(adc[i] - med);
   std::nth_element(dev.begin(), dev.begin() + m, dev.end());
   double sigma = 1.4826 * dev[m];
   return sigma < 1.0 ? 1.0 : sigma; // floor to avoid div-by-tiny
}

std::vector<int> AtPSAMultiFit::FindSeeds(const AtPad::trace &adc, double threshold, double noise, double tauTB) const
{
   const int N = (int)adc.size();
   // light 3-point smoothing to suppress single-bucket noise spikes
   std::vector<double> s(N);
   for (int i = 0; i < N; ++i) {
      double a = (i > 0) ? adc[i - 1] : adc[i];
      double b = adc[i];
      double c = (i < N - 1) ? adc[i + 1] : adc[i];
      s[i] = (a + b + c) / 3.0;
   }

   const bool absMode = fAbsProminence > 0;             // Spyral-style: absolute prominence on ALL peaks
   double promCut = absMode ? fAbsProminence : fSeedSigma * noise;
   int W = std::max(15, (int)(3 * tauTB));              // prominence search window

   // all local maxima above threshold (smoothed), strongest first
   std::vector<int> cand;
   for (int i = 21; i < 500; ++i)
      if (s[i] > threshold && s[i] >= s[i - 1] && s[i] > s[i + 1])
         cand.push_back(i);
   std::sort(cand.begin(), cand.end(), [&](int a, int b) { return s[a] > s[b]; });

   auto prominence = [&](int i) {
      double leftMin = s[i], rightMin = s[i];
      for (int j = i - 1; j >= std::max(1, i - W); --j) { if (s[j] > s[i]) break; leftMin = std::min(leftMin, s[j]); }
      for (int j = i + 1; j <= std::min(N - 2, i + W); ++j) { if (s[j] > s[i]) break; rightMin = std::min(rightMin, s[j]); }
      return s[i] - std::max(leftMin, rightMin);
   };

   std::vector<int> seeds;
   for (int c : cand) {
      bool spaced = true;
      for (int sd : seeds)
         if (std::abs(sd - c) < fMinSep) { spaced = false; break; }
      if (!spaced)
         continue;
      // Require prominence so broad baseline humps (low prominence) are rejected. With fPromPrimary the
      // PRIMARY peak needs it too (default); set false to restore "primary = threshold only" (matches AtPSAMax).
      bool needProm = fPromPrimary || absMode || !seeds.empty();
      if (!needProm || prominence(c) >= promCut)
         seeds.push_back(c);
      if ((int)seeds.size() >= fMaxPeaks)
         break;
   }
   std::sort(seeds.begin(), seeds.end());
   return seeds;
}

AtPSAMultiFit::HitVector AtPSAMultiFit::AnalyzePad(AtPad *pad)
{
   XYZPoint pos(pad->GetPadCoord().X(), pad->GetPadCoord().Y(), 0);
   if (pos.X() < -9000 || pos.Y() < -9000) {
      LOG(debug) << "Skipping pad, position is invalid";
      return {};
   }
   if (!(pad->IsPedestalSubtracted()))
      LOG(error) << "Pedestal should be subtracted to use this class!";

   const auto &adc = pad->GetADC();
   double threshold = getThreshold(pad->GetSizeID());
   double noise = EstimateNoise(adc);
   const double tauTB = fPeakingTime * 1000.0 / fTBTime;
   const double peakOff = fRespPeakRedT * tauTB;

   auto seeds = FindSeeds(adc, threshold, noise, tauTB);
   if (seeds.empty())
      return {};
   const int n = seeds.size();

   int lo = std::max(1, seeds.front() - 5);
   int hi = std::min(510, (int)(seeds.back() + 5 * tauTB));

   // sum-of-pulses model, params [A_i, t0_i]
   auto model = [this, n, tauTB](double *x, double *p) {
      double y = 0;
      for (int i = 0; i < n; ++i)
         y += p[2 * i] * Response(x[0], p[2 * i + 1], tauTB);
      return y;
   };
   TF1 tf("AtMFfit", model, lo, hi, 2 * n, 1, TF1::EAddToList::kNo);
   for (int i = 0; i < n; ++i) {
      double A0 = adc[seeds[i]] / fRespPeakVal;
      double t00 = seeds[i] - peakOff;
      tf.SetParameter(2 * i, A0);
      tf.SetParLimits(2 * i, 0, 5 * A0 + 10);
      tf.SetParameter(2 * i + 1, t00);
      tf.SetParLimits(2 * i + 1, t00 - tauTB, t00 + tauTB);
   }

   // fit a LOCAL, directory-detached histogram (no thread_local -> no teardown segfault)
   TH1D h("AtMFhist", "", fNumTbs, 0, fNumTbs);
   h.SetDirectory(nullptr);
   for (int i = 0; i < fNumTbs; ++i) {
      h.SetBinContent(i + 1, adc[i]);
      h.SetBinError(i + 1, noise);
   }
   auto r = h.Fit(&tf, "SQNR");
   bool good = r.Get() && r->IsValid();

   double ampCut = threshold; // keep gate matches AtPSAMax; prominence already vetted extra peaks
   HitVector hits;
   for (int i = 0; i < n; ++i) {
      double A = good ? tf.GetParameter(2 * i) : adc[seeds[i]] / fRespPeakVal;
      double t0 = good ? tf.GetParameter(2 * i + 1) : seeds[i] - peakOff;
      double peakTB = t0 + peakOff;
      double amp = A * fRespPeakVal;
      double Q = A * fRespIntegral * tauTB;

      if (amp < ampCut || peakTB < 20 || peakTB > 500) // amplitude + window gate
         continue;
      // FIT-QUALITY gate: reject pulses the fit could not constrain (A consistent with noise).
      double errA = good ? tf.GetParError(2 * i) : 0;
      if (good && errA > 0 && A / errA < fMinSignif) {
         LOG(debug) << "Reject pulse: A/sigma(A) = " << A / errA << " < " << fMinSignif;
         continue;
      }

      pos.SetZ(CalibrateZ(peakTB, pad->GetPadNum())); // base: per-pad offset + Spyral/Geo z
      auto hit = std::make_unique<AtHit>(pad->GetPadNum(), pos, amp);
      hit->SetTimeStamp(peakTB);
      hit->SetTimeStampCorr(peakTB);
      hit->SetTraceIntegral(Q);
      hits.push_back(std::move(hit));
   }
   return hits;
}

ClassImp(AtPSAMultiFit);
