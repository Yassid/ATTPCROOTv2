/// @file scan_corr.C
/// @brief Scan beam energy against the explorer's theta correction, and score each combination on
///        how FLAT the reference peak sits in theta_lab.
///
/// The correction is the one from the a2091 explorer, theta -> theta - slope*(KE - pivot). The two
/// parameters are not independent: theta - slope*(KE-pivot) = theta - slope*KE + slope*pivot, so
/// the pivot only adds a constant angle offset. Ebeam and pivot therefore trade against each other
/// and a scan over all three is largely degenerate -- which is exactly why the score here is NOT
/// "is the peak at zero" (any of them can arrange that) but "is the peak at the SAME place at every
/// angle". Flatness is the part that only a genuinely correct correction can achieve.
///
/// Score = RMS of the per-theta-slice peak positions, over slices with enough statistics. The peak
/// is located by the maximum of a smoothed Ex histogram inside a search window, NOT by a gaussian
/// fit, because at the wrong parameters there is no peak to fit and a failed fit would silently
/// drop the very points that should score badly.
///
///   root -b -q 'pp/scan_corr.C("~/a1975_panels/dt_kin.root",187,197,1)'                  // (d,t)
///   root -b -q 'pp/scan_corr.C("cache.root",190,200,1,"",0,0.25,0.025,5,15,55,5,\
///                16.0147013,2.0135532,3.01550072,15.0105993)'
///
/// Reports the best (Ebeam, slope) and the flatness achieved, plus the no-correction baseline so
/// the gain is visible. A best point at slope=0 means the correction buys nothing.

static double omS_(double x, double y, double z) { return sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static double exS_(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1*m1+m2*m2+2*m2*E1, u = m2*m2+m3*m3-2*m2*E3;
   double a = (cos(th)*omS_(s,m1*m1,m2*m2)*omS_(u,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   return a < 0 ? std::nan("") : sqrt(a) - m4;
}

void scan_corr(TString cache, double eLo = 187, double eHi = 197, double dE = 1.0,
               TString outTag = "", double slopeLo = 0.0, double slopeHi = 0.25,
               double dSlope = 0.025, double pivot = 5.0,
               double thLo = 15, double thHi = 55, double dth = 5,
               double mBeamAmu = 16.0147013, double mTargAmu = 2.0135532,
               double mEjectAmu = 3.01550072, double mResidAmu = 15.0105993,
               double chi2max = 5.0, double searchLo = -6, double searchHi = 14, int minN = 150)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = mBeamAmu*u, m2 = mTargAmu*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;
   TString dir = gSystem->DirName(__FILE__);
   // expand a leading ~ only; assigning ExpandPathName's result back onto its own argument
   // clobbers the string before it is read
   if (cache.BeginsWith("~")) cache.Replace(0, 1, gSystem->Getenv("HOME"));

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk`\n"); return; }
   float ke, th, c2;
   t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&th); t->SetBranchAddress("chi2ndf",&c2);
   std::vector<float> vk, vt;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ke <= 0 || c2 > chi2max) continue;
      if (th < thLo || th >= thHi) continue;
      vk.push_back(ke); vt.push_back(th);
   }
   f->Close();
   printf("\n=== scan_corr: %s ===\n%zu tracks in theta[%.0f,%.0f], chi2/ndf<%.1f\n",
          gSystem->BaseName(cache), vk.size(), thLo, thHi, chi2max);
   printf("theta correction: theta -> theta - slope*(KE - %.1f)\n", pivot);
   printf("score = RMS across theta slices of the peak position (lower = flatter)\n\n");

   const int nSl = (int)((thHi - thLo)/dth);
   auto evaluate = [&](double E, double slope, double &rms, double &mean, int &nUsed) {
      std::vector<double> pk;
      for (int s = 0; s < nSl; ++s) {
         double a0 = thLo + s*dth, a1 = a0 + dth;
         TH1D h("h","",200,searchLo,searchHi); h.SetDirectory(nullptr);
         for (size_t i = 0; i < vk.size(); ++i) {
            double thc = vt[i] - slope*(vk[i] - pivot);
            if (thc < a0 || thc >= a1) continue;
            double e = exS_(m1,m2,m3,m4,E,thc*TMath::DegToRad(),vk[i]);
            if (!std::isnan(e)) h.Fill(e);
         }
         if (h.GetEntries() < minN) continue;
         h.Smooth(2);
         pk.push_back(h.GetBinCenter(h.GetMaximumBin()));
      }
      nUsed = (int)pk.size();
      if (nUsed < 3) { rms = 1e9; mean = 0; return; }
      double s1 = 0, s2 = 0;
      for (double v : pk) { s1 += v; s2 += v*v; }
      mean = s1/nUsed;
      rms = sqrt(std::max(0.0, s2/nUsed - mean*mean));
   };

   double bestRms = 1e9, bestE = 0, bestSlope = 0, bestMean = 0;
   printf("%-8s", "Ebeam");
   for (double sl = slopeLo; sl <= slopeHi + 1e-9; sl += dSlope) printf("  s=%.3f", sl);
   printf("\n");
   for (double E = eLo; E <= eHi + 1e-9; E += dE) {
      printf("%-8.1f", E);
      for (double sl = slopeLo; sl <= slopeHi + 1e-9; sl += dSlope) {
         double rms, mean; int nU;
         evaluate(E, sl, rms, mean, nU);
         printf("  %6s", rms > 1e8 ? "-" : Form("%.3f", rms));
         if (rms < bestRms) { bestRms = rms; bestE = E; bestSlope = sl; bestMean = mean; }
      }
      printf("\n");
   }

   double rms0, mean0; int n0;
   evaluate(0.5*(eLo+eHi), 0.0, rms0, mean0, n0);
   printf("\nbaseline (no correction, Ebeam %.1f): RMS %.3f MeV, peak mean %+.2f\n",
          0.5*(eLo+eHi), rms0, mean0);
   printf("BEST: Ebeam %.1f, slope %.3f deg/MeV -> RMS %.3f MeV, peak mean %+.2f\n",
          bestE, bestSlope, bestRms, bestMean);
   if (bestSlope == 0)
      printf("  NOTE best slope is 0 -- the theta correction buys nothing on this sample.\n");
   printf("  (pivot fixed at %.1f; it only shifts theta by a constant, so it is degenerate with\n"
          "   Ebeam -- rescan with a different pivot only to re-centre, not to improve flatness)\n", pivot);

   // flatness profile at the best point, so the residual shape is visible
   printf("\nper-slice peak position at the best point:\n%-12s %10s %10s\n", "theta_lab", "no corr", "best");
   for (int s = 0; s < nSl; ++s) {
      double a0 = thLo + s*dth, a1 = a0 + dth;
      auto peak = [&](double E, double slope) -> double {
         TH1D h("h","",200,searchLo,searchHi); h.SetDirectory(nullptr);
         for (size_t i = 0; i < vk.size(); ++i) {
            double thc = vt[i] - slope*(vk[i] - pivot);
            if (thc < a0 || thc >= a1) continue;
            double e = exS_(m1,m2,m3,m4,E,thc*TMath::DegToRad(),vk[i]);
            if (!std::isnan(e)) h.Fill(e);
         }
         if (h.GetEntries() < minN) return std::nan("");
         h.Smooth(2);
         return h.GetBinCenter(h.GetMaximumBin());
      };
      double p0 = peak(0.5*(eLo+eHi), 0.0), pb = peak(bestE, bestSlope);
      printf("%-12s %10s %10s\n", Form("%.0f-%.0f",a0,a1),
             std::isnan(p0)?"-":Form("%+.2f",p0), std::isnan(pb)?"-":Form("%+.2f",pb));
   }
}
