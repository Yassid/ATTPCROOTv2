/// @file ex_eval.C
/// @brief Evaluate the 16C elastic-peak resolution for one UKF file (for param scans).
///
/// Same clean selection as ex_a1975.C (proton gate + IC gate + chi2 cut) on a single
/// <ukfFile> + <fribFile>, computes Ex, fits the elastic peak, and prints
/// "label nProton FWHM(MeV) peak(MeV)". Used by the UKF parameter scan.
///
///   root -b -q 'ex_eval.C("..._ukf_scan.root","..._FRIB.root","ms=1.0")'

#include <map>
#include <tuple>

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double kineEx(double m1, double m2, double m3, double m4, double Kp, double thlab, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4e = std::sqrt((std::cos(thlab) * omega2(s, m1 * m1, m2 * m2) * omega2(u, m2 * m2, m3 * m3) -
                           (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                            (2 * m2 * m2) +
                          s + u - m2 * m2);
   return m4e - m4;
}

void ex_eval(TString ukfFile, TString fribFile, TString label, TString gateFile = "proton_band.json",
             double Ebeam = 192.0, double icMin = 950, double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350,
             double chi2Cut = 5.0, double bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401, m_C16 = 16.0147 * u, m_p = 1.007825 * u;
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH1F *hex = new TH1F("hex", "", 240, -5, 10);
   TFile *fu = TFile::Open(ukfFile);
   TTree *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr;
   tu->SetBranchAddress("AtTrackingEvent", &te);
   TFile *fc = TFile::Open(fribFile);
   TTree *tc = (TTree *)fc->Get("cbmsim");
   TClonesArray *re = nullptr;
   tc->SetBranchAddress("AtRawEvent", &re);

   Long64_t N = std::min(tu->GetEntries(), tc->GetEntries());
   long nP = 0;
   for (Long64_t i = 0; i < N; ++i) {
      tc->GetEntry(i);
      double ic = -1;
      if (re->GetEntries() > 0) {
         auto *raw = (AtRawEvent *)re->At(0);
         if (raw && !raw->GetGenTraces().empty()) {
            auto &adc = raw->GetGenTraces()[0]->GetADC();
            double mx = -1e9;
            for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
               mx = std::max(mx, adc[b]);
            ic = mx;
         }
      }
      if (ic < icMin || ic > icMax)
         continue;
      tu->GetEntry(i);
      if (te->GetEntries() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev)
         continue;
      std::vector<AtTrack> orig = ev->GetTrackArray();
      std::map<int, AtTrack *> byID;
      for (auto &tr : orig)
         byID[tr.GetTrackID()] = &tr;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         auto &k = ft->GetKinematics();
         double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
         double c2n = ndf > 0 ? chi2 / ndf : 1e9, ke = k.kineticEnergy, thRad = k.theta;
         if (ke <= 0 || ke > 1000 || c2n > chi2Cut)
            continue;
         auto it = byID.find(ft->GetTrackID());
         if (it == byID.end())
            continue;
         auto r = spy.Estimate(*it->second);
         if (!r.valid || !pid.IsInside(r.sqrtdEdx, r.brho))
            continue;
         double ex = kineEx(m_C16, m_p, m_p, m_C16, Ebeam, thRad, ke);
         if (!std::isnan(ex)) {
            hex->Fill(ex);
            ++nP;
         }
      }
   }
   fu->Close();
   fc->Close();

   TF1 g("g", "gaus", -1.6, 1.6);
   g.SetParameters(hex->GetMaximum(), 0, 0.7);
   hex->Fit(&g, "QRN");
   double fwhm = 2.3548 * g.GetParameter(2);
   printf("SCAN %s  nProton=%ld  FWHM=%.3f MeV  peak=%.3f MeV\n", label.Data(), nP, fwhm, g.GetParameter(1));
}
