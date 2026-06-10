/// @file ex_eval_pd.C
/// @brief Evaluate the 15C ground-state peak resolution for a deuteron-fit param set.
///
/// Deuteron analog of ex_eval.C: same clean selection as ex_pd_a1975.C (deuteron
/// gate + IC gate + chi2) over a set of runs, builds the 15C Ex spectrum, fits the
/// g.s. peak, and prints "SCAN <label> nD FWHM peak". Used by the (p,d) UKF
/// parameter scan. Reads <run>_ukf<fitSuffix>.root from fitDir and <run>_FRIB.root
/// from inDir.
///
///   root -b -q 'pd/ex_eval_pd.C("run_0106,run_0107","ms=0.75","_d_ms075")'

#include <map>
static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double kineEx(double m1, double m2, double m3, double m4, double Kp, double thlab, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4e = std::sqrt((std::cos(thlab) * omega2(s, m1 * m1, m2 * m2) * omega2(uu, m2 * m2, m3 * m3) -
                           (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                            (2 * m2 * m2) +
                          s + uu - m2 * m2);
   return m4e - m4;
}

void ex_eval_pd(TString runsCSV, TString label, TString fitSuffix = "_d",
                TString inDir = "/mnt/f/a1975/reco/", TString fitDir = "/mnt/f/a1975/reco_pd/",
                TString gateFile = "pid/deuteron_band.json", double Ebeam = 192.0, double icMin = 950,
                double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350, double chi2Cut = 5.0, double bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401, m_C16 = 16.0147 * u, m_p = 1.007825 * u, m_d = 2.01410178 * u, m_C15 = 15.0105993 * u;
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH1F *hex = new TH1F("hex", "", 240, -10, 20);
   TObjArray *runs = runsCSV.Tokenize(",");
   long nD = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString uf = fitDir + run + "_ukf" + fitSuffix + ".root", ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(uf) || gSystem->AccessPathName(ff))
         continue;
      TFile *fu = TFile::Open(uf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      TFile *fc = TFile::Open(ff);
      TTree *tc = (TTree *)fc->Get("cbmsim");
      TClonesArray *re = nullptr;
      tc->SetBranchAddress("AtRawEvent", &re);
      Long64_t N = std::min(tu->GetEntries(), tc->GetEntries());
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
            double ex = kineEx(m_C16, m_p, m_d, m_C15, Ebeam, thRad, ke);
            if (!std::isnan(ex)) {
               hex->Fill(ex);
               ++nD;
            }
         }
      }
      fu->Close();
      fc->Close();
   }

   // fit the g.s. peak: seed at the max bin within [-2,2.5], fit gaus in a +/-1.3 window
   hex->GetXaxis()->SetRangeUser(-2, 2.5);
   double pk = hex->GetBinCenter(hex->GetMaximumBin()), amp = hex->GetMaximum();
   hex->GetXaxis()->SetRange(0, 0);
   TF1 g("g", "gaus", pk - 1.3, pk + 1.3);
   g.SetParameters(amp, pk, 0.5);
   hex->Fit(&g, "QRN");
   double fwhm = 2.3548 * g.GetParameter(2);
   printf("SCAN %s  nD=%ld  FWHM=%.3f MeV  peak=%.3f MeV\n", label.Data(), nD, fwhm, g.GetParameter(1));
}
