/// @file ex_overlay.C
/// @brief Overlay two 16C Ex spectra from two UKF files (e.g. different MeasurementSigma).
///   root -b -q 'ex_overlay.C("..._s05.root","..._s20.root","..._FRIB.root","measSigma=0.5","measSigma=2.0")'

#include <map>

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double kineEx(double m1, double m2, double m3, double m4, double Kp, double thlab, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4e = std::sqrt((std::cos(thlab) * omega2(s, m1 * m1, m2 * m2) * omega2(u, m2 * m2, m3 * m3) -
                           (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                            (2 * m2 * m2) +
                          s + u - m2 * m2);
   return m4e - m4;
}

static void fillEx(TString ukfFile, TString fribFile, TH1F *h, AtTools::AtParticleID &pid, AtTools::AtSpyralPID &spy,
                   double icMin, double icMax, int icTbLo, int icTbHi, double chi2Cut)
{
   const double u = 931.49401, m_C16 = 16.0147 * u, m_p = 1.007825 * u;
   TFile *fu = TFile::Open(ukfFile);
   TTree *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr;
   tu->SetBranchAddress("AtTrackingEvent", &te);
   TFile *fc = TFile::Open(fribFile);
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
         double ex = kineEx(m_C16, m_p, m_p, m_C16, 192.0, thRad, ke);
         if (!std::isnan(ex))
            h->Fill(ex);
      }
   }
   fu->Close();
   fc->Close();
}

void ex_overlay(TString ukfA, TString ukfB, TString fribFile, TString labelA, TString labelB,
                TString gateFile = "proton_band.json", double icMin = 950, double icMax = 1350, int icTbLo = 1000,
                int icTbHi = 1350, double chi2Cut = 5.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TH1F *hA = new TH1F("hA", "16C E_{x}: UKF MeasurementSigma comparison;E_{x} [MeV];protons (norm.)", 200, -4, 8);
   TH1F *hB = new TH1F("hB", "", 200, -4, 8);
   fillEx(ukfA, fribFile, hA, pid, spy, icMin, icMax, icTbLo, icTbHi, chi2Cut);
   fillEx(ukfB, fribFile, hB, pid, spy, icMin, icMax, icTbLo, icTbHi, chi2Cut);

   auto fwhm = [](TH1F *h) {
      TF1 g("g", "gaus", -1.6, 1.6);
      g.SetParameters(h->GetMaximum(), 0, 0.7);
      h->Fit(&g, "QRN");
      return 2.3548 * g.GetParameter(2);
   };
   double fA = fwhm(hA), fB = fwhm(hB);
   if (hA->Integral() > 0)
      hA->Scale(1.0 / hA->GetMaximum());
   if (hB->Integral() > 0)
      hB->Scale(1.0 / hB->GetMaximum());

   TCanvas *c = new TCanvas("c", "exov", 900, 680);
   hA->SetLineColor(kRed + 1);
   hA->SetLineWidth(2);
   hB->SetLineColor(kGray + 2);
   hB->SetLineWidth(2);
   hA->Draw("hist");
   hB->Draw("hist same");
   TLegend *leg = new TLegend(0.58, 0.72, 0.88, 0.88);
   leg->AddEntry(hA, Form("%s (FWHM %.3f)", labelA.Data(), fA), "l");
   leg->AddEntry(hB, Form("%s (FWHM %.3f)", labelB.Data(), fB), "l");
   leg->Draw();
   c->SaveAs("ex_overlay.png");
   printf("FWHM  %s=%.3f   %s=%.3f MeV\n", labelA.Data(), fA, labelB.Data(), fB);
}
