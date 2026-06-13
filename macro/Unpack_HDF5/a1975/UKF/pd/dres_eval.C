/// @file dres_eval.C
/// @brief Resolution metric for a (p,d) deuteron-fit param set — FAST (no IC read).
///
/// Reads <run>_ukf<suffix>.root deuteron fits for a subset of runs, applies the
/// deuteron PID gate + chi2 cut (NO FRIB/IC gate, for speed; resolution ranking is
/// unaffected), builds the 15C Ex spectrum, SELF-CALIBRATES (shifts the g.s. to 0),
/// and reports the resolution focused on the g.s./0.74 doublet:
///   - doublet peak/valley  (the g.s. vs 0.74 separation; higher = better)
///   - 3.10 MeV peak FWHM   (isolated single peak; robust resolution proxy)
///   - g.s. & 0.74 positions
/// Prints one "DRES <label> ..." line for the scan driver to collect.
///
///   root -b -q 'pd/dres_eval.C("run_0106,run_0107","baseline","_d")'

#include <map>
static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double exOf(double m1, double m2, double m3, double m4, double Kp, double th, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4e = std::sqrt((std::cos(th) * omega2(s, m1 * m1, m2 * m2) * omega2(uu, m2 * m2, m3 * m3) -
                           (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                            (2 * m2 * m2) +
                          s + uu - m2 * m2);
   return m4e - m4;
}

void dres_eval(TString runsCSV, TString label, TString fitSuffix = "_d", TString fitDir = "/mnt/f/a1975/reco_pd/",
               TString gateFile = "pid/deuteron_band.json", double Ebeam = 192.0, double chi2Cut = 5.0,
               double bField = 2.85, TString fitBase = "_ukf")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401, mC16 = 16.0147 * u, mp = 1.007825 * u, md = 2.01410178 * u, mC15 = 15.0105993 * u;
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH1F *h = new TH1F("h", "", 240, -6, 14);
   h->SetDirectory(nullptr);
   TObjArray *runs = runsCSV.Tokenize(",");
   long nD = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString uf = fitDir + run + fitBase + fitSuffix + ".root";
      if (gSystem->AccessPathName(uf))
         continue;
      TFile *fu = TFile::Open(uf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      if (!tu) {
         fu->Close();
         continue;
      }
      TClonesArray *te = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      for (Long64_t i = 0; i < tu->GetEntries(); ++i) {
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
            double c2n = ndf > 0 ? chi2 / ndf : 1e9, ke = k.kineticEnergy, thR = k.theta;
            if (ke <= 0 || ke > 1000 || c2n > chi2Cut)
               continue;
            auto it = byID.find(ft->GetTrackID());
            if (it == byID.end())
               continue;
            auto r = spy.Estimate(*it->second);
            if (!r.valid || !pid.IsInside(r.sqrtdEdx, r.brho))
               continue;
            double ex = exOf(mC16, mp, md, mC15, Ebeam, thR, ke);
            if (!std::isnan(ex))
               h->Fill(ex);
            ++nD;
         }
      }
      fu->Close();
   }

   h->Rebin(2); // 0.083 MeV bins
   // Doublet model: g.s.+0.74 (separation FIXED 0.74), COMMON sigma, linear bg.
   //   [0]=g.s.pos  [1]=sigma(=resolution)  [2]=A_gs  [3]=A_074  [4]=bg0  [5]=bg1
   TF1 dg("dg",
          "[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x", -1.5, 2.3);
   double mx = h->GetMaximum();
   dg.SetParameters(0.45, 0.32, 0.6 * mx, 0.6 * mx, 0.1 * mx, 0.0);
   dg.SetParLimits(0, -0.3, 1.0);
   dg.SetParLimits(1, 0.12, 0.55);
   dg.SetParLimits(2, 0, 2 * mx);
   dg.SetParLimits(3, 0, 2 * mx);
   int st = h->Fit(&dg, "QRNS");
   double gs = dg.GetParameter(0), sig = dg.GetParameter(1), fwhm = 2.3548 * sig;

   // fit-free doublet peak/valley (cross-check), on the self-calibrated smoothed hist
   TH1F *hs = new TH1F("hs", "", 200, -5, 15);
   hs->SetDirectory(nullptr);
   for (int b = 1; b <= h->GetNbinsX(); ++b)
      hs->Fill(h->GetBinCenter(b) - gs, h->GetBinContent(b));
   hs->Smooth(1);
   auto mIn = [&](double lo, double hi, bool mn) {
      double m = mn ? 1e18 : 0;
      for (int b = hs->FindBin(lo); b <= hs->FindBin(hi); ++b)
         m = mn ? std::min(m, hs->GetBinContent(b)) : std::max(m, hs->GetBinContent(b));
      return m;
   };
   double ptv = std::min(mIn(-0.25, 0.25, false), mIn(0.5, 1.0, false)) / std::max(1.0, mIn(0.25, 0.55, true));

   printf("DRES %-24s nD=%ld  gsFWHM=%.3f  doublet_PtV=%.3f  sig=%.3f  (gs@%.3f raw, fitstat=%d)\n", label.Data(), nD,
          fwhm, ptv, sig, gs, st);
}
