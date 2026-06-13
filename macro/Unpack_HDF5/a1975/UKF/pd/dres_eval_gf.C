/// @file dres_eval_gf.C
/// @brief Doublet-resolution metric for GENFIT (p,d) fits — same logic as
/// dres_eval.C but reads the PRA tracks from the parallel _reco.root.
///
/// The genfit output (AtFitterTaskOld) does NOT copy the PRA candidate tracks
/// into the AtTrackingEvent (only the new UKF task does), so the spy/PID gate has
/// nothing to gate on. Here we open BOTH the genfit fit file (kinematics) and the
/// original _reco.root (AtPatternEvent → PRA tracks for the gate), event-matched
/// (FairRunAna writes one tracking event per input event, same order). This keeps
/// the deuteron PID gate IDENTICAL to the UKF dres_eval, so only the fitter differs.
///
///   root -b -q 'pd/dres_eval_gf.C("run_0106,run_0107","genfit",5.0)'

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

void dres_eval_gf(TString runsCSV, TString label, double chi2Cut = 1e9, TString gfDir = "/mnt/f/a1975/reco_gf/",
                  TString recoDir = "/mnt/f/a1975/reco/", TString gateFile = "pid/deuteron_band.json",
                  double Ebeam = 192.0, double bField = 2.85, TString gfPat = "_genfit")
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
      TString gf = gfDir + run + gfPat + ".root"; // _genfit (AtGenfit) or _genfitter (AtGenfitter)
      TString pf = gfDir + run + "_pid.root"; // cached Spyral PID (AtPIDEvent), no spy recompute
      if (gSystem->AccessPathName(gf) || gSystem->AccessPathName(pf))
         continue;
      TFile *fg = TFile::Open(gf);
      TTree *tg = (TTree *)fg->Get("cbmsim");
      TFile *fr = TFile::Open(pf);
      TTree *tr = (TTree *)fr->Get("cbmsim");
      if (!tg || !tr) {
         fg->Close();
         fr->Close();
         continue;
      }
      TClonesArray *te = nullptr, *pe = nullptr;
      tg->SetBranchAddress("AtTrackingEvent", &te);
      tr->SetBranchAddress("AtPIDEvent", &pe);
      Long64_t N = std::min(tg->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         tg->GetEntry(i);
         if (te->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         tr->GetEntry(i);
         if (pe->GetEntries() == 0)
            continue;
         auto *pidev = (AtPIDEvent *)pe->At(0);
         if (!pidev)
            continue;
         std::map<int, const AtTools::AtSpyralResult *> byID; // keyed by stamped trackID
         for (auto &sr : pidev->GetSpyral())
            byID[sr.trackID] = &sr;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            // GenFit kinematics store theta in DEGREES (unlike the UKF's radians);
            // convert for exOf, which expects radians. (Confirmed by C16_pd_anaFit.)
            double c2n = ndf > 0 ? chi2 / ndf : 1e9, ke = k.kineticEnergy, thR = k.theta * TMath::DegToRad();
            // Only require a physical KE + the deuteron PID gate. chi2Cut defaults
            // to "off" (1e9): a fit is a fit — quality cuts belong to the resolution
            // study, not the yield/efficiency. Set chi2Cut explicitly to re-enable.
            if (ke <= 0 || c2n > chi2Cut)
               continue;
            auto it = byID.find(ft->GetTrackID());
            if (it == byID.end())
               continue;
            const auto &r = *it->second; // cached PID, looked up by trackID
            if (!r.valid || !pid.IsInside(r.sqrtdEdx, r.brho))
               continue;
            double ex = exOf(mC16, mp, md, mC15, Ebeam, thR, ke);
            if (!std::isnan(ex))
               h->Fill(ex);
            ++nD;
         }
      }
      fg->Close();
      fr->Close();
   }

   h->Rebin(2);
   TF1 dg("dg", "[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x", -1.5, 2.3);
   double mx = h->GetMaximum();
   dg.SetParameters(0.45, 0.32, 0.6 * mx, 0.6 * mx, 0.1 * mx, 0.0);
   dg.SetParLimits(0, -0.3, 1.0);
   dg.SetParLimits(1, 0.12, 0.55);
   dg.SetParLimits(2, 0, 2 * mx);
   dg.SetParLimits(3, 0, 2 * mx);
   int st = h->Fit(&dg, "QRNS");
   double gs = dg.GetParameter(0), sig = dg.GetParameter(1), fwhm = 2.3548 * sig;

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

   printf("DRES %-20s chi2<%-5g nD=%ld  gsFWHM=%.3f  doublet_PtV=%.3f  sig=%.3f  (gs@%.3f, fitstat=%d)\n", label.Data(),
          chi2Cut, nD, fwhm, ptv, sig, gs, st);
}
