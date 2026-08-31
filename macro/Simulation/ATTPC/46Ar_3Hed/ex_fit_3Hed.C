/// @file ex_fit_3Hed.C
/// @brief Ex(47K) from the GENFIT-FITTED tracks, against the pre-fit floor from kinematics_3Hed.C.
///
/// Same two-body inversion as kinematics_3Hed.C, same masses, same beam-energy-at-vertex
/// treatment -- the ONLY thing that changes is where (theta, T_d) comes from: the genfit
/// AtFittedTrack instead of the AtSpyralPID circle. That is deliberate, because a resolution
/// comparison is worthless if the two sides differ in anything else.
///
///     E_R = E_beam(z_vertex) + M(3He) - E_d,   p_R = p_beam - p_d,   Ex = sqrt(E_R^2 - p_R^2) - M(47K)
///
/// VERTEX COMES FROM GetVertex(), NOT FROM TrackProperties. EventFit::AtGenfitter leaves
/// TrackProperties::initialPositionXtr at zero -- measured, <z> = 0.00 cm over 5098 fitted tracks
/// -- while GetVertex() is filled and sensible (<z> = 46.9 cm, mid-chamber). Using the empty one
/// silently pins the beam energy to one end of the chamber and inflates Ex.
///
/// THE HANDEDNESS IS DECIDED FROM THE DATA, not inherited. AtSpyralPID's vertex comes back
/// mirrored (r = -1.000 against truth, kinematics_3Hed.C), but the fitter runs its own
/// extrapolation with SetZPadPlane(1000), so it need not share that convention. This macro
/// correlates the fitted vertex against truth, prints r, and mirrors only if r < 0.
///
/// GENFIT USES THE RAW KINEMATICS SLOT. GetKinematics() is the raw slot for genfit and the
/// corrected one for UKF -- do not swap this for GetKinematicsXtr() without re-reading which
/// fitter wrote the file. theta is in RADIANS here (EventFit::AtGenfitter).
///
///   root -b -q 'ex_fit_3Hed.C("gs_s3001")'
///   root -b -q 'ex_fit_3Hed.C("gs_s3001",0.0)'    // second arg = the state's true Ex

/// zShift (cm) is a DIAGNOSTIC, not a calibration: the fitted vertex sits ~3.5 cm from truth, and
/// at 0.957 MeV/cm that is 3.3 MeV of beam energy, the right size to explain the Ex offset. Pass
/// the measured bias to test that attribution. Do not leave it non-zero in production -- it tunes
/// the answer on the truth.
void ex_fit_3Hed(TString tag = "gs_s3001", Double_t exTrue = 0.0, TString dir = "/mnt/f/ar46_3hed_OLD_2.85T_placeholder",
                 TString png = "plots/ex_fit_3Hed.png", Double_t dThetaMax = 10.0, Double_t driftLength = 100.0,
                 Double_t zShift = 0.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957;

   TString fs = dir + "/" + tag + "_sim.root", ff = dir + "/" + tag + "_genfitter_d.root";
   if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(ff)) {
      printf("  missing %s or %s\n", fs.Data(), ff.Data());
      return;
   }
   TFile *Fs = TFile::Open(fs), *Ff = TFile::Open(ff);
   TTree *ts = (TTree *)Fs->Get("cbmsim"), *tf = (TTree *)Ff->Get("cbmsim");
   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   auto *hEx = new TH1D("hExF", TString::Format("E_{x}(^{47}K) from the fit, %s;E_{x} [MeV];tracks", tag.Data()),
                        200, -3, 6);
   auto *hdT = new TH1D("hdT", "fitted - true KE;#DeltaT_{d} [MeV];tracks", 200, -20, 20);
   auto *hdA = new TH1D("hdA", "fitted - true angle;#Delta#theta [deg];tracks", 200, -20, 20);
   auto *h2A = new TH2D("h2A", "angle: fitted vs true;#theta_{true} [deg];#theta_{fit} [deg]", 90, 50, 140, 90, 50, 140);
   auto *h2K = new TH2D("h2K", "energy: fitted vs true;T_{true} [MeV];T_{fit} [MeV]", 80, 0, 70, 80, 0, 70);
   auto *h2E = new TH2D("h2E", "E_{x} error vs angle -- where the width comes from;#theta_{lab} [deg];E_{x} [MeV]",
                        90, 50, 140, 90, -6, 8);
   double ssum = 0, ssum2 = 0, svtx = 0, sxtr = 0, cx = 0, cy = 0, cxy = 0, cxx = 0, cyy = 0;
   std::vector<double> vz, vzt, kz, kth, kT;
   bool mirror = false;
   long nv = 0, nTruth = 0, nFit = 0;

   Long64_t N = std::min(ts->GetEntries(), tf->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      tf->GetEntry(i);
      double thTrue = -1, Ttrue = 0, zTrue = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
            continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pp = std::sqrt(px * px + py * py + pz * pz);
         if (pp <= 0) break;
         thTrue = std::acos(pz / pp) * TMath::RadToDeg();
         Ttrue = std::sqrt(pp * pp + M_e * M_e) - M_e;
         zTrue = p->GetStartZ();
         break;
      }
      if (thTrue < 0) continue;
      ++nTruth;
      if (!te || !te->GetEntriesFast()) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;

      double bd = 1e9, bTh = 0, bT = 0, bZ = 0, bZv = 0;
      bool got = false;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue;
         auto &k = ft->GetKinematics();
         if (k.kineticEnergy <= 0) continue;
         double th = k.theta * TMath::RadToDeg(); // radians for EventFit::AtGenfitter
         double d = std::fabs(th - thTrue);
         if (d < bd) {
            bd = d; bTh = th; bT = k.kineticEnergy;
            bZ = ft->GetVertex().Z() / 10.0;                                    // mm -> cm, the filled slot
            bZv = ft->GetTrackPropertiesStruct().initialPositionXtr.Z() / 10.0; // empty for this fitter; kept as the check
            got = true;
         }
      }
      if (!got) continue;
      ssum += zTrue + bZ; ssum2 += (zTrue + bZ) * (zTrue + bZ); ++nv;
      svtx += bZv; sxtr += bZ;
      cx += zTrue; cy += bZ; cxy += zTrue * bZ; cxx += zTrue * zTrue; cyy += bZ * bZ;
      vz.push_back(bZ); vzt.push_back(zTrue);
      if (bd > dThetaMax) continue;
      ++nFit;
      hdT->Fill(bT - Ttrue);
      hdA->Fill(bTh - thTrue);
      h2A->Fill(thTrue, bTh);
      h2K->Fill(Ttrue, bT);

      kz.push_back(bZ); kth.push_back(bTh); kT.push_back(bT);
   }

   // ---- handedness of the FITTED vertex, decided before it is used
   double rV = (nv * cxy - cx * cy) / std::sqrt((nv * cxx - cx * cx) * (nv * cyy - cy * cy));
   mirror = (rV < 0);
   for (size_t j = 0; j < kz.size(); ++j) {
      double zUse = (mirror ? driftLength - kz[j] : kz[j]) + zShift;
      double Tb = Tb0 - dEdz * zUse;
      if (Tb < 50 || Tb > Tb0 + 20) continue;
      double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
      double Ed = kT[j] + M_e, pd = std::sqrt(kT[j] * (kT[j] + 2 * M_e));
      double th = kth[j] * TMath::DegToRad();
      double ER = Eb + M_t - Ed;
      double pRz = pb - pd * std::cos(th), pRt = pd * std::sin(th);
      double m2 = ER * ER - pRz * pRz - pRt * pRt;
      if (m2 <= 0) continue;
      hEx->Fill(std::sqrt(m2) - M_R);
      h2E->Fill(kth[j], std::sqrt(m2) - M_R);
   }

   double mSum = ssum / std::max(1L, nv);
   printf("\n  %s: truth %ld, fitted+matched %ld\n", tag.Data(), nTruth, nFit);
   printf("  vertex sources: GetVertex() <z> = %.2f cm (used),  initialPositionXtr <z> = %.2f cm%s\n",
          sxtr / std::max(1L, nv), svtx / std::max(1L, nv),
          std::fabs(svtx / std::max(1L, nv)) < 1e-6 ? "   <-- not filled by this fitter" : "");
   printf("  fitted vertex vs truth: r = %+.3f  ->  %s\n", rV,
          mirror ? "REVERSED, mirrored before use" : "SAME sense, used as is");
   printf("  z_true %s z_fit = %.2f cm%s\n", mirror ? "+" : "-", mirror ? mSum : (cx - cy) / std::max(1L, nv),
          mirror && std::fabs(mSum - driftLength) > 3.0 ? "   <-- expected the drift length; Ex not trustworthy" : "");
   printf("  dTheta: mean %+.2f  rms %.2f deg     dKE: mean %+.2f  rms %.2f MeV\n", hdA->GetMean(), hdA->GetRMS(),
          hdT->GetMean(), hdT->GetRMS());

   // FWHM BY GAUSSIAN FIT, not by walking bins from the maximum. The bin-walk estimator this
   // replaced truncates on a single downward fluctuation next to the peak, which made it report
   // 0.765 MeV and 1.800 MeV for two histograms that differ only by a CONSTANT shift of every
   // entry -- a shift cannot change a width, so the estimator was measuring noise. Smoothed
   // half-height crossings set the fit window; the fit itself runs on the unsmoothed histogram.
   auto stat = [&](TH1D *h, const char *lab, double truth) {
      if (h->GetEntries() < 50) { printf("  %-10s too few entries\n", lab); return; }
      TH1D *hs = (TH1D *)h->Clone(TString(h->GetName()) + "_s");
      hs->Smooth(3);
      int b = hs->GetMaximumBin();
      double peak = hs->GetBinCenter(b), half = hs->GetMaximum() / 2, lo = peak, hi = peak;
      for (int k = b; k > 1 && hs->GetBinContent(k) > half; --k) lo = hs->GetBinCenter(k);
      for (int k = b; k < hs->GetNbinsX() && hs->GetBinContent(k) > half; ++k) hi = hs->GetBinCenter(k);
      double w = std::max(0.15, hi - lo);
      TF1 g("g", "gaus", peak - w, peak + w);
      g.SetParameters(h->GetMaximum(), peak, w / 2.355);
      h->Fit(&g, "QNR");
      double mu = g.GetParameter(1), fw = 2.3548 * std::fabs(g.GetParameter(2));
      double qp[3] = {0.25, 0.50, 0.75}, qv[3];
      h->GetQuantiles(3, qv, qp);
      printf("  %-10s median %+6.3f (true %+.2f)   IQR %.3f MeV   gaus-peak %+6.3f FWHM %.3f   entries %.0f\n", lab,
             qv[1], truth, qv[2] - qv[0], mu, fw, h->GetEntries());
      delete hs;
   };
   printf("\n  FITTED Ex:\n");
   stat(hEx, "genfit", exTrue);

   TCanvas *c = new TCanvas("cF", "fitted Ex", 1500, 900);
   c->Divide(3, 2);
   gStyle->SetPalette(kBird);
   c->cd(1); gPad->SetRightMargin(0.13); gPad->SetLogz(); h2A->Draw("colz");
   c->cd(2); gPad->SetRightMargin(0.13); gPad->SetLogz(); h2K->Draw("colz");
   c->cd(3); gPad->SetRightMargin(0.13); gPad->SetLogz(); h2E->Draw("colz");
   c->cd(4); hdA->Draw("hist");
   c->cd(5); hdT->Draw("hist");
   c->cd(6); hEx->Draw("hist");
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());
}
