/// @file ex_weighted_3Hed.C
/// @brief The backward-slice Ex spectrum as the EXPERIMENT would see it: DWBA-weighted and
/// normalised to the proposal's own beam time.
///
/// The samples were generated uniform in cos(theta_cm) over 15-80 deg, which is the right thing
/// for measuring acceptance and resolution but is NOT the angular distribution nature produces.
/// This re-weights every event by dsigma/dOmega from ar46_dwba.txt -- the same table whose 15-80
/// integral reproduces the proposal's quoted 1.39 and 1.18 mb -- so no new simulation is needed.
///
/// theta_cm IS RECOMPUTED PER EVENT BY BOOSTING, not inferred from a fixed lookup. The beam loses
/// 0.957 MeV/cm, so the lab-to-CM mapping depends on where in the chamber the reaction happened;
/// a single conversion curve would misassign the weight by several degrees at the chamber ends.
/// The truth deuteron 4-vector is boosted into the CM frame with the beam energy at the TRUE
/// vertex, and theta_cm(DWBA) = 180 - theta_cm(deuteron w.r.t. beam), the convention the table and
/// the proposal both use.
///
/// NORMALISATION carries the acceptance with it. Each state's histogram is scaled by
///     N = L * sigma(15-80) * (sum of weights RECONSTRUCTED in the slice) / (sum of weights GENERATED)
/// so efficiency, the PID cut and the slice selection all enter through the weight ratio rather
/// than as separate factors. L = 500 pps * 7 d * 9.4e20 /cm2 = 2.84e29 /cm2, from the proposal.
///
/// ONLY s1/2 AND d3/2 ARE WEIGHTED. ar46_dwba.txt has no column for the 2.02 MeV 7/2- (l = 3), so
/// that state cannot be put on the same footing and is left out of the weighted spectrum entirely
/// rather than shown with a made-up shape.
///
/// ONE PANEL PER theta_lab SLICE, because the two things that decide whether the states separate
/// pull in opposite directions across the range. Resolution improves steeply toward backward lab
/// angles (pre-fit Ex IQR 3.7 MeV at 58-70 deg, 0.6 MeV at 110-140). Yield does the same thing,
/// since forward theta_cm maps to backward theta_lab in inverse kinematics and the DWBA is
/// forward-peaked. So the backward end wins twice -- but the d3/2 distribution is less
/// forward-peaked than the s1/2 one, so the MIXTURE changes with angle too, and that is visible
/// here as the blue component growing relative to the black one as the angle drops.
///
/// Each panel is normalised to its own expected counts for the proposal's beam time, so panel
/// areas are directly comparable: they are counts, not shapes.
///
///   root -b -q 'ex_weighted_3Hed.C()'

void ex_weighted_3Hed(TString dir = "/mnt/f/ar46_3hed", TString png = "plots/ex_weighted_slices.png",
                      Int_t minPoints = 15,
                      Double_t bField = 2.85, Double_t dThetaMax = 10.0, Double_t driftLength = 100.0,
                      Double_t rate = 500.0, Double_t days = 7.0, Double_t nTarget = 9.4e20, UInt_t seed = 20260813)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957, c = 299.792458;
   const double L = rate * days * 86400.0 * nTarget;

   // ---- DWBA table
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   std::vector<double> th, s12, d32;
   {
      std::ifstream in((here + "/ar46_dwba.txt").Data());
      std::string line;
      while (std::getline(in, line)) {
         std::istringstream is(line);
         double a, b, cc;
         if (!(is >> a >> b >> cc)) continue; // skips the header
         th.push_back(a); s12.push_back(b); d32.push_back(cc);
      }
   }
   if (th.size() < 10) { printf("  could not read ar46_dwba.txt\n"); return; }
   printf("  DWBA table: %zu angles, %.0f to %.0f deg\n", th.size(), th.front(), th.back());
   auto dsdo = [&](double t, int state) { // mb/sr, linear interpolation
      const std::vector<double> &y = (state == 0) ? s12 : d32;
      if (t <= th.front()) return y.front();
      if (t >= th.back()) return y.back();
      size_t i = (size_t)t;
      if (i + 1 >= th.size()) return y.back();
      return y[i] + (y[i + 1] - y[i]) * (t - th[i]) / (th[i + 1] - th[i]);
   };
   auto sigInt = [&](double lo, double hi, int state) { // mb
      double s = 0;
      for (double x = lo; x < hi; x += 0.05)
         s += 2 * TMath::Pi() * std::sin((x + 0.025) * TMath::DegToRad()) * dsdo(x + 0.025, state) *
              0.05 * TMath::DegToRad();
      return s;
   };
   const double sigTot[2] = {sigInt(15, 80, 0), sigInt(15, 80, 1)};
   printf("  sigma(15-80): s1/2 %.3f mb, d3/2 %.3f mb   L = %.3e /cm2\n", sigTot[0], sigTot[1], L);

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) spy.SetMinPoints(minPoints);

   const char *tags[4] = {"gs_s3001", "gs_s3002", "360_s3011", "360_s3012"};
   const int stateOf[4] = {0, 0, 1, 1};
   const int col[2] = {kBlack, kAzure + 2};
   const double edges[] = {58, 70, 83, 96, 110, 140};
   const int NS = 5;
   TH1D *h[2][NS];
   for (int s = 0; s < 2; ++s)
      for (int k = 0; k < NS; ++k) {
         h[s][k] = new TH1D(Form("hw%d_%d", s, k), "", 60, -4.0, 5.0);
         h[s][k]->SetLineColor(col[s]);
         h[s][k]->SetLineWidth(2);
      }
   double wGen[2] = {0, 0}, wSel[2][NS] = {{0}};

   for (int t = 0; t < 4; ++t) {
      int si = stateOf[t];
      TString fs = dir + "/" + tags[t] + "_sim.root", fr = dir + "/" + tags[t] + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) { printf("  skip %s\n", tags[t]); continue; }
      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = (TTree *)Fs->Get("cbmsim"), *tr = (TTree *)Fr->Get("cbmsim");
      TClonesArray *mc = nullptr, *pa = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);
      Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         tr->GetEntry(i);
         double thTrue = -1, Ttrue = 0, zTrue = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0) break;
            thTrue = std::acos(pz / pp) * TMath::RadToDeg();
            Ttrue = std::sqrt(pp * pp + M_e * M_e) - M_e;
            zTrue = p->GetStartZ();
            break;
         }
         if (thTrue < 0) continue;

         // ---- theta_cm by boosting the TRUE deuteron, with the beam energy at the TRUE vertex
         double Tb = Tb0 - dEdz * zTrue;
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double beta = pb / (Eb + M_t), gam = 1.0 / std::sqrt(1 - beta * beta);
         double Ed = Ttrue + M_e, pd = std::sqrt(Ttrue * (Ttrue + 2 * M_e));
         double pzl = pd * std::cos(thTrue * TMath::DegToRad()), ptl = pd * std::sin(thTrue * TMath::DegToRad());
         double pzc = gam * (pzl - beta * Ed);
         double thcm = 180.0 - std::atan2(ptl, pzc) * TMath::RadToDeg(); // DWBA convention
         double w = dsdo(thcm, si);
         wGen[si] += w;

         int sl = -1;
         for (int k = 0; k < NS; ++k)
            if (thTrue >= edges[k] && thTrue < edges[k + 1]) sl = k;
         if (sl < 0) continue;
         if (!pa || !pa->GetEntriesFast()) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         double bd = 1e9, bTh = 0, bT = 0, bZ = 0;
         bool got = false;
         for (auto &track : pe->GetTrackCand()) {
            auto res = spy.Estimate(const_cast<AtTrack &>(track));
            if (!res.valid) continue;
            double thr = 180.0 - res.polar * TMath::RadToDeg();
            double p = c * res.brho;
            if (std::fabs(thr - thTrue) < bd) {
               bd = std::fabs(thr - thTrue);
               bTh = thr; bT = std::sqrt(p * p + M_e * M_e) - M_e; bZ = res.vertex.Z() / 10.0; got = true;
            }
         }
         if (!got || bd > dThetaMax) continue;
         double TbR = Tb0 - dEdz * (driftLength - bZ);
         if (TbR < 50 || TbR > Tb0 + 20) continue;
         double EbR = TbR + M_b, pbR = std::sqrt(TbR * (TbR + 2 * M_b));
         double EdR = bT + M_e, pdR = std::sqrt(bT * (bT + 2 * M_e));
         double a = bTh * TMath::DegToRad();
         double ER = EbR + M_t - EdR, pRz = pbR - pdR * std::cos(a), pRt = pdR * std::sin(a);
         double m2 = ER * ER - pRz * pRz - pRt * pRt;
         if (m2 <= 0) continue;
         h[si][sl]->Fill(std::sqrt(m2) - M_R, w);
         wSel[si][sl] += w;
      }
      Fs->Close();
      Fr->Close();
      printf("  %-12s done\n", tags[t]);
   }

   // ---- normalise each slice to the proposal's beam time; acceptance rides in the weight ratio
   const char *nm[2] = {"s1/2 g.s.", "d3/2 0.36"};
   double Nexp[2][NS];
   printf("\n  %.0f pps x %.0f days,  L = %.3e /cm2\n", rate, days, L);
   printf("  theta_lab      s1/2 counts   IQR      d3/2 counts   IQR     sep/sigma\n");
   for (int k = 0; k < NS; ++k) {
      double iqr[2] = {0, 0};
      for (int s = 0; s < 2; ++s) {
         double acc = wSel[s][k] / std::max(1e-9, wGen[s]);
         Nexp[s][k] = L * sigTot[s] * 1e-27 * acc;
         if (h[s][k]->Integral() > 0) {
            double qp[3] = {0.25, 0.50, 0.75}, qv[3];
            h[s][k]->GetQuantiles(3, qv, qp);
            iqr[s] = qv[2] - qv[0];
            h[s][k]->Scale(Nexp[s][k] / h[s][k]->Integral());
         }
      }
      double sig = 0.5 * (iqr[0] + iqr[1]) / 1.349;
      printf("  %3.0f - %3.0f deg %10.1f %7.2f %13.1f %7.2f %10s\n", edges[k], edges[k + 1], Nexp[0][k], iqr[0],
             Nexp[1][k], iqr[1], sig > 0 ? Form("%.2f", 0.36 / sig) : "--");
   }

   TCanvas *cv = new TCanvas("cW", "weighted, per angle slice", 1500, 900);
   cv->Divide(3, 2);
   gRandom->SetSeed(seed);
   for (int k = 0; k < NS; ++k) {
      cv->cd(k + 1);
      auto *hs = new THStack(Form("hs%d", k),
                             Form("#theta_{lab} %.0f - %.0f deg;E_{x} [MeV];counts / %.0f keV", edges[k],
                                  edges[k + 1], 1000 * h[0][k]->GetBinWidth(1)));
      h[0][k]->SetFillColorAlpha(kBlack, 0.15);
      h[1][k]->SetFillColorAlpha(kAzure + 2, 0.25);
      hs->Add(h[1][k]);
      hs->Add(h[0][k]);
      // one Poisson realisation of the SUM: what a single run shows in this slice
      auto *toy = (TH1D *)h[0][k]->Clone(Form("toy%d", k));
      toy->Reset();
      toy->SetFillStyle(0);
      toy->SetLineColor(kRed + 1);
      toy->SetMarkerStyle(20);
      toy->SetMarkerSize(0.6);
      toy->SetMarkerColor(kRed + 1);
      for (int b = 1; b <= toy->GetNbinsX(); ++b) {
         double n = gRandom->Poisson(h[0][k]->GetBinContent(b) + h[1][k]->GetBinContent(b));
         toy->SetBinContent(b, n);
         toy->SetBinError(b, std::sqrt(n));
      }
      // The y-range must cover BOTH the smooth expectation and the integer realisation. In the
      // low-count slices the expectation peaks below 1 count/bin while a single run puts 2 or 3 in
      // some bin, and a range set by the stack alone turns those into spikes off the top.
      hs->Draw("hist");
      hs->SetMaximum(1.35 * std::max(hs->GetMaximum("nostack"), toy->GetMaximum()));
      hs->Draw("hist");
      toy->Draw("E same");
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.045);
      tx->DrawLatex(0.55, 0.84, Form("s1/2 %.0f", Nexp[0][k]));
      tx->DrawLatex(0.55, 0.78, Form("d3/2 %.0f", Nexp[1][k]));
   }
   cv->cd(6);
   auto *lg = new TLegend(0.05, 0.30, 0.95, 0.80);
   lg->SetBorderSize(0);
   lg->SetHeader("DWBA-weighted, 7 d at 500 pps");
   lg->AddEntry(h[0][0], "s1/2 g.s. (l = 0)", "f");
   lg->AddEntry(h[1][0], "d3/2 0.36 MeV (l = 2)", "f");
   auto *dum = new TH1D("dum", "", 1, 0, 1);
   dum->SetMarkerStyle(20);
   dum->SetMarkerColor(kRed + 1);
   dum->SetLineColor(kRed + 1);
   lg->AddEntry(dum, "one 7-day run (Poisson)", "pe");
   lg->Draw();
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   cv->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());
}
