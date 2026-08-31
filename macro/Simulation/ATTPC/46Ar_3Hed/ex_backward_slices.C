/// @file ex_backward_slices.C
/// @brief Backward-angle Ex spectra, DWBA-weighted, at a chosen slice width. Three states.
///
/// Two things this separates that the 13-30 deg-wide slices could not.
///
/// 1. RESOLUTION vs SLICE WIDTH. The Ex centroid drifts with angle (the g.s. runs +1.02 MeV at
///    58-70 deg to -0.26 at 110-140), so a wide slice adds that drift to the true resolution and
///    reports the sum as width. Run this at 10 and then at 5 deg: if the IQR shrinks, the earlier
///    number was inflated by the drift and the real resolution is better than it looked. If it
///    does not move, the width is genuine.
///
/// 2. THE 2.02 MeV 7/2- STATE, which the 15-80 deg picture left out entirely.
///
/// !! THE 2.02 NORMALISATION IS ASSUMED, NOT DERIVED !! ar46_dwba.txt has columns for s1/2 and
/// d3/2 only. For the 7/2- (l = 3) this macro uses the d3/2 ANGULAR SHAPE as a stand-in and the
/// d3/2 INTEGRATED CROSS SECTION as its normalisation, which is what "comparable statistics"
/// means here -- enough counts to see whether the peak separates, not a prediction of its yield.
/// A real l = 3 distribution is less forward-peaked than l = 2, so this proxy puts too much of the
/// 7/2- at backward angles: read its COUNTS as illustrative and its POSITION and WIDTH as real,
/// since those come from the simulation and not from the weight.
///
///   root -b -q 'ex_backward_slices.C(10)'    // 4 panels, 100-140 deg
///   root -b -q 'ex_backward_slices.C(5)'     // 8 panels, same range

void ex_backward_slices(Double_t width = 10.0, Double_t thLo = 100.0, Double_t thHi = 140.0,
                        TString dir = "/mnt/f/ar46_3hed_OLD_2.85T_placeholder", TString png = "", Int_t minPoints = 15,
                        Double_t bField = 2.85, Double_t dThetaMax = 10.0, Double_t driftLength = 100.0,
                        Double_t rate = 500.0, Double_t days = 7.0, Double_t nTarget = 9.4e20,
                        UInt_t seed = 20260813)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957, c = 299.792458;
   const double L = rate * days * 86400.0 * nTarget;
   const int NS = (int)std::lround((thHi - thLo) / width);
   if (NS < 1 || NS > 12) { printf("  %d slices is not a sensible number\n", NS); return; }
   if (png.IsNull()) png = TString::Format("plots/ex_backward_%.0fdeg.png", width);

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   std::vector<double> th, s12, d32;
   {
      std::ifstream in((here + "/ar46_dwba.txt").Data());
      std::string line;
      while (std::getline(in, line)) {
         std::istringstream is(line);
         double a, b, cc;
         if (!(is >> a >> b >> cc)) continue;
         th.push_back(a); s12.push_back(b); d32.push_back(cc);
      }
   }
   if (th.size() < 10) { printf("  could not read ar46_dwba.txt\n"); return; }
   auto dsdo = [&](double t, int col) {
      const std::vector<double> &y = (col == 0) ? s12 : d32;
      if (t <= th.front()) return y.front();
      if (t >= th.back()) return y.back();
      size_t i = (size_t)t;
      if (i + 1 >= th.size()) return y.back();
      return y[i] + (y[i + 1] - y[i]) * (t - th[i]) / (th[i + 1] - th[i]);
   };
   auto sigInt = [&](double lo, double hi, int col) {
      double s = 0;
      for (double x = lo; x < hi; x += 0.05)
         s += 2 * TMath::Pi() * std::sin((x + 0.025) * TMath::DegToRad()) * dsdo(x + 0.025, col) * 0.05 *
              TMath::DegToRad();
      return s;
   };
   // state -> which DWBA column weights it, and what integrated cross section normalises it
   const int wCol[3] = {0, 1, 1};                                        // 7/2- uses the d3/2 shape (proxy)
   const double sigTot[3] = {sigInt(15, 80, 0), sigInt(15, 80, 1), sigInt(15, 80, 1)};
   const char *nm[3] = {"s1/2 g.s.", "d3/2 0.36", "f7/2 2.02"};
   const int col[3] = {kBlack, kAzure + 2, kRed + 1};
   const double exTrue[3] = {0.0, 0.360, 2.020};
   printf("  sigma(15-80): s1/2 %.3f, d3/2 %.3f mb   (7/2- normalised to the d3/2 value -- assumed)\n", sigTot[0],
          sigTot[1]);

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) spy.SetMinPoints(minPoints);

   const char *tags[6] = {"gs_s3001", "gs_s3002", "360_s3011", "360_s3012", "2020_s3021", "2020_s3022"};
   std::vector<std::vector<TH1D *>> h(3, std::vector<TH1D *>(NS, nullptr));
   for (int s = 0; s < 3; ++s)
      for (int k = 0; k < NS; ++k) {
         h[s][k] = new TH1D(Form("hb%d_%d", s, k), "", 50, -2.0, 4.0);
         h[s][k]->SetLineColor(col[s]);
         h[s][k]->SetLineWidth(2);
      }
   std::vector<double> wGen(3, 0.0);
   std::vector<std::vector<double>> wSel(3, std::vector<double>(NS, 0.0));
   std::vector<std::vector<long>> nRaw(3, std::vector<long>(NS, 0)); // simulated tracks, for the noise check

   for (int t = 0; t < 6; ++t) {
      int si = t / 2;
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

         double Tb = Tb0 - dEdz * zTrue;
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double beta = pb / (Eb + M_t), gam = 1.0 / std::sqrt(1 - beta * beta);
         double Ed = Ttrue + M_e, pd = std::sqrt(Ttrue * (Ttrue + 2 * M_e));
         double pzl = pd * std::cos(thTrue * TMath::DegToRad()), ptl = pd * std::sin(thTrue * TMath::DegToRad());
         double thcm = 180.0 - std::atan2(ptl, gam * (pzl - beta * Ed)) * TMath::RadToDeg();
         double w = dsdo(thcm, wCol[si]);
         wGen[si] += w;

         if (thTrue < thLo || thTrue >= thHi) continue;
         int sl = (int)((thTrue - thLo) / width);
         if (sl < 0 || sl >= NS) continue;
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
         nRaw[si][sl]++;
      }
      Fs->Close();
      Fr->Close();
   }

   printf("\n  %.0f deg slices, %.0f-%.0f deg,  %.0f pps x %.0f d\n", width, thLo, thHi, rate, days);
   printf("  slice        state       counts   median     IQR   (sim tracks)\n");
   std::vector<std::vector<double>> Nexp(3, std::vector<double>(NS, 0.0)), iqr(3, std::vector<double>(NS, 0.0)),
      med(3, std::vector<double>(NS, 0.0));
   for (int k = 0; k < NS; ++k) {
      for (int s = 0; s < 3; ++s) {
         double acc = wSel[s][k] / std::max(1e-9, wGen[s]);
         Nexp[s][k] = L * sigTot[s] * 1e-27 * acc;
         if (h[s][k]->Integral() > 0) {
            double qp[3] = {0.25, 0.50, 0.75}, qv[3];
            h[s][k]->GetQuantiles(3, qv, qp);
            med[s][k] = qv[1];
            iqr[s][k] = qv[2] - qv[0];
            h[s][k]->Scale(Nexp[s][k] / h[s][k]->Integral());
         }
         printf("  %3.0f-%3.0f deg   %-10s %7.1f  %+7.3f %7.3f   %6ld%s\n", thLo + k * width, thLo + (k + 1) * width,
                nm[s], Nexp[s][k], med[s][k], iqr[s][k], nRaw[s][k], nRaw[s][k] < 100 ? "  <- noisy" : "");
      }
      printf("\n");
   }

   int nx = (NS <= 4) ? 2 : 4, ny = (NS + nx - 1) / nx;
   TCanvas *cv = new TCanvas("cB", "backward slices", 400 * nx, 380 * ny);
   cv->Divide(nx, ny);
   gRandom->SetSeed(seed);
   for (int k = 0; k < NS; ++k) {
      cv->cd(k + 1);
      auto *hs = new THStack(Form("hsb%d", k),
                             Form("#theta_{lab} %.0f - %.0f deg;E_{x} [MeV];counts / %.0f keV", thLo + k * width,
                                  thLo + (k + 1) * width, 1000 * h[0][k]->GetBinWidth(1)));
      h[0][k]->SetFillColorAlpha(kBlack, 0.15);
      h[1][k]->SetFillColorAlpha(kAzure + 2, 0.25);
      h[2][k]->SetFillColorAlpha(kRed + 1, 0.20);
      hs->Add(h[2][k]);
      hs->Add(h[1][k]);
      hs->Add(h[0][k]);
      auto *toy = (TH1D *)h[0][k]->Clone(Form("tb%d", k));
      toy->Reset();
      toy->SetFillStyle(0);
      toy->SetLineColor(kGray + 2);
      toy->SetMarkerStyle(20);
      toy->SetMarkerSize(0.6);
      toy->SetMarkerColor(kGray + 2);
      for (int b = 1; b <= toy->GetNbinsX(); ++b) {
         double n = gRandom->Poisson(h[0][k]->GetBinContent(b) + h[1][k]->GetBinContent(b) +
                                     h[2][k]->GetBinContent(b));
         toy->SetBinContent(b, n);
         toy->SetBinError(b, std::sqrt(n));
      }
      hs->Draw("hist");
      hs->SetMaximum(1.35 * std::max(hs->GetMaximum("nostack"), toy->GetMaximum()));
      hs->Draw("hist");
      toy->Draw("E same");
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.05);
      for (int s = 0; s < 3; ++s) {
         tx->SetTextColor(col[s]);
         tx->DrawLatex(0.60, 0.84 - 0.06 * s, Form("%.0f  (%.2f)", Nexp[s][k], iqr[s][k]));
      }
      tx->SetTextColor(kBlack);
      tx->SetTextSize(0.038);
      tx->DrawLatex(0.60, 0.63, "counts (IQR)");
   }
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   cv->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());
}
