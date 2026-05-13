/// @file phi_at_xref.C
/// @brief Compute σ_φ when the UKF angle is evaluated at a FIXED reference
/// x along the beam, instead of at the (event-dependent) first cluster.
///
/// Method:
///   - UKF tangent: linear fit of y_smoothed vs x_smoothed on points within
///     ±halfWin of x_ref → slope dy/dx → φ_UKF = atan(slope)
///   - MC tangent: AtMCPoint in drift_volume closest to x_ref →
///     φ_MC = atan2(Py, Px)
///
/// If σ_φ_ref << σ_φ at first cluster, the broadening is caused by the
/// spread in first-cluster x (× local curvature).
///
/// Output: data/phi_at_xref_HYDRA.png + per-p σ to stdout.

void phi_at_xref(double x_ref = 50., double halfWin = 30.,
                 const char *outPng = "data/phi_at_xref_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.055, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.135);
   gStyle->SetPadBottomMargin(0.135);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.085);

   const std::vector<int> plist = {200, 400, 600, 800, 1000, 1200};
   const int nP = plist.size();
   const double radToMrad = 1000.;

   std::vector<double> pV(nP), sigPhiRef(nP), biasPhiRef(nP);
   std::vector<double> sigPhiFC(nP);   // first-cluster baseline for comparison

   auto *c = new TCanvas("c", "σ_φ vs anchor", 1200, 700);
   c->Divide(2, 1, 0.012, 0.015);

   // Panel A: σ_φ vs p, two curves (at first cluster vs at fixed x_ref)
   auto *gFC = new TGraph(nP);
   auto *gRef = new TGraph(nP);
   gFC->SetMarkerStyle(20); gFC->SetMarkerSize(1.5);
   gFC->SetMarkerColor(kRed + 1); gFC->SetLineColor(kRed + 1); gFC->SetLineWidth(2);
   gRef->SetMarkerStyle(20); gRef->SetMarkerSize(1.5);
   gRef->SetMarkerColor(kBlue + 1); gRef->SetLineColor(kBlue + 1); gRef->SetLineWidth(2);

   // Panel B: residual histograms at x_ref overlaid
   const int colors[] = {kBlue + 1, kRed + 1, kGreen + 3, kMagenta + 1, kOrange + 7, kCyan + 2};
   std::vector<TH1F *> hList(nP);

   for (int k = 0; k < nP; ++k) {
      int P = plist[k];
      pV[k] = P;
      TFile fS(Form("data/HYDRAsim_p%d.root", P));
      TFile fU(Form("data/output_ukf_HYDRA_p%d.root", P));
      auto *tS = (TTree *)fS.Get("cbmsim");
      auto *tU = (TTree *)fU.Get("cbmsim");
      auto *trks = new TClonesArray("AtMCTrack");
      auto *pts  = new TClonesArray("AtMCPoint");
      auto *te = new TClonesArray("AtTrackingEvent");
      tS->SetBranchAddress("MCTrack", &trks);
      tS->SetBranchAddress("AtTpcPoint", &pts);
      tU->SetBranchAddress("AtTrackingEvent", &te);

      TH1F hRef(Form("hRef_%d", P), "", 80, -150., 150.);
      TH1F hFC(Form("hFC_%d", P), "", 80, -150., 150.);
      hRef.SetDirectory(nullptr); hFC.SetDirectory(nullptr);

      Long64_t n = std::min(tS->GetEntries(), tU->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         tS->GetEntry(i);
         tU->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;

         // MC angle: AtMCPoint nearest x_ref
         double phMC = std::numeric_limits<double>::quiet_NaN();
         double phMCfirst = std::numeric_limits<double>::quiet_NaN();
         double dMin = 1e9, xFirstMin = 1e9;
         for (int j = 0; j < pts->GetEntries(); ++j) {
            auto *p = (AtMCPoint *)pts->At(j);
            if (p->GetVolName() != TString("drift_volume")) continue;
            if (p->GetTrackID() != 0) continue;
            double x_mm = p->GetX() * 10.;
            double pp = std::hypot(p->GetPx(), std::hypot(p->GetPy(), p->GetPz()));
            if (pp <= 0) continue;
            double phi = std::atan2(p->GetPy(), p->GetPx());
            if (std::abs(x_mm - x_ref) < dMin) {
               dMin = std::abs(x_mm - x_ref);
               phMC = phi;
            }
            if (x_mm < xFirstMin) { xFirstMin = x_mm; phMCfirst = phi; }
         }
         if (std::isnan(phMC) || std::isnan(phMCfirst)) continue;

         if (te->GetEntries() == 0) continue;
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (!best) continue;

         // UKF angle at first cluster (baseline)
         double phUKFfc = best->GetKinematics().phi;

         // UKF angle at x_ref: local linear fit of y(x) on smoothed positions
         auto &sp = best->GetSmoothedPositions();
         double sumW = 0, sumX = 0, sumY = 0, sumXX = 0, sumXY = 0;
         for (auto &pp : sp) {
            double x = pp.X(), y = pp.Y();
            if (std::abs(x - x_ref) > halfWin) continue;
            sumW += 1.; sumX += x; sumY += y;
            sumXX += x * x; sumXY += x * y;
         }
         if (sumW < 4) continue; // need enough points near x_ref
         double xbar = sumX / sumW, ybar = sumY / sumW;
         double Sxx = sumXX - sumW * xbar * xbar;
         double Sxy = sumXY - sumW * xbar * ybar;
         if (Sxx <= 0) continue;
         double slope = Sxy / Sxx;
         // φ = direction of motion in xy. For motion mostly +x, φ ≈ atan(slope).
         double phUKFref = std::atan(slope);

         auto wrap = [](double d) {
            while (d > M_PI) d -= 2 * M_PI;
            while (d < -M_PI) d += 2 * M_PI;
            return d;
         };
         hRef.Fill(wrap(phUKFref - phMC) * radToMrad);
         hFC.Fill(wrap(phUKFfc - phMCfirst) * radToMrad);
      }

      auto fitCore = [](TH1F &h, double &mu, double &sig) {
         double rms = h.GetRMS();
         TFitResultPtr fr = h.Fit("gaus", "SQR0", "", -3 * rms, 3 * rms);
         mu = fr.Get() ? fr->Parameter(1) : h.GetMean();
         sig = fr.Get() ? fr->Parameter(2) : rms;
      };
      double mu, sig, muFC, sigFC;
      fitCore(hRef, mu, sig);
      fitCore(hFC, muFC, sigFC);
      sigPhiRef[k] = sig; biasPhiRef[k] = mu;
      sigPhiFC[k] = sigFC;
      gFC->SetPoint(k, P, sigFC);
      gRef->SetPoint(k, P, sig);

      auto *h = (TH1F *)hRef.Clone(Form("hRef_keep_%d", P));
      h->SetDirectory(nullptr);
      h->SetLineColor(colors[k % 6]); h->SetLineWidth(2);
      hList[k] = h;

      std::cout << "p=" << P << " MeV/c"
                << "  σ_φ@first_cluster=" << sigFC << " mrad (bias=" << muFC << ")"
                << "  σ_φ@x_ref=" << x_ref << "=" << sig << " mrad (bias=" << mu << ")\n";
   }

   // Panel A
   c->cd(1); gPad->SetGrid();
   gFC->SetTitle(Form("#sigma_{#varphi} vs p: anchor comparison;p_{MC} (MeV/c);#sigma_{#varphi} (mrad)"));
   gFC->Draw("APL");
   double mx = std::max(*std::max_element(sigPhiFC.begin(), sigPhiFC.end()),
                         *std::max_element(sigPhiRef.begin(), sigPhiRef.end()));
   gFC->GetYaxis()->SetRangeUser(0., 1.3 * mx);
   gFC->GetXaxis()->SetLimits(0., 1400.);
   gRef->Draw("PL same");
   auto *leg = new TLegend(0.40, 0.65, 0.93, 0.88);
   leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
   leg->AddEntry(gFC, "at first cluster (UKF default)", "lp");
   leg->AddEntry(gRef, Form("at x_{ref} = %.0f mm (#pm%.0f mm window)", x_ref, halfWin), "lp");
   leg->Draw();

   // Panel B
   c->cd(2); gPad->SetGrid();
   double maxY = 0;
   for (auto *h : hList) if (h && h->Integral() > 0) {
      h->Scale(1.0 / h->Integral());
      maxY = std::max(maxY, h->GetMaximum());
   }
   for (int k = 0; k < nP; ++k) {
      if (!hList[k]) continue;
      hList[k]->SetMaximum(1.3 * maxY);
      hList[k]->SetTitle(Form("#Delta#varphi at x_{ref} = %.0f mm;#Delta#varphi (mrad);fraction", x_ref));
      hList[k]->Draw(k == 0 ? "HIST" : "HIST same");
   }
   auto *leg2 = new TLegend(0.62, 0.55, 0.93, 0.92);
   leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.04);
   for (int k = 0; k < nP; ++k) if (hList[k]) leg2->AddEntry(hList[k], Form("%d MeV/c", plist[k]), "l");
   leg2->Draw();

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
