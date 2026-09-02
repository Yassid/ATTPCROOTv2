/// @file dp_kinematics_C17.C
/// @brief Analytic 17C(d,p)18C kinematics: where the protons go, how fast, and how much leverage
///        the excitation energy has on each observable. No simulation input -- this is the map the
///        simulation is measured against.
///
///   root -b -q 'dp_kinematics_C17.C'
///   root -b -q 'dp_kinematics_C17.C(135.0)'    // at the mean vertex beam energy
///
/// WHY IT MATTERS FOR THE BINNING. ex_res_C14_hf.C hard-codes theta_lab slices 20-90 deg, because
/// it was written for 14C(p,p') where the recoil proton follows theta_lab = (180-theta_cm)/2 and
/// 90 deg IS the physical limit. A (d,p) proton has no such limit. Before trusting any per-slice
/// resolution table for this channel, the mapping below has to say which lab angles the transfer
/// peak actually occupies.
///
/// THE COMPARISON THAT MATTERS. 14C(d,p)15C has Q = -1.007 MeV on an 11.5 MeV/u beam; this channel
/// has Q = +1.959 MeV on an 8.37 MeV/u beam. The higher Q pushes the proton energy UP and the
/// slower beam pushes it DOWN, so the two channels are not related by a simple shift and the
/// 14C(d,p) result (backward sigma(Ex) flat at ~0.2 MeV, all the field/pitch gain forward) cannot
/// be carried over by assumption.
///
/// dEx/dKE and dEx/dtheta are the leverage terms: they say how an error in the measured proton
/// energy or angle propagates into the excitation energy, which is what sets the resolution
/// budget. They are computed numerically from the same inversion the analysis uses.
/// @param outDir  where to write kinematics_C17dp.png. Empty disables plotting.
/// @param simDir  campaign output to overlay as MC TRUTH points. Empty disables the overlay.
///                The overlay is the point of the figure: analytic curves drawn alone are a
///                restatement of the formula, and prove nothing about what the simulation did.
///                Drawn together, a curve that misses its own generated events is visible at once.
void dp_kinematics_C17(Double_t Ebeam = 135.0, Double_t mBeamAmu = 17.02257865, Double_t mTgtAmu = 2.0141017778,
                       Double_t mResAmu = 18.02675193, Double_t mEjAmu = 1.0078250322,
                       TString outDir = "./plots/", TString simDir = "/mnt/f/C17dp/b285_attpc/")
{
   const double u = 931.49401;
   const double mB = mBeamAmu * u, mT = mTgtAmu * u, mR0 = mResAmu * u, mE = mEjAmu * u;

   // the four BOUND states of 18C (S_n = 4.184 MeV)
   const int nL = 4;
   const double Ex[nL] = {0.0, 1.588, 2.515, 3.972};
   const char *nm[nL] = {"0+ g.s.", "2+ 1.588", "(2+) 2.515", "(2,3)+ 3.972"};

   printf("\n\033[1;33m===== 17C(d,p)18C kinematics, Ebeam = %.2f MeV (%.3f MeV/u) =====\033[0m\n", Ebeam,
          Ebeam / 17);
   printf("  Q(g.s.) = %+.4f MeV     S_n(18C) = 4.184 MeV\n", mB + mT - mR0 - mE);

   const double Eb = Ebeam + mB;
   const double pb = std::sqrt(Eb * Eb - mB * mB);
   TLorentzVector Lb(0, 0, pb, Eb), Lt(0, 0, 0, mT);
   TLorentzVector W = Lb + Lt;
   TVector3 bst = W.BoostVector();
   const double s = W.M2();

   // lab (theta, KE) of the ejectile at a given cm angle and residual excitation.
   //
   // ANGLE CONVENTION -- GET THIS WRONG AND EVERY TABLE BELOW IS SILENTLY REVERSED.
   // theta_cm here is the STANDARD (d,p) angle, measured from the DEUTERON direction in the
   // entrance channel, which in inverse kinematics is the SUPPLEMENT of the proton's cm polar
   // angle about the beam. That is exactly what acceptance_C14.C uses -- it computes
   // theta_cm = pi - acos(...) (acceptance_C14.C:44) -- so this macro, the acceptance and the
   // resolution tables all mean the same thing by "theta_cm", and small theta_cm is where a
   // stripping angular distribution has its yield.
   // Verified against the 14C(d,p) reference campaign: at Ebeam 159.75 this reproduces its
   // documented "theta_cm 60 -> theta_lab 76 deg, 14.3 MeV" and "theta_cm 20 -> 125 deg, 3.1 MeV".
   auto labOf = [&](double thcmDeg, double ex, double &thLab, double &ke) {
      const double mR = mR0 + ex;
      const double arg = (s - (mR + mE) * (mR + mE)) * (s - (mR - mE) * (mR - mE));
      if (arg <= 0) {
         thLab = ke = -1;
         return false;
      }
      const double pcm = std::sqrt(arg) / (2 * std::sqrt(s));
      const double th = (180.0 - thcmDeg) * TMath::DegToRad();
      TLorentzVector L(pcm * std::sin(th), 0, pcm * std::cos(th), std::sqrt(pcm * pcm + mE * mE));
      L.Boost(bst);
      thLab = L.Vect().Theta() * TMath::RadToDeg();
      ke = L.E() - mE;
      return true;
   };

   for (int l = 0; l < nL; ++l) {
      printf("\n\033[1;36m--- 18C  %s  (Ex = %.3f MeV) ---\033[0m\n", nm[l], Ex[l]);
      printf("  %8s %10s %10s %12s %12s\n", "theta_cm", "theta_lab", "KE [MeV]", "dEx/dKE", "dEx/dtheta");
      printf("  %8s %10s %10s %12s %12s\n", "[deg]", "[deg]", "", "", "[MeV/deg]");
      for (double tcm : {2., 5., 8., 10., 15., 20., 30., 40., 60., 80., 100., 120., 150.}) {
         double th, ke;
         if (!labOf(tcm, Ex[l], th, ke))
            continue;
         // leverage: how much Ex moves when the MEASURED (KE, theta_lab) move, at fixed beam
         // energy. Computed by inverting the two-body relation numerically around this point.
         auto exOf = [&](double keM, double thM) {
            const double E = keM + mE;
            const double p = std::sqrt(std::max(0., E * E - mE * mE));
            const double t = thM * TMath::DegToRad();
            TLorentzVector Le(p * std::sin(t), 0, p * std::cos(t), E);
            return (W - Le).M() - mR0;
         };
         const double dKE = 0.01, dTh = 0.01;
         const double dExdKE = (exOf(ke + dKE, th) - exOf(ke - dKE, th)) / (2 * dKE);
         const double dExdTh = (exOf(ke, th + dTh) - exOf(ke, th - dTh)) / (2 * dTh);
         printf("  %8.0f %10.1f %10.2f %12.3f %12.4f\n", tcm, th, ke, dExdKE, dExdTh);
      }
   }

   // Level separation in the two observables, at the angles where a transfer distribution has its
   // yield. If two levels differ by less in KE than the detector resolves, no amount of statistics
   // separates them.
   printf("\n\033[1;33m===== level separation in proton KE (the measured quantity) =====\033[0m\n");
   printf("  %8s %10s", "theta_cm", "th_lab(gs)");
   for (int l = 0; l < nL; ++l)
      printf(" %11.3f", Ex[l]);
   printf("   [KE, MeV]\n");
   for (double tcm : {2., 5., 8., 10., 15., 20., 30., 40., 60., 80., 100., 120., 150.}) {
      double th0, ke0;
      labOf(tcm, 0.0, th0, ke0);
      printf("  %8.0f %10.1f", tcm, th0);
      for (int l = 0; l < nL; ++l) {
         double th, ke;
         printf(" %11.2f", labOf(tcm, Ex[l], th, ke) ? ke : -1);
      }
      printf("\n");
   }
   printf("\n  The gap between adjacent columns, divided by dEx/dKE, is the Ex separation the\n"
          "  detector has to beat. The tightest one in 18C is 1.588 -> 2.515, i.e. 927 keV.\n\n");

   if (outDir.IsNull())
      return;

   // ---------------------------------------------------------------------------------------------
   // THE FIGURE. Four panels, each framed on its own content -- a range chosen in advance hides
   // the shape of what actually came out.
   //   A  theta_lab vs theta_cm          where the protons go
   //   B  proton KE vs theta_lab         the locus a real spectrum sits on, with MC truth under it
   //   C  proton KE vs theta_cm          the level separation in the MEASURED quantity
   //   D  |dEx/dKE| vs theta_lab         the leverage that turns a KE error into an Ex error
   // The transfer-peak band (theta_cm 2-40) is shaded in A and C and marked in B, because every
   // number worth quoting comes from inside it and it is NOT where the flat generator puts most
   // of its events.
   gStyle->SetOptStat(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gSystem->mkdir(outDir, kTRUE);

   const int col[nL] = {kBlack, kRed + 1, kAzure + 2, kGreen + 3};
   const char *tag[nL] = {"gs", "ex1588", "ex2515", "ex3972"};
   const double cmPeakLo = 2.0, cmPeakHi = 40.0;

   TGraph *gThTh[nL], *gKeTh[nL], *gKeCm[nL], *gLev[nL];
   double keMax = 0, thLabAtPeakLo = 0, thLabAtPeakHi = 0;
   for (int l = 0; l < nL; ++l) {
      gThTh[l] = new TGraph();
      gKeTh[l] = new TGraph();
      gKeCm[l] = new TGraph();
      gLev[l] = new TGraph();
      for (double tcm = 1.0; tcm <= 179.0; tcm += 0.5) {
         double th, ke;
         if (!labOf(tcm, Ex[l], th, ke))
            continue;
         gThTh[l]->SetPoint(gThTh[l]->GetN(), tcm, th);
         gKeTh[l]->SetPoint(gKeTh[l]->GetN(), th, ke);
         gKeCm[l]->SetPoint(gKeCm[l]->GetN(), tcm, ke);
         keMax = std::max(keMax, ke);
         // leverage, same numeric inversion as the tables above
         auto exOf = [&](double keM, double thM) {
            const double E = keM + mE;
            const double p = std::sqrt(std::max(0., E * E - mE * mE));
            const double t = thM * TMath::DegToRad();
            TLorentzVector Le(p * std::sin(t), 0, p * std::cos(t), E);
            return (W - Le).M() - mR0;
         };
         gLev[l]->SetPoint(gLev[l]->GetN(), th, std::fabs((exOf(ke + 0.01, th) - exOf(ke - 0.01, th)) / 0.02));
      }
      if (l == 0) {
         double ke;
         labOf(cmPeakLo, Ex[0], thLabAtPeakLo, ke);
         labOf(cmPeakHi, Ex[0], thLabAtPeakHi, ke);
      }
   }

   TCanvas *ck = new TCanvas("ckin17", "17C(d,p)18C kinematics", 1400, 1000);
   ck->Divide(2, 2);

   auto peakBand = [&](double lo, double hi, double ymin, double ymax) {
      TBox *b = new TBox(lo, ymin, hi, ymax);
      b->SetFillColorAlpha(kOrange - 4, 0.20);
      b->SetLineWidth(0);
      b->Draw();
   };
   auto legend = [&](TGraph *g[nL], double x1, double y1, double x2, double y2) {
      TLegend *lg = new TLegend(x1, y1, x2, y2);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);
      for (int l = 0; l < nL; ++l)
         lg->AddEntry(g[l], nm[l], "l");
      lg->Draw();
      return lg;
   };

   // --- A: theta_lab vs theta_cm
   ck->cd(1);
   gPad->SetGrid();
   for (int l = 0; l < nL; ++l) {
      gThTh[l]->SetLineColor(col[l]);
      gThTh[l]->SetLineWidth(2);
      if (l == 0) {
         gThTh[l]->SetTitle(TString::Format("^{17}C(d,p)^{18}C at %.1f MeV;#theta_{cm} [deg];"
                                            "#theta_{lab} (proton) [deg]",
                                            Ebeam));
         gThTh[l]->GetYaxis()->SetRangeUser(0, 180);
         gThTh[l]->GetXaxis()->SetLimits(0, 180);
         gThTh[l]->Draw("AL");
         peakBand(cmPeakLo, cmPeakHi, 0, 180);
         gThTh[l]->Draw("L SAME");
      } else
         gThTh[l]->Draw("L SAME");
   }
   legend(gThTh, 0.50, 0.62, 0.88, 0.87);
   {
      TLatex *tx = new TLatex(0.16, 0.20, "#color[801]{transfer peak}");
      tx->SetNDC();
      tx->SetTextSize(0.038);
      tx->Draw();
   }

   // --- B: KE vs theta_lab, with MC TRUTH under the curves
   ck->cd(2);
   gPad->SetGrid();
   TH2D *frame = new TH2D("kinframe", "proton locus, curves = analytic, points = MC truth;"
                                      "#theta_{lab} (proton) [deg];proton KE [MeV]",
                          180, 0, 180, 200, 0, 1.08 * keMax);
   frame->Draw();
   long nOverlay = 0;
   if (!simDir.IsNull()) {
      for (int l = 0; l < nL; ++l) {
         TString found = gSystem->GetFromPipe(TString("ls -1 ") + simDir + "exres_" + tag[l] +
                                              "_s*_b285_attpc.root 2>/dev/null | head -1");
         found = found.Strip(TString::kBoth);
         if (found.IsNull())
            continue;
         TFile *f = TFile::Open(found);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (!t)
            continue;
         double thT, keT;
         t->SetBranchAddress("thTrue", &thT);
         t->SetBranchAddress("keTrue", &keT);
         TGraph *gp = new TGraph();
         for (Long64_t i = 0; i < t->GetEntries(); i += 3) { // thinned, it is only a visual check
            t->GetEntry(i);
            gp->SetPoint(gp->GetN(), thT, keT);
         }
         gp->SetMarkerStyle(1);
         gp->SetMarkerColorAlpha(col[l], 0.30);
         gp->Draw("P SAME");
         nOverlay += gp->GetN();
      }
   }
   for (int l = 0; l < nL; ++l) {
      gKeTh[l]->SetLineColor(col[l]);
      gKeTh[l]->SetLineWidth(3);
      gKeTh[l]->Draw("L SAME");
   }
   {
      // the transfer peak in LAB angle -- it runs backwards, hence hi/lo swapped
      TLine *l1 = new TLine(thLabAtPeakHi, 0, thLabAtPeakHi, 1.08 * keMax);
      TLine *l2 = new TLine(thLabAtPeakLo, 0, thLabAtPeakLo, 1.08 * keMax);
      for (TLine *ln : {l1, l2}) {
         ln->SetLineColor(kOrange + 7);
         ln->SetLineStyle(2);
         ln->SetLineWidth(2);
         ln->Draw();
      }
      // bottom-left is the empty corner here (forward lab angles carry the FAST protons), so the
      // label goes there rather than beside the legend
      TLatex *tx = new TLatex(0.15, 0.20,
                              TString::Format("#color[801]{transfer peak: #theta_{lab} %.0f-%.0f#circ}",
                                              thLabAtPeakHi, thLabAtPeakLo));
      tx->SetNDC();
      tx->SetTextSize(0.036);
      tx->Draw();
   }
   legend(gKeTh, 0.60, 0.55, 0.90, 0.78);

   // --- C: KE vs theta_cm
   ck->cd(3);
   gPad->SetGrid();
   for (int l = 0; l < nL; ++l) {
      gKeCm[l]->SetLineColor(col[l]);
      gKeCm[l]->SetLineWidth(2);
      if (l == 0) {
         gKeCm[l]->SetTitle("level separation in the measured quantity;#theta_{cm} [deg];proton KE [MeV]");
         gKeCm[l]->GetYaxis()->SetRangeUser(0, 1.08 * keMax);
         gKeCm[l]->GetXaxis()->SetLimits(0, 180);
         gKeCm[l]->Draw("AL");
         peakBand(cmPeakLo, cmPeakHi, 0, 1.08 * keMax);
         gKeCm[l]->Draw("L SAME");
      } else
         gKeCm[l]->Draw("L SAME");
   }
   legend(gKeCm, 0.55, 0.62, 0.88, 0.87);

   // --- D: leverage
   ck->cd(4);
   gPad->SetGrid();
   gPad->SetLogy();
   for (int l = 0; l < nL; ++l) {
      gLev[l]->SetLineColor(col[l]);
      gLev[l]->SetLineWidth(2);
      if (l == 0) {
         gLev[l]->SetTitle("leverage: a KE error becomes this much E_{x} error;"
                           "#theta_{lab} (proton) [deg];|dE_{x}/dKE|");
         gLev[l]->GetXaxis()->SetLimits(0, 180);
         gLev[l]->Draw("AL");
      } else
         gLev[l]->Draw("L SAME");
   }
   legend(gLev, 0.18, 0.62, 0.52, 0.87);

   ck->SaveAs(outDir + "kinematics_C17dp.png");
   printf("wrote %skinematics_C17dp.png", outDir.Data());
   if (nOverlay)
      printf("   (%ld MC truth protons overlaid in panel B)", nOverlay);
   printf("\n\n");
}
