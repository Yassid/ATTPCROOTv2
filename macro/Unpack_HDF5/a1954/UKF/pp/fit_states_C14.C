/// @file fit_states_C14.C
/// @brief Fit the a1954 14C(p,p') excitation spectrum with the known level scheme.
///
/// Level scheme and identifications follow Ayyad et al., "Direct reactions with the AT-TPC",
/// Front. Phys. 13:1539148 (2025), section 3.1, which analysed this same experiment:
///     6.091  1-      first excited state, well isolated, the one with a published ang. dist.
///     6.728  3-_1  \
///     7.012  2+_1   |  the "group of states around 7 MeV", per Lozowski
///     7.341  2-_1  /
///     8.317  2+_1    resolved in the paper
///     ~9.1           NOT 14C: an excited state in 14N above the proton emission threshold,
///                    populated by the (p,n) charge-exchange channel. Fitted with a FREE
///                    centroid and flagged, because forcing it onto a 14C level would silently
///                    absorb charge-exchange yield into a 14C state.
///
/// Strategy. The g.s. is fitted first, alone: it carries ~6e4 counts and sets both the energy
/// zero-point and the resolution. Its centroid offset and sigma are then CARRIED INTO the
/// excited-state fit, where the level energies are held at their literature values shifted by
/// that offset and only the amplitudes float. This is deliberate: with sigma ~0.17 and levels
/// 0.28-0.33 MeV apart, letting the centroids float lets neighbouring peaks trade yield and the
/// fit stops meaning anything. sigmaScale > 1 allows the excited peaks to be broader than the
/// g.s. (kinematic broadening grows with Ex), and is reported so it is never silent.
///
/// The paper quotes 150 keV (sigma) for the g.s. resolution with 30 keV accuracy; this fit
/// should reproduce that, and it is a check on the whole chain if it does not.
///
///   root -b -q 'fit_states_C14.C()'
///   root -b -q 'fit_states_C14.C("plots/proton_kin_300_ukf.root",5.0,5.6,10.2)'

void fit_states_C14(TString cache = "plots/proton_kin_300_ukf.root", Double_t chi2Cut = 5.0, Double_t exLo = 5.6,
                    Double_t exHi = 10.2, Double_t sigmaScale = -1, Int_t nbins = 260, Double_t axLo = -2,
                    Double_t axHi = 24)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *f = TFile::Open(here + "/" + cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) {
      printf("\033[1;31mno ntuple\033[0m\n");
      return;
   }

   auto *h = new TH1D("hEx", "14C(p,p') excitation energy;E_{x} [MeV];counts", nbins, axLo, axHi);
   h->Sumw2();
   t->Draw("ex>>hEx", TString::Format("chi2ndf<%g", chi2Cut), "goff");

   // ---- 1. ground state alone: sets the zero-point and the resolution -------------------
   TF1 gs("gs", "gaus(0)+pol1(3)", -1.2, 1.2);
   gs.SetParameters(h->GetMaximum(), 0, 0.17, 0, 0);
   gs.SetParNames("A_gs", "mu_gs", "sigma_gs", "bg0", "bg1");
   h->Fit(&gs, "QRN");
   const double off = gs.GetParameter(1);   // energy zero-point offset
   const double sig = std::fabs(gs.GetParameter(2));
   printf("\n===== ground state =====\n");
   printf("  centroid %+.4f MeV   sigma %.4f MeV (FWHM %.4f)   area %.0f\n", off, sig, 2.355 * sig,
          gs.GetParameter(0) * sig * std::sqrt(2 * TMath::Pi()) / h->GetBinWidth(1));
   printf("  paper (Front.Phys.13:1539148) quotes sigma = 0.150 MeV, accuracy 30 keV\n");

   // ---- 2. excited states: literature energies, shifted by the g.s. offset --------------
   struct Lvl { double e; const char *jp; bool c14; };
   // 6.589 (0+) and 6.903 (0-) are NuDat levels inside the same window. The paper names only
   // 6.728/7.012/7.341 for the "group around 7 MeV", but omitting the other two leaves the fit
   // with no way to describe the yield between them and drives chi2/ndf to ~5. Included here,
   // and their yields reported separately so the paper's three-level grouping stays recoverable.
   std::vector<Lvl> lv = {{6.091, "1-", true},   {6.589, "0+_2", true}, {6.728, "3-_1", true},
                          {6.903, "0-_1", true}, {7.012, "2+_1", true}, {7.341, "2-_1", true},
                          {8.317, "2+_1", true}, {9.100, "14N (p,n)", false}};
   const int nL = lv.size();

   // ONE SHARED WIDTH, FREE. Fixing sigma made the yields depend entirely on that choice:
   // scaling it 1.0 -> 1.15 -> 1.3 moved chi2/ndf 5.09 -> 4.13 -> 3.43 and the 1- yield
   // 280 -> 325 -> 385, a 30 % swing with no way to choose between them. Letting the data set
   // one common width removes that arbitrariness. Level energies are embedded as LITERALS
   // (fixed at literature + the measured zero-point) so only amplitudes, the shared width, the
   // 14N centroid and the background float -- neighbouring peaks cannot trade yield by sliding.
   TString form;
   int nPar = 0;
   std::vector<int> iAmp(nL);
   const int iSig = 0;      // [0] = shared sigma
   nPar = 1;
   for (int i = 0; i < nL; ++i) {
      iAmp[i] = nPar++;
      if (lv[i].c14)
         form += TString::Format("%s[%d]*exp(-0.5*((x-%.6f)/[%d])^2)", i ? "+" : "", iAmp[i], lv[i].e + off, iSig);
      else {
         int iMu = nPar++;
         form += TString::Format("%s[%d]*exp(-0.5*((x-[%d])/[%d])^2)", i ? "+" : "", iAmp[i], iMu, iSig);
         form += ""; // 14N centroid floats: it is a different nucleus, not a 14C level
      }
   }
   form += TString::Format("+[%d]+[%d]*x", nPar, nPar + 1);
   TF1 fx("fx", form, exLo, exHi);
   fx.SetParameter(iSig, sig * (sigmaScale > 0 ? sigmaScale : 1.0));
   fx.SetParLimits(iSig, 0.5 * sig, 4.0 * sig);
   if (sigmaScale > 0)
      fx.FixParameter(iSig, sig * sigmaScale); // opt-in: reproduce the old fixed-width behaviour
   fx.SetParName(iSig, "sigma");
   int p = 1;
   for (int i = 0; i < nL; ++i) {
      int b = h->FindBin(lv[i].e + off);
      fx.SetParameter(p, std::max(1.0, h->GetBinContent(b) * 0.5));
      fx.SetParLimits(p, 0, 1e7);
      fx.SetParName(p, TString::Format("A_%s", lv[i].jp));
      ++p;
      if (!lv[i].c14) {
         fx.SetParameter(p, lv[i].e + off);
         fx.SetParLimits(p, 8.7, 9.6);
         fx.SetParName(p, "mu_14N");
         ++p;
      }
   }
   fx.SetParameter(nPar, 10);
   fx.SetParameter(nPar + 1, 0);
   h->Fit(&fx, "QRN");
   const double sigFit = std::fabs(fx.GetParameter(iSig));

   const double bw = h->GetBinWidth(1), k = std::sqrt(2 * TMath::Pi()) / bw;
   printf("\n===== excited states (E fixed at literature %+.3f MeV; shared sigma %s = %.4f) =====\n", off,
          sigmaScale > 0 ? "FIXED" : "FITTED", sigFit);
   printf("  %-12s %8s   %10s %10s\n", "state", "E_x [MeV]", "yield", "err");
   double tot = 0;
   int q = 1;
   for (int i = 0; i < nL; ++i) {
      double A = fx.GetParameter(q), dA = fx.GetParError(q);
      ++q;
      double mu = lv[i].e + off;
      if (!lv[i].c14) { mu = fx.GetParameter(q); ++q; }
      double Y = A * sigFit * k, dY = dA * sigFit * k;
      tot += Y;
      printf("  %-12s %8.3f   %10.0f %10.0f%s\n", lv[i].jp, mu, Y, dY, lv[i].c14 ? "" : "   <-- NOT 14C");
   }
   // GROUPED totals. The 6.5-7.4 levels sit 0.14-0.28 MeV apart against a width of ~0.16-0.26,
   // so they are NOT resolvable and the fit trades yield freely between them: individual
   // amplitudes there swing from 0 to 2000 depending on the width, while their SUM is stable.
   // Report the sum, which is what the resolution supports, and never quote the split.
   {
      double y1 = 0, ygrp = 0, y83 = 0, yN = 0;
      int qq = 1;
      for (int i = 0; i < nL; ++i) {
         double A = fx.GetParameter(qq); ++qq;
         if (!lv[i].c14) ++qq;
         double Y = A * sigFit * k;
         if (!lv[i].c14) yN += Y;
         else if (lv[i].e < 6.3) y1 += Y;
         else if (lv[i].e < 7.6) ygrp += Y;
         else y83 += Y;
      }
      printf("  ---- GROUPED (what the resolution supports) ----\n");
      printf("    6.091 1-        %8.0f\n", y1);
      printf("    6.5-7.4 group   %8.0f   (do NOT quote the split)\n", ygrp);
      printf("    8.317 2+        %8.0f\n", y83);
      printf("    ~9.1 14N (p,n)  %8.0f\n", yN);
   }
   printf("  chi2/ndf = %.2f / %d = %.2f\n", fx.GetChisquare(), fx.GetNDF(),
          fx.GetNDF() > 0 ? fx.GetChisquare() / fx.GetNDF() : 0);
   printf("  summed peak yield %.0f\n", tot);

   // ---- draw ----------------------------------------------------------------------------
   TCanvas *c = new TCanvas("c", "states", 1400, 560);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogy();
   h->GetXaxis()->SetRangeUser(-2, 12);
   h->SetLineColor(kBlack);
   h->Draw("hist");
   gs.SetLineColor(kAzure + 2);
   gs.SetRange(-1.2, 1.2);
   gs.Draw("same");
   fx.SetLineColor(kRed + 1);
   fx.Draw("same");
   c->cd(2);
   auto *hz = (TH1D *)h->Clone("hz");
   hz->GetXaxis()->SetRangeUser(exLo - 0.6, exHi + 0.4);
   hz->SetTitle("excited states;E_{x} [MeV];counts");
   hz->Draw("hist");
   fx.Draw("same");
   // individual components
   int r = 1;
   for (int i = 0; i < nL; ++i) {
      double A = fx.GetParameter(r); ++r;
      double mu = lv[i].e + off;
      if (!lv[i].c14) { mu = fx.GetParameter(r); ++r; }
      auto *g1 = new TF1(TString::Format("c%d", i), "gaus", exLo, exHi);
      g1->SetParameters(A, mu, sigFit);
      g1->SetLineColor(lv[i].c14 ? kAzure + 2 : kGreen + 2);
      g1->SetLineStyle(2);
      g1->SetLineWidth(1);
      g1->Draw("same");
   }
   TString png = here + "/plots/states_C14.png";
   c->SaveAs(png);
   printf("wrote %s\n\n", png.Data());
}
