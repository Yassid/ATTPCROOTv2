/// @file exc8_vs_fresco_C14.C
/// @brief Test the multipolarity of the two structures above the multiplet against DWBA.
///
/// The two peaks at 8.53 and 9.36 MeV were separated using their ENERGIES alone, so their angular
/// SHAPES carry independent information about the transferred angular momentum. This macro
/// compares each to collective (type=11) DWBA for L = 1, 2, 3, 4 and ranks them by chi2.
///
/// A NOTE ON THE FRESCO INPUTS, because it is an easy and silent mistake. In a type=11 coupling
/// the five numbers p(1:5) are the deformation lengths PER MULTIPOLE -- p(1) is delta_1, p(2) is
/// delta_2, and so on. Changing only jt and ptyt on the &STATES line while leaving delta in the
/// p(2) slot leaves the state uncoupled, and FRESCO then returns a cross section of exactly zero
/// without any error. The symptom is several multipoles reporting an IDENTICAL chi2 (the
/// data-against-zero value), which is what happened on the first pass here. L = 0 is dropped
/// altogether: there is no p(0) slot, monopole is not a type=11 collective coupling.
///
/// The published analysis of these data assigns the lower structure to the 8.317 MeV 2+, i.e.
/// L = 2. That is a prediction this comparison can confirm or contradict; it is NOT used as an
/// input. The calculations are run at the LITERATURE energies 8.317 and 9.363 MeV, since the
/// exit-channel momentum depends on the excitation energy and the shape is sensitive to it here.
///
/// NEAR-THRESHOLD CAVEAT, and it is a serious one. E_cm is only 10.80 MeV for a 161 MeV beam, so
/// the exit channel retains 2.5 MeV at 8.317 and 1.4 MeV at 9.363. The outgoing wave number is
/// therefore small, the diffraction pattern is stretched to large angles, and successive L differ
/// less than they do for the 6-7 MeV multiplet. The discrimination between multipoles is
/// correspondingly weaker, and for the upper structure it is weak enough that no assignment
/// should be quoted from the shape alone.
///
/// Each curve carries ONE free normalisation, fitted over the angular range where the data are
/// significant. The comparison is of shape only.
///
/// NO RESOLUTION FOLDING is applied, following the choice made for the multiplet: the correction
/// proved smaller than the statistical errors and folding obscured the raw comparison.
///
///   root -b -q 'exc8_vs_fresco_C14.C()'

void exc8_vs_fresco_C14(TString dataFile = "plots/exc8_angdist_hi.root", TString frDir = "",
                        Double_t normLo = 30.0, Double_t normHi = 130.0, TString outTag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frDir.IsNull())
      frDir = here + "/../fresco/outputs/";

   const int NP = 2;
   const char *PNAME[NP] = {"8.53 MeV  (8.317, 2^{+}?)", "9.36 MeV  (blend?)"};
   const char *ETAG[NP] = {"8317", "9363"};
   const int LMIN = 1, NL = 4; // L = 1..4
   const int LCOL[NL] = {kAzure + 2, kRed + 1, kGreen + 3, kMagenta + 2};

   TFile *fd = TFile::Open(here + "/" + dataFile);
   if (!fd || fd->IsZombie()) {
      printf("\033[1;31mmissing %s -- run exc8_angdist_C14.C first\033[0m\n", dataFile.Data());
      return;
   }
   TH1D *dat[NP];
   for (int i = 0; i < NP; ++i) {
      auto *h = (TH1D *)fd->Get(TString::Format("dsdo_%d", i));
      if (!h) {
         printf("\033[1;31mno dsdo_%d in %s\033[0m\n", i, dataFile.Data());
         return;
      }
      dat[i] = (TH1D *)h->Clone(TString::Format("D%d", i));
      dat[i]->SetDirectory(nullptr);
   }
   fd->Close();

   auto loadFR = [&](TString fn) -> TGraph * {
      auto *g = new TGraph();
      std::ifstream in((frDir + fn).Data());
      double a, x;
      int n = 0;
      while (in >> a >> x)
         g->SetPoint(n++, a, x);
      if (g->GetN() == 0) {
         printf("\033[1;31mcannot read %s\033[0m\n", fn.Data());
         return nullptr;
      }
      return g;
   };

   TCanvas *c = new TCanvas("c8f", "upper states vs DWBA", 1400, 600);
   c->Divide(2, 1);

   for (int i = 0; i < NP; ++i) {
      c->cd(i + 1);
      gPad->SetLogy();

      // which bins carry information at all
      std::vector<int> use;
      for (int b = 1; b <= dat[i]->GetNbinsX(); ++b) {
         double x = dat[i]->GetBinCenter(b), y = dat[i]->GetBinContent(b), e = dat[i]->GetBinError(b);
         if (x >= normLo && x <= normHi && y > 0 && e > 0 && y > 2 * e)
            use.push_back(b);
      }
      printf("\n===== %s : %zu usable bins =====\n", PNAME[i], use.size());
      if (use.size() < 3) {
         printf("  too few significant bins to discriminate between multipoles -- no ranking\n");
      }

      double mx = dat[i]->GetMaximum();
      auto *frm = new TH1D(TString::Format("F%d", i),
                           TString::Format("%s;#theta_{cm} [deg];d#sigma/d#Omega [arb]", PNAME[i]), 1, 15, 145);
      frm->SetMinimum(std::max(1.0, mx * 5e-3));
      frm->SetMaximum(mx * 30);
      frm->Draw();

      auto *lg = new TLegend(0.14, 0.66, 0.55, 0.89);
      lg->SetNColumns(2);
      lg->SetTextSize(0.033);

      int best = -1;
      double bestC2 = 1e30;
      double c2v[NL];
      for (int j = 0; j < NL; ++j) {
         const int L = LMIN + j;
         TGraph *fr = loadFR(TString::Format("p14C_inel_161_%s_L%d_dsdo_ex2.dat", ETAG[i], L));
         c2v[j] = -1;
         if (!fr)
            continue;
         // one free normalisation, weighted least squares on the usable bins
         double sn = 0, sd = 0;
         for (int b : use) {
            double y = dat[i]->GetBinContent(b), e = dat[i]->GetBinError(b), f = fr->Eval(dat[i]->GetBinCenter(b));
            if (f <= 0)
               continue;
            sn += y * f / (e * e);
            sd += f * f / (e * e);
         }
         double k = sd > 0 ? sn / sd : 0;
         double c2 = 0;
         int nd = 0;
         for (int b : use) {
            double y = dat[i]->GetBinContent(b), e = dat[i]->GetBinError(b),
                   f = fr->Eval(dat[i]->GetBinCenter(b)) * k;
            c2 += (y - f) * (y - f) / (e * e);
            ++nd;
         }
         c2v[j] = nd > 1 ? c2 / (nd - 1) : -1; // one parameter fitted: the normalisation
         if (c2v[j] >= 0 && c2v[j] < bestC2) {
            bestC2 = c2v[j];
            best = j;
         }
         auto *g2 = new TGraph(*fr);
         for (int n = 0; n < g2->GetN(); ++n)
            g2->SetPointY(n, g2->GetPointY(n) * k);
         g2->SetLineColor(LCOL[j]);
         g2->SetLineWidth(L == 2 ? 4 : 2); // L=2 is the published assignment for the lower peak
         g2->SetLineStyle(1);
         g2->Draw("L same");
         lg->AddEntry(g2, TString::Format("L=%d  #chi^{2}/n %.1f", L, c2v[j]), "l");
      }

      dat[i]->SetMarkerStyle(20);
      dat[i]->SetMarkerColor(kBlack);
      dat[i]->SetLineColor(kBlack);
      dat[i]->SetLineWidth(2);
      dat[i]->SetMarkerSize(1.3);
      dat[i]->Draw("E1 same");
      lg->AddEntry(dat[i], "data", "lp");
      lg->Draw();

      printf("      L :  chi2/ndf   (one free normalisation)\n");
      for (int j = 0; j < NL; ++j)
         printf("      %d :  %8.2f %s\n", LMIN + j, c2v[j], j == best ? "  <-- best" : "");
      if (best >= 0) {
         // is the best actually distinguishable from the others?
         double second = 1e30;
         for (int j = 0; j < NL; ++j)
            if (j != best && c2v[j] >= 0)
               second = std::min(second, c2v[j]);
         double sep = second - bestC2;
         printf("\n      best L = %d at chi2/ndf %.2f; next best is %.2f, a separation of %.2f.\n", LMIN + best,
                bestC2, second, sep);
         // Two separate questions, and both must be answered before quoting an assignment:
         // does ANY multipole describe the data, and is the best one distinguishable from the rest?
         const double kBadFit = 3.0;
         if (bestC2 > kBadFit) {
            printf("      \033[1;31mNO MULTIPOLE DESCRIBES THIS DISTRIBUTION.\033[0m Even the best is at chi2/ndf\n"
                   "      %.1f, so ranking them is meaningless -- 'best' here is only the least bad. This is\n"
                   "      what a blend of levels with different multipolarities looks like, and it is\n"
                   "      consistent with the measured width being twice the resolution.\n",
                   bestC2);
         } else if (sep < 1.0 || use.size() < 4) {
            printf("      \033[1;33mTHAT IS NOT A DISCRIMINATION.\033[0m With %zu usable bins and a separation of\n"
                   "      only %.2f in chi2/ndf, these multipoles are not distinguishable in these data.\n"
                   "      No spin-parity assignment should be quoted from this shape.\n",
                   use.size(), sep);
         } else {
            printf("      The separation is large enough to prefer L = %d over the alternatives.\n", LMIN + best);
         }
      }
   }

   TString png = here + "/plots/exc8_vs_fresco_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
