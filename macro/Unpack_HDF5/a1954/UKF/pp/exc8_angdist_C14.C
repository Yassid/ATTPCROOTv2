/// @file exc8_angdist_C14.C
/// @brief Angular distributions of the two structures ABOVE the 6-7.5 MeV multiplet.
///
/// Two peaks sit over a clear minimum at 8.05-8.15 MeV: one at 8.533 MeV with sigma 0.162
/// (matching the instrumental resolution, so a single level -- identified in the published
/// analysis of these data as the 8.317 MeV 2+), and one at 9.363 MeV with sigma 0.299 (twice the
/// resolution, so probably a blend, and not discussed in that analysis).
///
/// Both are fitted together per angular bin, because the tail of the upper one reaches the lower
/// one and fitting either alone lets the background absorb the difference. Centroids and widths
/// are fixed at the values measured on the angle-integrated spectrum; only the two amplitudes and
/// a linear background float.
///
/// THE WIDTH OF THE UPPER PEAK IS FIXED AT ITS MEASURED VALUE, NOT AT THE RESOLUTION. For the
/// lower peak the two agree (0.162 measured, 0.162 from the resolution model), which is the
/// evidence that it is a single level. For the upper one they do not (0.299 against 0.172), so
/// fixing it at the resolution would force the fit to explain a blend with a single narrow
/// component. Its yield is therefore the yield of whatever group sits there, not of one level.
///
/// STATISTICS ARE THE LIMIT HERE: about 180 counts in the lower peak over the whole angular
/// range, against 3300 in the multiplet. 20 deg bins are used rather than 10, and even so the
/// lower peak carries only ~30 counts per bin.
///
/// THE Ex BINNING WAS SCANNED (50, 100, 150 keV) AND THE YIELDS ARE STABLE: the 8.53 peak gives
/// 0/16/39/73/50/3 counts per angular bin at 50 keV and 0/16/40/74/50/3 at 100 keV. Only the
/// reported chi2/ndf moves (1.92, 1.41, 1.47 angle-integrated), which is expected: the fit is a
/// likelihood fit, so at ~1 count per 50 keV bin the Pearson chi2 printed alongside it is itself
/// a noisy statistic and should not be read as a goodness-of-fit at that binning. 100 keV is
/// used, being the best-behaved of the three. 150 keV starts to lose the peak shape.
///
/// THE ACCEPTANCE IS NOW MEASURED AT THIS EXCITATION ENERGY. A dedicated simulation was run at
/// Ex = 8.317 MeV (three seeds, diagnostics/acc_level_genfit.sh) rather than reusing the 6.094
/// one. That matters more at the ends of the range than in the middle: against the 6.094
/// acceptance the yields move by -6.7 % at 30 deg and -9.1 % at 130 deg, but by only 1 % or so
/// between 50 and 110 deg. Extrapolating the gs-to-6.094 gradient would have predicted the
/// middle to 0.4 % and the ends to only 5-8 %, so the measurement was worth making.
///
/// The 9.363 MeV structure still uses it, being 1 MeV further up with no simulation of its own;
/// that remains the leading systematic on the upper distribution alone.
///
///   root -b -q 'exc8_angdist_C14.C()'

namespace
{
// METHOD CHOICE, and the two alternatives that were tried and rejected. Both peaks are fitted
// together over 7.9-10.3 with a linear background. The alternatives:
//
//   * extend the window down to 7.45 and add a component at 7.36 for the 7.27 multiplet member,
//     so that the multiplet tail leaking in below 8.2 is modelled rather than absorbed by the
//     background. WORSE: chi2/ndf 3.25 against 2.99 in the worst bin, 1.75 against 1.18
//     angle-integrated. The data below 7.6 MeV are the multiplet's entire right flank, far
//     broader than any single gaussian, so the added component absorbed multiplet strength
//     (941 counts) instead of describing a tail, and pulled the 8.53 yield down to 145 +- 18.
//
//   * fit each peak alone in a narrow window bounded by the valleys either side of it
//     (8.15-8.95 and 8.95-10.10). WORST: over so short a range a free linear background is
//     degenerate with the peak amplitude, and the yields came back meaningless -- 0 +- 57 and
//     4 +- 16 in bins where this method gives 0 +- 4 and 73 +- 10.
//
// This method and the first alternative agree on the angular SHAPE within about 1.5 sigma (both
// peak at 80-100 deg); their difference is quoted as the method systematic. The second is
// discarded rather than averaged in, because its failure mode is understood.
const int NP = 2;
const char *PNAME[NP] = {"8.53 (8.317?)", "9.36 (blend?)"};
const int PCOL[NP] = {kAzure + 2, kRed + 1};
double gMu[NP] = {8.533, 9.363};
double gSg[NP] = {0.162, 0.299};

const double FITLO = 7.9, FITHI = 10.3;

double twoPeak(double *x, double *p)
{
   double y = p[NP] + p[NP + 1] * x[0];
   for (int i = 0; i < NP; ++i)
      y += p[i] / (gSg[i] * std::sqrt(2 * TMath::Pi())) * std::exp(-0.5 * std::pow((x[0] - gMu[i]) / gSg[i], 2));
   return y;
}
} // namespace

void exc8_angdist_C14(TString cache = "plots/proton_kin_300gfx_ex.root",
                      TString accDir = "/mnt/f/a1954_C14_acc_gf_nochi2/", TString accLevel = "ex8",
                      Double_t zMin = -1e9, Double_t zMax = 1e9,
                      Double_t cmMin = 20.0,
                      Double_t cmMax = 140.0, Double_t dcm = 20.0, Int_t minN = 25, Int_t nEx = 27, // 100 keV bins; see the binning note above
                      
                      TString tag = "hi")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + cache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_" + accLevel + ".root");
   if (!fd || fd->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mcannot open the cache or the acceptance\033[0m\n");
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_" + accLevel + "_sum");
   if (!t || !acc)
      return;

   // ---- angle-integrated first: confirm the fixed parameters on this cache ----
   auto *hAll = new TH1D("hAll8", "", nEx + 3, 7.7, 10.4);
   // see exc_angdist_C14.C: the z window exists so the elastic normalisation cancels the flux,
   // the target thickness and the trigger efficiency
   TString zcut = (zMin > -1e8 || zMax < 1e8) ? TString::Format("&&vertexz>%g&&vertexz<%g", zMin, zMax) : TString("");
   t->Draw("ex>>hAll8", TString::Format("thcm>=%g&&thcm<%g", cmMin, cmMax) + zcut, "goff");
   hAll->SetDirectory(nullptr);
   double bwA = hAll->GetBinWidth(1);
   TF1 fA("fA", twoPeak, FITLO, FITHI, NP + 2);
   for (int i = 0; i < NP; ++i) {
      fA.SetParameter(i, 0.2 * hAll->Integral() * bwA);
      fA.SetParLimits(i, 0, 20 * hAll->Integral() * bwA);
   }
   fA.SetParameter(NP, 0.05 * hAll->Integral());
   fA.SetParameter(NP + 1, 0);
   hAll->Fit(&fA, "QNR");
   printf("\n===== angle-integrated, %.0f counts in 7.7-10.4 MeV, chi2/ndf %.2f =====\n", hAll->Integral(),
          fA.GetNDF() > 0 ? fA.GetChisquare() / fA.GetNDF() : -1);
   for (int i = 0; i < NP; ++i)
      printf("  %-14s mu %6.3f (fixed)  sigma %5.3f (fixed)  area %6.0f +- %.0f\n", PNAME[i], gMu[i], gSg[i],
             fA.GetParameter(i) / bwA, fA.GetParError(i) / bwA);

   // ---- per angle bin ----
   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   auto mk = [&](const char *n) {
      auto *h = new TH1D(n, "", NB, cmMin, cmMax);
      h->Sumw2();
      h->SetDirectory(nullptr);
      return h;
   };
   TH1D *yld[NP], *dsd[NP];
   for (int i = 0; i < NP; ++i) {
      yld[i] = mk(TString::Format("y8_%d", i));
      dsd[i] = mk(TString::Format("d8_%d", i));
   }

   printf("\n  theta_cm |  N  |");
   for (int i = 0; i < NP; ++i)
      printf(" %17s |", PNAME[i]);
   printf(" chi2/ndf\n");
   TCanvas *cf = new TCanvas("cf8", "fits", 1600, 900);
   cf->Divide(3, 2);
   for (int b = 1; b <= NB; ++b) {
      double lo = cmMin + (b - 1) * dcm, hi = lo + dcm;
      auto *h = new TH1D(TString::Format("h8_%d", b), "", nEx, 7.7, 10.4);
      t->Draw(TString::Format("ex>>h8_%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, hi) + zcut, "goff");
      h->SetDirectory(nullptr);
      double N = h->Integral();
      if (N < minN) {
         printf("  %3.0f-%3.0f  | %4.0f | only %.0f counts -- DROPPED\n", lo, hi, N, N);
         delete h;
         continue;
      }
      double bw = h->GetBinWidth(1);
      TF1 fn(TString::Format("f8_%d", b), twoPeak, FITLO, FITHI, NP + 2);
      for (int i = 0; i < NP; ++i) {
         fn.SetParameter(i, 0.2 * N * bw);
         fn.SetParLimits(i, 0, 20 * N * bw);
      }
      fn.SetParameter(NP, 0.05 * N);
      fn.SetParameter(NP + 1, 0);
      int st = h->Fit(&fn, "QNRL");
      double c2n = fn.GetNDF() > 0 ? fn.GetChisquare() / fn.GetNDF() : -1;
      printf("  %3.0f-%3.0f  | %4.0f |", lo, hi, N);
      for (int i = 0; i < NP; ++i) {
         double a = fn.GetParameter(i) / bw, e = fn.GetParError(i) / bw;
         printf(" %8.0f +-%5.0f |", a, e);
         // a yield consistent with zero is reported, not silently stored as a measurement
         if (st == 0 && a > 2 * e) {
            yld[i]->SetBinContent(b, a);
            yld[i]->SetBinError(b, e);
         }
      }
      printf(" %6.2f\n", c2n);
      if (b <= 6) {
         cf->cd(b);
         h->SetTitle(TString::Format("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", lo, hi));
         h->SetLineColor(kBlack);
         h->Draw("hist");
         auto *fc = (TF1 *)fn.DrawCopy("same");
         fc->SetLineColor(kRed + 1);
         for (int i = 0; i < NP; ++i) {
            auto *g = new TF1(TString::Format("p8_%d_%d", b, i), "gaus", FITLO, FITHI);
            g->SetParameters(fn.GetParameter(i) / (gSg[i] * std::sqrt(2 * TMath::Pi())), gMu[i], gSg[i]);
            g->SetLineColor(PCOL[i]);
            g->SetLineStyle(2);
            g->SetLineWidth(2);
            g->Draw("same");
         }
      } else
         delete h;
   }
   cf->SaveAs(here + "/plots/exc8_fits_" + tag + ".png");

   for (int i = 0; i < NP; ++i)
      for (int b = 1; b <= NB; ++b) {
         double c = yld[i]->GetBinCenter(b), y = yld[i]->GetBinContent(b);
         if (y <= 0)
            continue;
         double A = acc->GetBinContent(acc->FindBin(c)), s = std::sin(c * TMath::DegToRad());
         if (A <= 0.05 || s <= 1e-3)
            continue;
         dsd[i]->SetBinContent(b, y / A / s);
         dsd[i]->SetBinError(b, yld[i]->GetBinError(b) / A / s);
      }

   printf("\n  summed over the bins that are significant:");
   for (int i = 0; i < NP; ++i)
      printf("  %s %.0f", PNAME[i], yld[i]->Integral());
   printf("\n");

   TCanvas *c1 = new TCanvas("c8", "upper states", 1300, 560);
   c1->Divide(2, 1);
   auto style = [&](TH1D *h, int i) {
      h->SetMarkerStyle(20 + i);
      h->SetMarkerColor(PCOL[i]);
      h->SetLineColor(PCOL[i]);
      h->SetLineWidth(2);
      h->SetMarkerSize(1.3);
   };
   c1->cd(1);
   double mx = 0;
   for (int i = 0; i < NP; ++i)
      mx = std::max(mx, yld[i]->GetMaximum());
   auto *fr1 = new TH1D("fr8", "raw yield;#theta_{cm} [deg];counts / bin", 1, cmMin, cmMax);
   fr1->SetMinimum(0);
   fr1->SetMaximum(mx * 1.6);
   fr1->Draw();
   auto *lg = new TLegend(0.55, 0.72, 0.89, 0.88);
   for (int i = 0; i < NP; ++i) {
      style(yld[i], i);
      yld[i]->Draw("E1 same");
      lg->AddEntry(yld[i], PNAME[i], "lp");
   }
   lg->Draw();
   c1->cd(2);
   mx = 0;
   for (int i = 0; i < NP; ++i)
      mx = std::max(mx, dsd[i]->GetMaximum());
   auto *fr2 = new TH1D("fr82", "acceptance-corrected d#sigma/d#Omega (shape);#theta_{cm} [deg];arb.", 1, cmMin,
                        cmMax);
   fr2->SetMinimum(0);
   fr2->SetMaximum(mx * 1.6);
   fr2->Draw();
   for (int i = 0; i < NP; ++i) {
      style(dsd[i], i);
      dsd[i]->Draw("E1 same");
   }
   lg->Draw();
   TString png = here + "/plots/exc8_angdist_" + tag + ".png";
   c1->SaveAs(png);

   TFile fo(here + "/plots/exc8_angdist_" + tag + ".root", "RECREATE");
   for (int i = 0; i < NP; ++i) {
      yld[i]->Write(TString::Format("yield_%d", i));
      dsd[i]->Write(TString::Format("dsdo_%d", i));
   }
   hAll->Write("ex_hi");
   fo.Close();
   printf("\nwrote %s\n       plots/exc8_fits_%s.png\n\n", png.Data(), tag.Data());
}
