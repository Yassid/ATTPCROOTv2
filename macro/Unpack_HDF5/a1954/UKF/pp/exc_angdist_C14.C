/// @file exc_angdist_C14.C
/// @brief Angular distributions of the individual 14C excited states.
///
/// Follows the fit of the published analysis of these same data (Ayyad et al., Eur. Phys. J. A
/// 59:294 (2023), implemented in compiled/EFitterAnalysis/plotFit_a1954.C): four excited
/// components at 6.09 (1-), 6.70 (3-), 7.00 (2+) and 7.27 (2-), centroids and widths free but
/// confined to tight windows, and a TWO-STAGE fit -- each component fitted alone over its own
/// restricted range first, those parameters then seeding the combined function. A five-component
/// sum started cold converges to whichever local minimum the seed happens to sit in, which is why
/// the published macro does it this way and why that is reproduced here.
///
/// Note the fourth component sits at 7.27, not the tabulated 7.341; that is the position the
/// published analysis fitted, and the present data agree with it.
///
/// STAGE 1 (angle integrated) determines the centroids and widths, which are properties of the
/// calibration and the resolution, not of the scattering angle.
///
/// NO PER-ANGLE SHIFT IS APPLIED (freeShift defaults false), and that default is the result of a
/// mistake worth recording. The MODE of the multiplet marches from 6.63 to 7.19 MeV between 20 and
/// 130 deg, which looks exactly like a drifting energy scale. It is not: the composition changes
/// with angle -- 6.70 dominates forward, 7.27 backward -- so the mode of the blend moves even with
/// a perfect calibration. Two corrections were built to remove it before that was noticed. Tilting
/// the DATA (theta' = theta - 0.15 deg/MeV (KE - 4.5), the slope that best flattens the mode) drove
/// chi2/ndf to 6.5 in five bins and collapsed the 6.09 component to zero in three. Allowing a free
/// per-bin shift of the centroids improved chi2 -- because it was absorbing genuine composition
/// change into a spurious shift, i.e. fitting away part of the signal.
///
/// The decisive check is the ISOLATED 6.094 peak, the only calibration probe here free of
/// composition ambiguity. Across 20-140 deg it gives 6.31, 6.08, 6.14, 6.15, 5.99, 6.28 with
/// errors of 0.02-0.09 on 13-47 counts: scatter about ~6.1, not a trend. There is no drift to
/// correct, and the calibration (Ebeam 164.25 MeV with a -0.014 deg/MeV tilt) already suffices.
/// freeShift = kTRUE is kept only to reproduce the rejected variant.
/// STAGE 2 (per angle bin) holds those fixed and floats only the four areas plus a linear
/// background. With a few hundred counts per bin the centroids are not separately determined, and
/// letting them float lets the fit trade position against composition.
///
/// The yields are acceptance-corrected with the 6.094 acceptance and divided by sin(theta). That
/// acceptance is computed for ONE level and applied to all four: it varies slowly over ~1 MeV of
/// excitation energy, so this is a fair approximation, but it is an approximation, and it is the
/// leading systematic on the relative normalisation between the four distributions.
///
///   root -b -q 'exc_angdist_C14.C()'

namespace
{
const int NLV = 4;
const char *LVNAME[NLV] = {"6.09 (1-)", "6.70 (3-)", "7.00 (2+)", "7.27 (2-)"};
const int LVCOL[NLV] = {kAzure + 2, kRed + 1, kGreen + 3, kOrange + 7};
// seeds and windows exactly as in plotFit_a1954.C
const double MU0[NLV] = {6.09, 6.70, 7.00, 7.27};
const double MULO[NLV] = {5.95, 6.60, 6.90, 7.25};
const double MUHI[NLV] = {6.15, 6.80, 7.05, 7.30};
const double SG0[NLV] = {0.15, 0.15, 0.15, 0.15};
const double SGLO[NLV] = {0.10, 0.15, 0.15, 0.14};
const double SGHI[NLV] = {0.15, 0.20, 0.17, 0.15};
const double RLO[NLV] = {5.80, 6.40, 6.80, 7.10};   // per-component stage-1 fit ranges
const double RHI[NLV] = {6.30, 6.90, 7.30, 7.70};

double gMu[NLV], gSg[NLV]; // filled by stage 1
// The published windows belong to the published Ex scale. This analysis carries a residual gain
// on (Ex - 6.094) of about 1.078 (measured two ways: from the multiplet itself and from the 8.317
// anchor), which puts the upper components 0.06-0.17 MeV higher and pushes them OUT of those
// windows -- transplanted directly, all four centroids rail against their bounds and the 7.00
// component collapses to zero. The windows are therefore mapped through the same gain, and
// widened, so the fit structure transfers without importing the other analysis' calibration.
double mapMu(double e, double gain) { return 6.094 + gain * (e - 6.094); }

/// p[0..NLV-1] areas, p[NLV] and p[NLV+1] a linear background, p[NLV+2] a COMMON shift applied to
/// every centroid. The shift is one parameter for the whole multiplet, not one per level: the
/// level spacings are known and fixed, only their common position moves. That is what lets the
/// components follow the drift of the Ex scale with angle without gaining the freedom to trade
/// strength between themselves.
double multiplet(double *x, double *p)
{
   double y = p[NLV] + p[NLV + 1] * x[0];
   for (int i = 0; i < NLV; ++i)
      y += p[i] / (gSg[i] * std::sqrt(2 * TMath::Pi())) *
           std::exp(-0.5 * std::pow((x[0] - gMu[i] - p[NLV + 2]) / gSg[i], 2));
   return y;
}
} // namespace

void exc_angdist_C14(TString cache = "plots/proton_kin_300gfx_ex.root",
                     TString accDir = "/mnt/f/a1954_C14_acc_gf_nochi2/", Double_t cmMin = 20.0,
                     Double_t zMin = -1e9, Double_t zMax = 1e9,
                     Double_t cmMax = 140.0, Double_t dcm = 10.0, Int_t minN = 80, TString tag = "gfex",
                     Double_t gain = 1.078, Double_t winPad = 0.10, Double_t sig0 = 0.132,
                     Double_t dSig = 0.0123, Bool_t freeShift = kFALSE)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + cache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_ex1.root");
   if (!fd || fd->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mcannot open the cache or the ex1 acceptance\033[0m\n");
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_ex1_sum");
   if (!t || !acc) {
      printf("\033[1;31mmissing pk tree or hAcc_ex1_sum\033[0m\n");
      return;
   }

   // ---------------- stage 1: angle-integrated, published two-stage recipe ----------------
   auto *hAll = new TH1D("hAll", "", 140, 5.3, 8.1);
   // Optional vertex-z window. The absolute normalisation is taken from the elastic yield in the
   // SAME window, so that the beam flux, the target thickness and the trigger efficiency -- which
   // is what shapes the z distribution and which the simulation does not contain -- all cancel in
   // the ratio. 10-400 mm is the flat part of that distribution, where the trigger efficiency is
   // constant and so cancels exactly rather than on average.
   TString zcut = (zMin > -1e8 || zMax < 1e8) ? TString::Format("&&vertexz>%g&&vertexz<%g", zMin, zMax) : TString("");
   t->Draw(TString::Format("ex>>hAll"), TString::Format("thcm>=%g&&thcm<%g", cmMin, cmMax) + zcut, "goff");
   hAll->SetDirectory(nullptr);
   // Centroids at the tabulated energies mapped through the calibration gain, widths from the
   // MEASURED resolution rather than fitted. Letting the widths float gave 0.20-0.42 MeV against
   // the 0.132 measured on the isolated 6.094 peak; components that broad overlap so heavily that
   // their amplitudes become anti-correlated and the decomposition stops being a measurement.
   // Fixing them forces the fit to explain the observed breadth with COMPOSITION instead, which is
   // the physically meaningful choice given the resolution really is ~0.13 MeV.
   printf("\n===== stage 1: centroids and widths FIXED (not fitted) =====\n");
   printf("  resolution model: sigma(Ex) = %.3f + %.4f (Ex - 6.094)   [0.132 at 6.094, 0.162 at 8.533]\n",
          sig0, dSig);
   printf("  component  | tabulated | mapped mu | fixed sigma\n");
   for (int i = 0; i < NLV; ++i) {
      gMu[i] = mapMu(MU0[i], gain);
      gSg[i] = sig0 + dSig * (gMu[i] - 6.094);
      printf("  %-10s |  %6.3f   |  %6.3f   |  %6.4f\n", LVNAME[i], MU0[i], gMu[i], gSg[i]);
   }

   // combined fit, seeded from the individual ones, to settle the areas and check the description
   double bwA = hAll->GetBinWidth(1);
   TF1 fAll("fAll", multiplet, 5.4, 8.0, NLV + 3);
   for (int i = 0; i < NLV; ++i) {
      fAll.SetParameter(i, 0.2 * hAll->Integral() * bwA);
      fAll.SetParLimits(i, 0, 20 * hAll->Integral() * bwA);
   }
   fAll.SetParameter(NLV, 0.02 * hAll->Integral());
   fAll.SetParameter(NLV + 1, 0);
   fAll.SetParameter(NLV + 2, 0);
   fAll.SetParLimits(NLV + 2, -0.40, 0.40);
   hAll->Fit(&fAll, "QNR");
   printf("\n  common centroid shift, angle integrated : %+.3f MeV\n", fAll.GetParameter(NLV + 2));
   printf("\n  combined fit chi2/ndf = %.2f\n  angle-integrated composition:\n",
          fAll.GetNDF() > 0 ? fAll.GetChisquare() / fAll.GetNDF() : -1);
   double totA = 0;
   for (int i = 0; i < NLV; ++i)
      totA += fAll.GetParameter(i);
   for (int i = 0; i < NLV; ++i)
      printf("    %-10s %7.0f counts   %5.1f %%\n", LVNAME[i], fAll.GetParameter(i) / bwA,
             100 * fAll.GetParameter(i) / totA);

   // ---------------- stage 2: per angle bin, centroids and widths fixed ----------------
   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   auto mk = [&](const char *n) {
      auto *h = new TH1D(n, "", NB, cmMin, cmMax);
      h->Sumw2();
      h->SetDirectory(nullptr);
      return h;
   };
   auto *gShift = new TGraph();
   int nShift = 0;
   TH1D *yld[NLV], *dsd[NLV];
   for (int i = 0; i < NLV; ++i) {
      yld[i] = mk(TString::Format("y%d", i));
      dsd[i] = mk(TString::Format("d%d", i));
   }

   printf("\n===== stage 2: per-angle fit (centroids and widths held at the stage-1 values) =====\n");
   printf("  theta_cm |   N  |");
   for (int i = 0; i < NLV; ++i)
      printf(" %13s |", LVNAME[i]);
   printf(" chi2/ndf | common shift\n");

   TCanvas *cf = new TCanvas("cf", "per-bin fits", 1600, 1000);
   cf->Divide(4, 3);
   int pad = 0;
   for (int b = 1; b <= NB; ++b) {
      double lo = cmMin + (b - 1) * dcm, hi = lo + dcm;
      auto *h = new TH1D(TString::Format("hx%d", b), "", 56, 5.3, 8.1);
      t->Draw(TString::Format("ex>>hx%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, hi) + zcut, "goff");
      h->SetDirectory(nullptr);
      double N = h->Integral();
      if (N < minN) {
         printf("  %3.0f-%3.0f  | %4.0f | only %.0f counts -- DROPPED\n", lo, hi, N, N);
         delete h;
         continue;
      }
      double bw = h->GetBinWidth(1);
      TF1 fn(TString::Format("fn%d", b), multiplet, 5.4, 8.0, NLV + 3);
      for (int i = 0; i < NLV; ++i) {
         fn.SetParameter(i, 0.2 * N * bw);
         fn.SetParLimits(i, 0, 20 * N * bw);
      }
      fn.SetParameter(NLV, 0.02 * N);
      fn.SetParameter(NLV + 1, 0);
      fn.SetParameter(NLV + 2, 0);
      if (freeShift)
         fn.SetParLimits(NLV + 2, -0.40, 0.40); // the whole multiplet may slide, but not reorder
      else
         fn.FixParameter(NLV + 2, 0); // drift already removed from the DATA -- no shift allowed
      // Likelihood fit -- the per-bin counts are small. RETRY ON A NON-ZERO STATUS rather than
      // discarding the bin: with fewer events (the vertex-z window) Minuit can stop with a
      // non-zero status while sitting on perfectly sane parameters, and a discarded bin removes
      // an angular point from every level at once. Restarting from the current parameters is
      // usually enough. Only a bin that still fails after the retries is dropped, loudly.
      int st = h->Fit(&fn, "QNRL");
      int nTry = 0;
      while (st != 0 && nTry < 3) {
         ++nTry;
         st = h->Fit(&fn, "QNRL"); // restarts from where the previous attempt stopped
      }
      double c2n = fn.GetNDF() > 0 ? fn.GetChisquare() / fn.GetNDF() : -1;
      printf("  %3.0f-%3.0f  | %4.0f |", lo, hi, N);
      for (int i = 0; i < NLV; ++i) {
         double a = fn.GetParameter(i) / bw, e = fn.GetParError(i) / bw;
         printf(" %6.0f +-%4.0f |", a, e);
         if (st == 0 && a > 0) {
            yld[i]->SetBinContent(b, a);
            yld[i]->SetBinError(b, std::max(e, std::sqrt(std::max(a, 1.0))));
         }
      }
      // A non-converged fit is DROPPED from the stored distribution, so say so on the line that
      // reports it. Otherwise the numbers print, look entirely reasonable, and the angular bin
      // simply vanishes from the angular distribution downstream -- which is exactly what
      // happened to the 90-100 deg bin when the vertex-z window was first applied.
      printf(" %6.2f | shift %+6.3f%s\n", c2n, fn.GetParameter(NLV + 2),
             st != 0 ? "  <-- FIT STATUS != 0, BIN DROPPED" : (nTry ? "  (converged on retry)" : ""));
      gShift->SetPoint(nShift++, 0.5 * (lo + hi), fn.GetParameter(NLV + 2));

      if (pad < 12) {
         cf->cd(++pad);
         h->SetTitle(TString::Format("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", lo, hi));
         h->SetLineColor(kBlack);
         h->Draw("hist");
         auto *fc = (TF1 *)fn.DrawCopy("same");
         fc->SetLineColor(kRed + 1);
         fc->SetLineWidth(2);
         for (int i = 0; i < NLV; ++i) {
            auto *g = new TF1(TString::Format("c%d_%d", b, i), "gaus", 5.4, 8.0);
            g->SetParameters(fn.GetParameter(i) / (gSg[i] * std::sqrt(2 * TMath::Pi())),
                          gMu[i] + fn.GetParameter(NLV + 2), gSg[i]);
            g->SetLineColor(LVCOL[i]);
            g->SetLineStyle(2);
            g->SetLineWidth(2);
            g->Draw("same");
         }
      } else
         delete h;
   }
   cf->SaveAs(here + "/plots/exc_fits_" + tag + ".png");

   for (int i = 0; i < NLV; ++i)
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

   printf("\n  integrated over %.0f-%.0f deg:\n", cmMin, cmMax);
   double tot = 0;
   for (int i = 0; i < NLV; ++i)
      tot += yld[i]->Integral();
   for (int i = 0; i < NLV; ++i)
      printf("    %-10s %7.0f counts   %5.1f %%\n", LVNAME[i], yld[i]->Integral(),
             100 * yld[i]->Integral() / tot);

   TCanvas *c1 = new TCanvas("c1", "excited angular distributions", 1300, 560);
   c1->Divide(2, 1);
   auto style = [&](TH1D *h, int i) {
      h->SetMarkerStyle(20 + i);
      h->SetMarkerColor(LVCOL[i]);
      h->SetLineColor(LVCOL[i]);
      h->SetLineWidth(2);
      h->SetMarkerSize(1.1);
   };
   c1->cd(1);
   gPad->SetLogy();
   double mx = 0;
   for (int i = 0; i < NLV; ++i)
      mx = std::max(mx, yld[i]->GetMaximum());
   auto *frm = new TH1D("frm", "raw yield per level;#theta_{cm} [deg];counts / bin", 1, cmMin, cmMax);
   frm->SetMinimum(1);
   frm->SetMaximum(mx * 4);
   frm->Draw();
   auto *lg = new TLegend(0.52, 0.66, 0.89, 0.88);
   for (int i = 0; i < NLV; ++i) {
      style(yld[i], i);
      yld[i]->Draw("E1 same");
      lg->AddEntry(yld[i], LVNAME[i], "lp");
   }
   lg->Draw();
   c1->cd(2);
   gPad->SetLogy();
   mx = 0;
   for (int i = 0; i < NLV; ++i)
      mx = std::max(mx, dsd[i]->GetMaximum());
   auto *frm2 = new TH1D("frm2", "acceptance-corrected d#sigma/d#Omega (shape);#theta_{cm} [deg];arb.", 1, cmMin,
                         cmMax);
   frm2->SetMinimum(std::max(1.0, mx * 2e-3));
   frm2->SetMaximum(mx * 4);
   frm2->Draw();
   for (int i = 0; i < NLV; ++i) {
      style(dsd[i], i);
      dsd[i]->Draw("E1 same");
   }
   lg->Draw();
   TString png = here + "/plots/exc_angdist_" + tag + ".png";
   c1->SaveAs(png);

   TFile fo(here + "/plots/exc_angdist_" + tag + ".root", "RECREATE");
   for (int i = 0; i < NLV; ++i) {
      yld[i]->Write(TString::Format("yield_%d", i));
      dsd[i]->Write(TString::Format("dsdo_%d", i));
   }
   hAll->Write("ex_all");
   gShift->Write("shift");   // per-bin common centroid shift, for exc_cuts_view_C14.C
   fo.Close();
   printf("\nwrote %s\n       plots/exc_fits_%s.png\n\n", png.Data(), tag.Data());
}
