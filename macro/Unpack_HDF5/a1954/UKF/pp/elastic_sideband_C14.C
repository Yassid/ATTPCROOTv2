/// @file elastic_sideband_C14.C
/// @brief Elastic angular distribution by locus-tracking integration with sideband subtraction.
///
/// Replaces both earlier attempts. elastic_yield_fit_C14.C fitted a gaussian per bin and collapsed
/// onto the sigma limit in sparse bins; elastic_yield2_C14.C smoothed the locus with a pol3/pol2
/// which over-constrained the broad bins. The deeper problem with both is the gaussian itself: past
/// theta_cm ~ 70 deg the elastic peak is visibly non-gaussian (see plots/ex_slices_C14.png), so a
/// gaussian AREA is a shape assumption the data does not support.
///
/// What this does instead:
///
/// 1. LOCUS. Per theta_cm bin, smooth the E_x spectrum, take the mode within a band around the
///    running seed, then iterate a truncated centroid inside +-1.2 w of that mode. The width w is
///    measured from the HALF-MAXIMUM CROSSINGS (w = FWHM/2.355), which is insensitive to the tails
///    that break an RMS or a gaussian fit. Bins failing a quality test are left empty here.
/// 2. SMOOTHING. mu(theta) and w(theta) are smoothed by a 3-point moving average over the bins
///    that passed, then linearly interpolated (and held flat outside) onto every bin. Deliberately
///    NOT a global polynomial: the locus is monotone-ish but not low-order, and a pol3 forced the
///    80-110 deg region to the wrong place.
/// 3. INTEGRATION. Signal = raw counts in [mu - nSig*w, mu + nSig*w]. Background = a straight line
///    through the two sideband densities, integrated over the signal window:
///        B = W[(1-t) dL + t dR],  t = (Xmid - xL)/(xR - xL),  dL,R = counts/width
///        Var(B) = [W(1-t)/wL]^2 L + [W t/wR]^2 R,   Var(yield) = S + Var(B)
///    No peak shape is assumed anywhere. Yield, S/B and the windows are all reported per bin so a
///    bad bin is visible rather than silently plausible.
///
/// Then acceptance-correct, divide out sin(theta), and normalise to FRESCO over a forward window.
///
///   root -b -q 'elastic_sideband_C14.C("plots/proton_kin_300_ukf_nc.root","/mnt/f/a1954_C14_acc_nochi2/","ukf")'
///   root -b -q 'elastic_sideband_C14.C("plots/proton_kin_300gfx_nc.root","/mnt/f/a1954_C14_acc_gf_nochi2/","gf")'

namespace
{

/// half-maximum width of the structure containing bin bm, in x units; -1 if it cannot be measured
double fwhmAround(TH1D *h, int bm)
{
   double half = 0.5 * h->GetBinContent(bm);
   if (half <= 0)
      return -1;
   int i = bm;
   while (i > 1 && h->GetBinContent(i) > half)
      --i;
   double xL = h->GetBinCenter(i);
   int j = bm;
   while (j < h->GetNbinsX() && h->GetBinContent(j) > half)
      ++j;
   double xR = h->GetBinCenter(j);
   // refuse if either side ran into the histogram edge
   if (i <= 1 || j >= h->GetNbinsX())
      return -1;
   return xR - xL;
}

/// truncated centroid: iterate the mean inside +-halfWin of the current centre
double truncCentroid(TH1D *h, double c0, double halfWin, int nIter = 3)
{
   double c = c0;
   for (int it = 0; it < nIter; ++it) {
      double sw = 0, sxw = 0;
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double x = h->GetBinCenter(b), y = h->GetBinContent(b);
         if (y <= 0 || std::fabs(x - c) > halfWin)
            continue;
         sw += y;
         sxw += x * y;
      }
      if (sw <= 0)
         return c;
      c = sxw / sw;
   }
   return c;
}

} // namespace

/// `level` selects the acceptance ("gs" or "ex1") and `seed0` the E_x the locus search starts
/// from (0 for the elastic peak, ~6.1 for the 6.094 MeV state). Everything else is identical, so
/// the inelastic group goes through exactly the same extraction as the elastic one.
void elastic_sideband_C14(TString dataCache = "plots/proton_kin_300_ukf_nc.root",
                          TString accDir = "/mnt/f/a1954_C14_acc_nochi2/", TString tag = "ukf",
                          Double_t cmMin = 20.0, Double_t cmMax = 145.0, Double_t nSig = 2.0,
                          Double_t sbGap = 1.0, Double_t sbWidth = 2.0, Int_t minN = 60,
                          TString frFile = "",
                          TString level = "gs", Double_t seed0 = 0.0, Double_t exLo = -6.0,
                          Double_t exHi = 4.0, Double_t zMin = -1e9, Double_t zMax = 1e9)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   // FRESCO products live beside this analysis, in ../fresco; resolve them relative to
   // this macro so the comparison is reproducible outside one machine.
   if (frFile.IsNull())
      frFile = here + "/../fresco/outputs/p14C_el_161_dsdo.dat";

   TFile *fd = TFile::Open(here + "/" + dataCache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_" + level + ".root");
   if (!fd || fd->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mcannot open %s or the acceptance in %s\033[0m\n", dataCache.Data(), accDir.Data());
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_" + level + "_sum");
   if (!t || !acc) {
      printf("\033[1;31mmissing pk tree or hAcc_%s_sum\033[0m\n", level.Data());
      return;
   }

   const int nb = acc->GetNbinsX();
   auto mk = [&](const char *n) {
      auto *h = new TH1D(n, "", nb, acc->GetXaxis()->GetXmin(), acc->GetXaxis()->GetXmax());
      h->Sumw2();
      h->SetDirectory(nullptr);
      return h;
   };
   TH1D *yield = mk("yield"), *dsdo = mk("dsdo"), *win = mk("win"), *winD = mk("winD"), *sbFrac = mk("sbFrac");

   // ---- Ex spectra, one per angular bin ----
   std::vector<TH1D *> hx(nb + 1, nullptr);
   std::vector<double> Nent(nb + 1, 0.0);
   for (int b = 1; b <= nb; ++b) {
      double ctr = acc->GetBinCenter(b), lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      if (ctr < cmMin || ctr > cmMax)
         continue;
      auto *h = new TH1D(TString::Format("hsb_%d", b), "", 200, exLo, exHi);
      // The vertex-z slab, when one is asked for. It has to be applied HERE, on the same
      // spectra the locus and the sidebands are measured from, so that signal and background are
      // both taken from the selection the normalisation will use. Cutting afterwards would
      // subtract a background measured on a different event sample.
      TString cut = TString::Format("thcm>=%g&&thcm<%g", lo, lo + wid);
      if (zMin > -1e8 || zMax < 1e8)
         cut += TString::Format("&&vertexz>%g&&vertexz<%g", zMin, zMax);
      t->Draw(TString::Format("ex>>hsb_%d", b), cut, "goff");
      h->SetDirectory(nullptr);
      hx[b] = h;
      Nent[b] = h->Integral();
   }

   // ---- pass 1: locus and width where they can be measured ----
   std::vector<double> mu(nb + 1, NAN), wd(nb + 1, NAN);
   printf("\n===== %s : pass 1, locus by mode + truncated centroid, width by FWHM =====\n", tag.Data());
   printf("  theta_cm   entries |   mode     mu      w=FWHM/2.355 | status\n");
   double seed = seed0;
   for (int b = 1; b <= nb; ++b) {
      if (!hx[b])
         continue;
      double lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      if (Nent[b] < 200) {
         printf("  %3.0f-%3.0f  %8.0f | (too thin for pass 1, will interpolate)\n", lo, lo + wid, Nent[b]);
         continue;
      }
      auto *s = (TH1D *)hx[b]->Clone(TString::Format("sm_%d", b));
      s->SetDirectory(nullptr);
      s->Smooth(2);
      int bm = 0;
      double vm = -1;
      for (int i = 1; i <= s->GetNbinsX(); ++i) {
         double x = s->GetBinCenter(i);
         if (std::fabs(x - seed) > 1.3)
            continue;
         if (s->GetBinContent(i) > vm) {
            vm = s->GetBinContent(i);
            bm = i;
         }
      }
      if (bm == 0) {
         printf("  %3.0f-%3.0f  %8.0f | no structure within +-1.3 of %.2f -- skipped\n", lo, lo + wid, Nent[b], seed);
         delete s;
         continue;
      }
      double mode = s->GetBinCenter(bm);
      double fw = fwhmAround(s, bm);
      double w = (fw > 0) ? fw / 2.355 : NAN;
      if (!std::isfinite(w) || w < 0.06 || w > 1.6) {
         printf("  %3.0f-%3.0f  %8.0f | %+6.2f            FWHM unusable (%.2f) -- skipped\n", lo, lo + wid, Nent[b],
                mode, fw);
         delete s;
         continue;
      }
      double c = truncCentroid(hx[b], mode, 1.2 * w);
      mu[b] = c;
      wd[b] = w;
      seed = c;
      printf("  %3.0f-%3.0f  %8.0f | %+6.2f  %+6.3f        %6.3f | ok\n", lo, lo + wid, Nent[b], mode, c, w);
      delete s;
   }

   // ---- smoothing + interpolation of the locus (moving average, then linear interp) ----
   auto smoothInterp = [&](std::vector<double> &v) {
      std::vector<int> idx;
      for (int b = 1; b <= nb; ++b)
         if (std::isfinite(v[b]))
            idx.push_back(b);
      if (idx.size() < 3)
         return false;
      std::vector<double> sm(idx.size());
      for (size_t i = 0; i < idx.size(); ++i) {
         double s = 0;
         int n = 0;
         for (int d = -1; d <= 1; ++d) {
            int j = (int)i + d;
            if (j < 0 || j >= (int)idx.size())
               continue;
            s += v[idx[j]];
            ++n;
        }
         sm[i] = s / n;
      }
      for (size_t i = 0; i < idx.size(); ++i)
         v[idx[i]] = sm[i];
      // linear interpolation between measured bins, flat extrapolation outside
      for (int b = 1; b <= nb; ++b) {
         if (std::isfinite(v[b]))
            continue;
         int lo = -1, hi = -1;
         for (int j = b - 1; j >= 1; --j)
            if (std::isfinite(v[j])) { lo = j; break; }
         for (int j = b + 1; j <= nb; ++j)
            if (std::isfinite(v[j])) { hi = j; break; }
         if (lo > 0 && hi > 0) {
            double f = double(b - lo) / double(hi - lo);
            v[b] = v[lo] + f * (v[hi] - v[lo]);
         } else if (lo > 0)
            v[b] = v[lo];
         else if (hi > 0)
            v[b] = v[hi];
      }
      return true;
   };
   if (!smoothInterp(mu) || !smoothInterp(wd)) {
      printf("\033[1;31mtoo few usable bins to build a locus\033[0m\n");
      return;
   }

   // ---- pass 2: sideband-subtracted integration ----
   auto *gMu = new TGraph(), *gW = new TGraph();
   int ng = 0;
   printf("\n===== %s : pass 2, sideband-subtracted yield "
          "(signal +-%.1f w, sidebands %.1f-%.1f w) =====\n",
          tag.Data(), nSig, sbGap + nSig, sbGap + nSig + sbWidth);
   printf("  theta_cm  |   mu      w    | signal    bkg     yield   +-err |  S/B  | fixed-window\n");
   for (int b = 1; b <= nb; ++b) {
      if (!hx[b])
         continue;
      double lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      double m = mu[b], w = wd[b];
      double wy = hx[b]->Integral(hx[b]->FindBin(seed0 - 0.6), hx[b]->FindBin(seed0 + 0.6));
      win->SetBinContent(b, wy);
      win->SetBinError(b, std::sqrt(std::max(wy, 1.0)));
      if (Nent[b] < minN || !std::isfinite(m) || !std::isfinite(w)) {
         printf("  %3.0f-%3.0f   | only %.0f entries -- DROPPED\n", lo, lo + wid, Nent[b]);
         continue;
      }
      gMu->SetPoint(ng, acc->GetBinCenter(b), m);
      gW->SetPoint(ng, acc->GetBinCenter(b), w);
      ++ng;

      double x1 = m - nSig * w, x2 = m + nSig * w;
      double lA = m - (nSig + sbGap + sbWidth) * w, lB = m - (nSig + sbGap) * w;
      double rA = m + (nSig + sbGap) * w, rB = m + (nSig + sbGap + sbWidth) * w;
      auto integ = [&](double a, double c) {
         int ba = hx[b]->FindBin(a), bc = hx[b]->FindBin(c);
         return hx[b]->Integral(ba, bc);
      };
      double S = integ(x1, x2), L = integ(lA, lB), R = integ(rA, rB);
      double wL = lB - lA, wR = rB - rA, W = x2 - x1;
      double xL = 0.5 * (lA + lB), xR = 0.5 * (rA + rB), Xm = 0.5 * (x1 + x2);
      double dL = wL > 0 ? L / wL : 0, dR = wR > 0 ? R / wR : 0;
      double tt = (xR != xL) ? (Xm - xL) / (xR - xL) : 0.5;
      double B = W * ((1 - tt) * dL + tt * dR);
      double cL = W * (1 - tt) / std::max(wL, 1e-9), cR = W * tt / std::max(wR, 1e-9);
      double varB = cL * cL * L + cR * cR * R;
      double Y = S - B;
      double eY = std::sqrt(std::max(S, 0.0) + varB);
      if (Y <= 0 || eY <= 0 || Y < 2 * eY) {
         printf("  %3.0f-%3.0f   | %+6.3f %5.3f | %7.0f %7.1f %8.1f %7.1f | %5.2f | not significant -- DROPPED\n", lo,
                lo + wid, m, w, S, B, Y, eY, B > 0 ? Y / B : 0);
         continue;
      }
      yield->SetBinContent(b, Y);
      yield->SetBinError(b, eY);
      sbFrac->SetBinContent(b, B > 0 ? Y / B : 0);
      printf("  %3.0f-%3.0f   | %+6.3f %5.3f | %7.0f %7.1f %8.1f %7.1f | %5.2f | %8.0f\n", lo, lo + wid, m, w, S, B, Y,
             eY, B > 0 ? Y / B : 0, wy);
   }

   // ---- acceptance + solid angle ----
   for (int b = 1; b <= nb; ++b) {
      double A = acc->GetBinContent(b), s = std::sin(acc->GetBinCenter(b) * TMath::DegToRad());
      if (A <= 0.05 || s <= 1e-3)
         continue;
      if (yield->GetBinContent(b) > 0) {
         dsdo->SetBinContent(b, yield->GetBinContent(b) / A / s);
         dsdo->SetBinError(b, yield->GetBinError(b) / A / s);
      }
      if (win->GetBinContent(b) > 0) {
         winD->SetBinContent(b, win->GetBinContent(b) / A / s);
         winD->SetBinError(b, win->GetBinError(b) / A / s);
      }
   }

   // ---- FRESCO ----
   auto *fr = new TGraph();
   {
      std::ifstream in(frFile.Data());
      double th, xs;
      int n = 0;
      while (in >> th >> xs)
         fr->SetPoint(n++, th, xs);
   }
   auto norm = [&](TH1D *h) {
      double sn = 0, sd = 0;
      for (int b = 1; b <= nb; ++b) {
         double c = h->GetBinCenter(b), y = h->GetBinContent(b), e = h->GetBinError(b);
         if (c < 20 || c > 50 || y <= 0 || e <= 0)
            continue;
         double f = fr->Eval(c);
         sn += y * f / (e * e);
         sd += f * f / (e * e);
      }
      return sn > 0 ? sd / sn : 1.0;
   };
   double kP = norm(dsdo), kW = norm(winD);
   printf("\nnormalisation to FRESCO over 20-50 deg: sideband x %.5g, fixed window x %.5g\n", kP, kW);
   printf("\n  theta_cm     FRESCO   sideband  ratio |   window  ratio\n");
   double cP = 0, cW = 0;
   int nP = 0, nW = 0;
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b);
      if (c < cmMin || c > cmMax)
         continue;
      double f = fr->Eval(c), d = dsdo->GetBinContent(b) * kP, w = winD->GetBinContent(b) * kW;
      if (d <= 0 && w <= 0)
         continue;
      printf("  %3.0f-%3.0f %10.4g %9.4g %6.2f | %8.4g %6.2f\n", dsdo->GetBinLowEdge(b),
             dsdo->GetBinLowEdge(b) + dsdo->GetBinWidth(b), f, d, f > 0 ? d / f : 0, w, f > 0 ? w / f : 0);
      if (f > 0 && d > 0) {
         cP += std::pow(std::log(d / f), 2);
         ++nP;
      }
      if (f > 0 && w > 0) {
         cW += std::pow(std::log(w / f), 2);
         ++nW;
      }
   }
   if (nP && nW)
      printf("\nrms of ln(data/FRESCO):  sideband %.3f over %d bins   |   fixed window %.3f over %d bins\n",
             std::sqrt(cP / nP), nP, std::sqrt(cW / nW), nW);

   // ---- figure ----
   TCanvas *c1 = new TCanvas("c1", "sideband elastic", 1500, 950);
   c1->Divide(2, 2);
   auto style = [](TH1D *h, int col, int m) {
      h->SetMarkerStyle(m);
      h->SetMarkerColor(col);
      h->SetLineColor(col);
      h->SetLineWidth(2);
   };
   c1->cd(1);
   gMu->SetMarkerStyle(20);
   gMu->SetMarkerColor(kAzure + 2);
   gMu->SetLineColor(kAzure + 2);
   gMu->SetLineWidth(2);
   gMu->SetTitle(TString::Format("%s: tracked locus and width;#theta_{cm} [deg];MeV", tag.Data()));
   gMu->GetYaxis()->SetRangeUser(-2.5, 1.2);
   gMu->Draw("ALP");
   gW->SetMarkerStyle(24);
   gW->SetMarkerColor(kOrange + 7);
   gW->SetLineColor(kOrange + 7);
   gW->SetLineWidth(2);
   gW->Draw("LP same");
   auto *z = new TLine(cmMin, 0, cmMax, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   auto *lg1 = new TLegend(0.58, 0.74, 0.89, 0.89);
   lg1->AddEntry(gMu, "locus #mu", "lp");
   lg1->AddEntry(gW, "width w", "lp");
   lg1->Draw();

   c1->cd(2);
   gPad->SetLogy();
   style(yield, kAzure + 2, 20);
   style(win, kRed + 1, 21);
   yield->SetTitle(TString::Format("%s: raw yield;#theta_{cm} [deg];counts", tag.Data()));
   yield->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   yield->SetMinimum(0.5);
   yield->Draw("E1");
   win->Draw("E1 same");
   auto *lg2 = new TLegend(0.50, 0.72, 0.89, 0.88);
   lg2->AddEntry(yield, "sideband-subtracted", "lp");
   lg2->AddEntry(win, "fixed |E_{x}|<0.6", "lp");
   lg2->Draw();

   c1->cd(3);
   gPad->SetLogy();
   auto *frm = new TH1D("frm", TString::Format("%s: d#sigma/d#Omega vs FRESCO;#theta_{cm} [deg];mb/sr (shape)",
                                               tag.Data()),
                        1, cmMin - 5, cmMax + 5);
   frm->SetMinimum(0.2);
   frm->SetMaximum(3e3);
   frm->Draw();
   fr->SetLineWidth(3);
   fr->SetLineColor(kBlack);
   fr->Draw("L same");
   auto *dN = (TH1D *)dsdo->Clone("dN");
   dN->Scale(kP);
   style(dN, kAzure + 2, 20);
   dN->Draw("E1 same");
   auto *wN = (TH1D *)winD->Clone("wN");
   wN->Scale(kW);
   style(wN, kRed + 1, 21);
   wN->Draw("E1 same");
   auto *lg3 = new TLegend(0.40, 0.70, 0.89, 0.88);
   lg3->AddEntry(fr, "FRESCO KD03, E_{lab}(p)=11.58 MeV", "l");
   lg3->AddEntry(dN, "sideband-subtracted, acc-corrected", "lp");
   lg3->AddEntry(wN, "fixed window, acc-corrected", "lp");
   lg3->SetTextSize(0.029);
   lg3->Draw();

   c1->cd(4);
   gPad->SetLogy();
   auto *rP = (TH1D *)dsdo->Clone("rP");
   rP->Reset();
   auto *rW = (TH1D *)winD->Clone("rW");
   rW->Reset();
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b), f = fr->Eval(c);
      if (c < cmMin || c > cmMax || f <= 0)
         continue;
      if (dsdo->GetBinContent(b) > 0) {
         rP->SetBinContent(b, dsdo->GetBinContent(b) * kP / f);
         rP->SetBinError(b, dsdo->GetBinError(b) * kP / f);
      }
      if (winD->GetBinContent(b) > 0) {
         rW->SetBinContent(b, winD->GetBinContent(b) * kW / f);
         rW->SetBinError(b, winD->GetBinError(b) * kW / f);
      }
   }
   style(rP, kAzure + 2, 20);
   style(rW, kRed + 1, 21);
   rP->SetTitle(TString::Format("%s: data / FRESCO;#theta_{cm} [deg];ratio", tag.Data()));
   rP->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   rP->SetMinimum(0.03);
   rP->SetMaximum(30);
   rP->Draw("E1");
   rW->Draw("E1 same");
   auto *one = new TLine(cmMin - 5, 1, cmMax + 5, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();
   lg2->Draw();

   TString png = here + "/plots/elastic_sideband_" + level + "_" + tag + ".png";
   c1->SaveAs(png);

   // ---- window diagnostic: six representative bins with the bands drawn ----
   TCanvas *c2 = new TCanvas("c2", "windows", 1500, 900);
   c2->Divide(3, 2);
   int shown = 0;
   for (int b = 1; b <= nb && shown < 6; ++b) {
      if (!hx[b] || yield->GetBinContent(b) <= 0)
         continue;
      double ctr = acc->GetBinCenter(b);
      if (!(ctr > 22 && (shown == 0 || ctr > 35) && (shown < 2 || ctr > 60) && (shown < 3 || ctr > 80) &&
            (shown < 4 || ctr > 100) && (shown < 5 || ctr > 115)))
         continue;
      c2->cd(++shown);
      double m = mu[b], w = wd[b];
      hx[b]->SetTitle(TString::Format("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", acc->GetBinLowEdge(b),
                                      acc->GetBinLowEdge(b) + acc->GetBinWidth(b)));
      hx[b]->GetXaxis()->SetRangeUser(m - 5 * w, m + 5 * w);
      hx[b]->SetLineColor(kBlack);
      hx[b]->Draw("hist");
      double ymax = hx[b]->GetMaximum() * 1.1;
      auto band = [&](double a, double c, int col) {
         auto *bx = new TBox(a, 0, c, ymax);
         bx->SetFillColorAlpha(col, 0.22);
         bx->SetLineColor(col);
         bx->Draw("same");
      };
      band(m - nSig * w, m + nSig * w, kAzure + 2);
      band(m - (nSig + sbGap + sbWidth) * w, m - (nSig + sbGap) * w, kOrange + 7);
      band(m + (nSig + sbGap) * w, m + (nSig + sbGap + sbWidth) * w, kOrange + 7);
      hx[b]->Draw("hist same");
   }
   TString png2 = here + "/plots/elastic_sideband_windows_" + level + "_" + tag + ".png";
   c2->SaveAs(png2);

   TFile fo(here + "/plots/elastic_sideband_" + level + "_" + tag + ".root", "RECREATE");
   yield->Write("yield_sideband");
   win->Write("yield_window");
   dsdo->Write("dsigma_dOmega");
   winD->Write("dsigma_dOmega_window");
   gMu->Write("locus");
   gW->Write("width");
   fo.Close();
   printf("\nwrote %s\n       %s\n\n", png.Data(), png2.Data());
}
