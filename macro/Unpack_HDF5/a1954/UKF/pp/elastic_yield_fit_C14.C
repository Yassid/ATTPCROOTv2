/// @file elastic_yield_fit_C14.C
/// @brief Elastic angular distribution from a PEAK FIT per theta_cm bin, not a fixed Ex window.
///
/// ex_slices_C14.C showed the elastic locus is not at Ex = 0 beyond ~50 deg: it walks to about
/// -1.4 MeV by theta_cm ~ 100 and broadens from sigma 0.15 to 0.6-0.9. A window fixed at zero
/// therefore slides off the peak, which is why both fitters undershoot the FRESCO secondary
/// maximum by ~10x while their TOTAL track counts per angle bin stay equal. Counting the peak
/// where it actually is removes that artefact and is the honest way to get the yield.
///
/// Per 5 deg bin: gaussian + linear background, seeded from the previous bin's centroid so the
/// fit follows the ridge instead of jumping onto a background bump. Yield = gaussian area / bin
/// width (counts), so it is directly comparable to the window method. Bins whose fit does not
/// converge, or whose peak lands outside the search band, are reported and dropped rather than
/// silently filled.
///
/// Then: divide by the per-level acceptance (same file the window method used), divide out
/// sin(theta) for the dsigma/dOmega shape, and normalise to FRESCO over a forward window.
///
///   root -b -q 'elastic_yield_fit_C14.C("plots/proton_kin_300_ukf_nc.root","/mnt/f/a1954_C14_acc_nochi2/","ukf")'
///   root -b -q 'elastic_yield_fit_C14.C("plots/proton_kin_300gfx_nc.root","/mnt/f/a1954_C14_acc_gf_nochi2/","gf")'

void elastic_yield_fit_C14(TString dataCache = "plots/proton_kin_300_ukf_nc.root",
                           TString accDir = "/mnt/f/a1954_C14_acc_nochi2/", TString tag = "ukf",
                           Double_t cmMin = 20.0, Double_t cmMax = 150.0, Double_t exHalf = 1.6,
                           TString frFile = "/home/yassid/a1954_C14_fresco/outputs/p14C_el_161_dsdo.dat")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fd = TFile::Open(here + "/" + dataCache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   if (!fd || fd->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mcannot open data cache or acceptance\033[0m\n");
      return;
   }
   TTree *t = (TTree *)fd->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   if (!t || !acc) {
      printf("\033[1;31mmissing pk tree or hAcc_gs_sum\033[0m\n");
      return;
   }

   const int nb = acc->GetNbinsX();
   auto mk = [&](const char *n) {
      auto *h = new TH1D(n, "", nb, acc->GetXaxis()->GetXmin(), acc->GetXaxis()->GetXmax());
      h->Sumw2();
      h->SetDirectory(nullptr);
      return h;
   };
   TH1D *yield = mk("yield"), *cor = mk("cor"), *dsdo = mk("dsdo"), *win = mk("win");
   auto *gMu = new TGraphErrors(), *gSg = new TGraphErrors();

   printf("\n===== %s : elastic yield by peak fit =====\n", tag.Data());
   printf("  theta_cm   entries |     mu      sigma   |  peak yield   window yield  peak/window | chi2/ndf\n");
   double seed = 0.0;
   int np = 0;
   for (int b = 1; b <= nb; ++b) {
      double ctr = acc->GetBinCenter(b), lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      if (ctr < cmMin || ctr > cmMax)
         continue;
      auto *h = new TH1D(TString::Format("hx%d", b), "", 140, seed - 3.5, seed + 3.5);
      // TTree::Draw looks ">>name" up in the current directory -- detaching before the Draw
      // makes it fill a different, throwaway histogram and leaves this one empty.
      t->Draw(TString::Format("ex>>hx%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, lo + wid), "goff");
      h->SetDirectory(nullptr);
      double N = h->Integral();
      if (N < 40) {
         printf("  %3.0f-%3.0f  %8.0f | too few entries -- DROPPED\n", lo, lo + wid, N);
         delete h;
         continue;
      }
      // seed the centroid from the tallest bin within exHalf of the running seed
      int bm = 0;
      double vm = -1;
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         double c = h->GetBinCenter(i);
         if (std::fabs(c - seed) > exHalf)
            continue;
         if (h->GetBinContent(i) > vm) {
            vm = h->GetBinContent(i);
            bm = i;
         }
      }
      if (bm == 0) {
         printf("  %3.0f-%3.0f  %8.0f | no peak within +-%.1f of %.2f -- DROPPED\n", lo, lo + wid, N, exHalf, seed);
         delete h;
         continue;
      }
      double c0 = h->GetBinCenter(bm);
      // width guess grows with angle; let the fit move it
      double s0 = ctr < 55 ? 0.20 : 0.55;
      TF1 fn(TString::Format("fn%d", b), "gaus(0)+pol1(3)", c0 - 3 * s0, c0 + 3 * s0);
      fn.SetParameters(vm, c0, s0, 0.1 * vm, 0.0);
      fn.SetParLimits(0, 0, 10 * vm);
      fn.SetParLimits(1, c0 - 0.6, c0 + 0.6);
      fn.SetParLimits(2, 0.08, 1.2);
      int st = h->Fit(&fn, "QNR");
      double mu = fn.GetParameter(1), sg = std::fabs(fn.GetParameter(2));
      double A = fn.GetParameter(0);
      double bw = h->GetBinWidth(1);
      double area = A * sg * std::sqrt(2 * TMath::Pi()) / bw;
      double dA = fn.GetParError(0), dS = fn.GetParError(2);
      double relE = (A > 0 && sg > 0) ? std::sqrt((dA / A) * (dA / A) + (dS / sg) * (dS / sg)) : 1.0;
      double c2n = fn.GetNDF() > 0 ? fn.GetChisquare() / fn.GetNDF() : -1;

      // the fixed-window yield, for comparison
      double wy = h->Integral(h->FindBin(-0.6), h->FindBin(0.6));

      if (st != 0 || area <= 0 || !std::isfinite(area)) {
         printf("  %3.0f-%3.0f  %8.0f | fit failed (status %d) -- DROPPED\n", lo, lo + wid, N, st);
         delete h;
         continue;
      }
      yield->SetBinContent(b, area);
      yield->SetBinError(b, area * std::max(relE, 1.0 / std::sqrt(std::max(area, 1.0))));
      win->SetBinContent(b, wy);
      win->SetBinError(b, std::sqrt(std::max(wy, 1.0)));
      gMu->SetPoint(np, ctr, mu);
      gMu->SetPointError(np, 0, fn.GetParError(1));
      gSg->SetPoint(np, ctr, sg);
      gSg->SetPointError(np, 0, dS);
      ++np;
      printf("  %3.0f-%3.0f  %8.0f | %+6.3f    %5.3f   | %10.0f   %10.0f   %8.2f   | %6.2f\n", lo, lo + wid, N, mu, sg,
             area, wy, wy > 0 ? area / wy : 0, c2n);
      seed = mu; // follow the ridge
      delete h;
   }

   // acceptance correction + solid angle
   for (int b = 1; b <= nb; ++b) {
      double A = acc->GetBinContent(b), y = yield->GetBinContent(b);
      if (y <= 0 || A <= 0.05)
         continue;
      double C = y / A, relS = yield->GetBinError(b) / y;
      cor->SetBinContent(b, C);
      cor->SetBinError(b, C * relS);
      double s = std::sin(cor->GetBinCenter(b) * TMath::DegToRad());
      if (s > 1e-3) {
         dsdo->SetBinContent(b, C / s);
         dsdo->SetBinError(b, C * relS / s);
      }
   }

   // FRESCO overlay, normalised over 20-50 deg
   auto *fr = new TGraph();
   {
      std::ifstream in(frFile.Data());
      double th, xs;
      int n = 0;
      while (in >> th >> xs)
         fr->SetPoint(n++, th, xs);
   }
   double sn = 0, sd = 0;
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b), y = dsdo->GetBinContent(b), e = dsdo->GetBinError(b);
      if (c < 20 || c > 50 || y <= 0 || e <= 0)
         continue;
      double f = fr->Eval(c);
      sn += y * f / (e * e);
      sd += f * f / (e * e);
   }
   double s2f = sd > 0 ? sd / sn : 1.0; // multiply data by this to sit on FRESCO
   printf("\nnormalisation to FRESCO over 20-50 deg: x %.5g\n", s2f);
   printf("\n  theta_cm    FRESCO   data(peak fit)  ratio  |  data(window)  ratio\n");
   auto *winC = (TH1D *)win->Clone("winC");
   for (int b = 1; b <= nb; ++b) {
      double A = acc->GetBinContent(b);
      double c = win->GetBinCenter(b), s = std::sin(c * TMath::DegToRad());
      if (win->GetBinContent(b) > 0 && A > 0.05 && s > 1e-3)
         winC->SetBinContent(b, win->GetBinContent(b) / A / s);
      else
         winC->SetBinContent(b, 0);
   }
   double wn = 0, wd = 0;
   for (int b = 1; b <= nb; ++b) {
      double c = winC->GetBinCenter(b), y = winC->GetBinContent(b);
      if (c < 20 || c > 50 || y <= 0)
         continue;
      double f = fr->Eval(c);
      wn += y * f;
      wd += f * f;
   }
   double w2f = wn > 0 ? wd / wn : 1.0;
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b);
      if (c < cmMin || c > cmMax || dsdo->GetBinContent(b) <= 0)
         continue;
      double f = fr->Eval(c), d = dsdo->GetBinContent(b) * s2f, w = winC->GetBinContent(b) * w2f;
      printf("  %3.0f-%3.0f %9.4g %13.4g %7.2f  | %11.4g %7.2f\n", dsdo->GetBinLowEdge(b),
             dsdo->GetBinLowEdge(b) + dsdo->GetBinWidth(b), f, d, f > 0 ? d / f : 0, w, f > 0 ? w / f : 0);
   }

   // ---------- figure ----------
   TCanvas *c1 = new TCanvas("c1", "peak-fit elastic", 1500, 950);
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
   gMu->SetTitle(TString::Format("%s: elastic locus and width;#theta_{cm} [deg];E_{x} peak / #sigma [MeV]",
                                 tag.Data()));
   gMu->GetYaxis()->SetRangeUser(-2.2, 1.2);
   gMu->Draw("ALP");
   gSg->SetMarkerStyle(24);
   gSg->SetMarkerColor(kOrange + 7);
   gSg->SetLineColor(kOrange + 7);
   gSg->Draw("LP same");
   auto *z = new TLine(cmMin, 0, cmMax, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   auto *lg1 = new TLegend(0.58, 0.74, 0.89, 0.89);
   lg1->AddEntry(gMu, "peak centroid", "lp");
   lg1->AddEntry(gSg, "peak #sigma", "lp");
   lg1->Draw();

   c1->cd(2);
   gPad->SetLogy();
   style(yield, kAzure + 2, 20);
   style(win, kRed + 1, 21);
   yield->SetTitle(TString::Format("%s: yield, peak fit vs fixed window;#theta_{cm} [deg];counts", tag.Data()));
   yield->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   yield->SetMinimum(0.5);
   yield->Draw("E1");
   win->Draw("E1 same");
   auto *lg2 = new TLegend(0.55, 0.72, 0.89, 0.88);
   lg2->AddEntry(yield, "gaussian peak area", "lp");
   lg2->AddEntry(win, "fixed |E_{x}|<0.6 window", "lp");
   lg2->Draw();

   c1->cd(3);
   gPad->SetLogy();
   auto *frm = new TH1D("frm", TString::Format("%s: d#sigma/d#Omega vs FRESCO;#theta_{cm} [deg];mb/sr (shape)",
                                               tag.Data()),
                        1, cmMin - 5, cmMax + 5);
   frm->SetMinimum(0.3);
   frm->SetMaximum(3e3);
   frm->Draw();
   fr->SetLineWidth(3);
   fr->SetLineColor(kBlack);
   fr->Draw("L same");
   auto *dsN = (TH1D *)dsdo->Clone("dsN");
   dsN->Scale(s2f);
   style(dsN, kAzure + 2, 20);
   dsN->Draw("E1 same");
   auto *wN = (TH1D *)winC->Clone("wN");
   wN->Scale(w2f);
   style(wN, kRed + 1, 21);
   wN->Draw("E1 same");
   auto *lg3 = new TLegend(0.45, 0.70, 0.89, 0.88);
   lg3->AddEntry(fr, "FRESCO KD03", "l");
   lg3->AddEntry(dsN, "peak fit, acc-corrected", "lp");
   lg3->AddEntry(wN, "fixed window, acc-corrected", "lp");
   lg3->Draw();

   c1->cd(4);
   gPad->SetLogy();
   auto *rat = (TH1D *)dsdo->Clone("rat");
   rat->Reset();
   auto *ratW = (TH1D *)dsdo->Clone("ratW");
   ratW->Reset();
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b), f = fr->Eval(c);
      if (c < cmMin || c > cmMax || f <= 0)
         continue;
      if (dsdo->GetBinContent(b) > 0) {
         rat->SetBinContent(b, dsdo->GetBinContent(b) * s2f / f);
         rat->SetBinError(b, dsdo->GetBinError(b) * s2f / f);
      }
      if (winC->GetBinContent(b) > 0)
         ratW->SetBinContent(b, winC->GetBinContent(b) * w2f / f);
   }
   style(rat, kAzure + 2, 20);
   style(ratW, kRed + 1, 21);
   rat->SetTitle(TString::Format("%s: data / FRESCO;#theta_{cm} [deg];ratio", tag.Data()));
   rat->GetXaxis()->SetRangeUser(cmMin - 5, cmMax + 5);
   rat->SetMinimum(0.03);
   rat->SetMaximum(30);
   rat->Draw("E1");
   ratW->Draw("E1 same");
   auto *one = new TLine(cmMin - 5, 1, cmMax + 5, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();
   lg2->Draw();

   TString png = here + "/plots/elastic_peakfit_" + tag + ".png";
   c1->SaveAs(png);

   TFile fo(here + "/plots/elastic_peakfit_" + tag + ".root", "RECREATE");
   yield->Write("yield_peak");
   win->Write("yield_window");
   cor->Write("corrected");
   dsdo->Write("dsigma_dOmega");
   gMu->Write("locus");
   gSg->Write("width");
   fo.Close();
   printf("\nwrote %s and .root\n\n", png.Data());
}
