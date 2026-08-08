/// @file elastic_yield2_C14.C
/// @brief Elastic angular distribution by a TWO-PASS peak fit that follows the E_x locus.
///
/// Supersedes the single-pass elastic_yield_fit_C14.C, whose per-bin fit could collapse onto the
/// sigma lower limit in low-statistics bins and return a near-zero area (55-60, 115-135 deg all
/// came back with sigma pinned at 0.080). The failure mode matters: a collapsed fit does not look
/// like a failure, it looks like a physical dip.
///
/// Pass 1 -- locate. In every bin with enough entries, smooth the E_x spectrum and take the mode
/// within a band around the running seed, then fit a loose gaussian. Fit a smooth pol3 through
/// the resulting centroids and a pol2 through the widths, using only bins that passed.
/// Pass 2 -- integrate. Refit every bin with the centroid constrained to mu_smooth +- 0.3 and the
/// width to sigma_smooth x [0.6, 1.8], over mu_smooth +- 3 sigma_smooth. Yield = gaussian area.
///
/// The smooth locus is the point: the E_x drift is a property of the reconstruction, not of the
/// statistics in one angular bin, so it should be modelled globally and used to guide the sparse
/// bins rather than re-derived independently in each.
///
///   root -b -q 'elastic_yield2_C14.C("plots/proton_kin_300_ukf_nc.root","/mnt/f/a1954_C14_acc_nochi2/","ukf")'
///   root -b -q 'elastic_yield2_C14.C("plots/proton_kin_300gfx_nc.root","/mnt/f/a1954_C14_acc_gf_nochi2/","gf")'

void elastic_yield2_C14(TString dataCache = "plots/proton_kin_300_ukf_nc.root",
                        TString accDir = "/mnt/f/a1954_C14_acc_nochi2/", TString tag = "ukf",
                        Double_t cmMin = 20.0, Double_t cmMax = 145.0, Int_t minN = 60,
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
   if (!t || !acc)
      return;

   const int nb = acc->GetNbinsX();
   const double axLo = acc->GetXaxis()->GetXmin(), axHi = acc->GetXaxis()->GetXmax();
   auto mk = [&](const char *n) {
      auto *h = new TH1D(n, "", nb, axLo, axHi);
      h->Sumw2();
      h->SetDirectory(nullptr);
      return h;
   };
   TH1D *yield = mk("yield"), *dsdo = mk("dsdo"), *win = mk("win"), *winD = mk("winD");

   // one Ex spectrum per angular bin, kept for both passes
   std::vector<TH1D *> hx(nb + 1, nullptr);
   std::vector<double> Nent(nb + 1, 0.0);
   for (int b = 1; b <= nb; ++b) {
      double ctr = acc->GetBinCenter(b), lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      if (ctr < cmMin || ctr > cmMax)
         continue;
      auto *h = new TH1D(TString::Format("hx_%d", b), "", 160, -5.0, 3.0);
      t->Draw(TString::Format("ex>>hx_%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, lo + wid), "goff");
      h->SetDirectory(nullptr);
      hx[b] = h;
      Nent[b] = h->Integral();
   }

   // ---------------- pass 1: locate ----------------
   auto *gMu1 = new TGraph(), *gSg1 = new TGraph();
   int n1 = 0;
   double seed = 0.0;
   printf("\n===== %s : pass 1, locating the elastic locus =====\n", tag.Data());
   printf("  theta_cm   entries |    mode     mu      sigma\n");
   for (int b = 1; b <= nb; ++b) {
      if (!hx[b] || Nent[b] < 200)
         continue; // pass 1 only uses well-populated bins
      auto *h = (TH1D *)hx[b]->Clone(TString::Format("sm_%d", b));
      h->SetDirectory(nullptr);
      h->Smooth(2);
      int bm = 0;
      double vm = -1;
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         double c = h->GetBinCenter(i);
         if (std::fabs(c - seed) > 1.4)
            continue;
         if (h->GetBinContent(i) > vm) {
            vm = h->GetBinContent(i);
            bm = i;
         }
      }
      if (bm == 0) {
         delete h;
         continue;
      }
      double mode = h->GetBinCenter(bm);
      TF1 f1("f1", "gaus(0)+pol1(3)", mode - 0.9, mode + 0.9);
      f1.SetParameters(vm, mode, 0.3, 0.05 * vm, 0.0);
      f1.SetParLimits(1, mode - 0.5, mode + 0.5);
      f1.SetParLimits(2, 0.08, 1.3);
      if (hx[b]->Fit(&f1, "QNR") == 0 && f1.GetParameter(2) > 0.09 && f1.GetParameter(2) < 1.25) {
         gMu1->SetPoint(n1, acc->GetBinCenter(b), f1.GetParameter(1));
         gSg1->SetPoint(n1, acc->GetBinCenter(b), std::fabs(f1.GetParameter(2)));
         ++n1;
         seed = f1.GetParameter(1);
         printf("  %3.0f-%3.0f  %8.0f | %+6.2f  %+6.3f    %5.3f\n", acc->GetBinLowEdge(b),
                acc->GetBinLowEdge(b) + acc->GetBinWidth(b), Nent[b], mode, f1.GetParameter(1),
                std::fabs(f1.GetParameter(2)));
      }
      delete h;
   }
   if (n1 < 6) {
      printf("\033[1;31mpass 1 found only %d usable bins -- cannot build a locus\033[0m\n", n1);
      return;
   }
   TF1 fMu("fMu", "pol3", cmMin, cmMax);
   TF1 fSg("fSg", "pol2", cmMin, cmMax);
   gMu1->Fit(&fMu, "QN");
   gSg1->Fit(&fSg, "QN");
   printf("\nlocus   mu(theta)    = %+.4f %+.5f t %+.7f t^2 %+.9f t^3\n", fMu.GetParameter(0), fMu.GetParameter(1),
          fMu.GetParameter(2), fMu.GetParameter(3));
   printf("width   sigma(theta) = %+.4f %+.5f t %+.7f t^2\n", fSg.GetParameter(0), fSg.GetParameter(1),
          fSg.GetParameter(2));

   // ---------------- pass 2: integrate ----------------
   auto *gMu = new TGraphErrors(), *gSg = new TGraphErrors();
   int np = 0;
   printf("\n===== %s : pass 2, yield =====\n", tag.Data());
   printf("  theta_cm   entries | mu_smooth sg_smooth |    mu      sigma  |   peak yield    window | chi2/ndf\n");
   for (int b = 1; b <= nb; ++b) {
      if (!hx[b])
         continue;
      double ctr = acc->GetBinCenter(b), lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b);
      double ms = fMu.Eval(ctr), ss = std::max(0.12, std::min(1.1, fSg.Eval(ctr)));
      // the fixed |Ex|<0.6 window, for comparison, always from the same spectrum
      double wy = hx[b]->Integral(hx[b]->FindBin(-0.6), hx[b]->FindBin(0.6));
      win->SetBinContent(b, wy);
      win->SetBinError(b, std::sqrt(std::max(wy, 1.0)));
      if (Nent[b] < minN) {
         printf("  %3.0f-%3.0f  %8.0f | only %d entries -- DROPPED\n", lo, lo + wid, Nent[b], (int)Nent[b]);
         continue;
      }
      TF1 f2("f2", "gaus(0)+pol1(3)", ms - 3 * ss, ms + 3 * ss);
      double a0 = hx[b]->GetBinContent(hx[b]->FindBin(ms));
      f2.SetParameters(std::max(a0, 1.0), ms, ss, 0.05 * a0, 0.0);
      f2.SetParLimits(0, 0, 20 * std::max(a0, 1.0));
      f2.SetParLimits(1, ms - 0.3, ms + 0.3);
      f2.SetParLimits(2, 0.6 * ss, 1.8 * ss);
      int st = hx[b]->Fit(&f2, "QNR");
      double A = f2.GetParameter(0), mu = f2.GetParameter(1), sg = std::fabs(f2.GetParameter(2));
      double bw = hx[b]->GetBinWidth(1);
      double area = A * sg * std::sqrt(2 * TMath::Pi()) / bw;
      double rel = (A > 0) ? std::sqrt(std::pow(f2.GetParError(0) / A, 2) + std::pow(f2.GetParError(2) / sg, 2)) : 1;
      double c2n = f2.GetNDF() > 0 ? f2.GetChisquare() / f2.GetNDF() : -1;
      if (st != 0 || area <= 0 || !std::isfinite(area)) {
         printf("  %3.0f-%3.0f  %8.0f | fit failed (status %d) -- DROPPED\n", lo, lo + wid, Nent[b], st);
         continue;
      }
      yield->SetBinContent(b, area);
      yield->SetBinError(b, area * std::max(rel, 1.0 / std::sqrt(std::max(area, 1.0))));
      gMu->SetPoint(np, ctr, mu);
      gMu->SetPointError(np, 0, f2.GetParError(1));
      gSg->SetPoint(np, ctr, sg);
      gSg->SetPointError(np, 0, f2.GetParError(2));
      ++np;
      printf("  %3.0f-%3.0f  %8.0f | %+8.3f %9.3f | %+6.3f   %5.3f  | %10.0f %9.0f | %6.2f\n", lo, lo + wid, Nent[b],
             ms, ss, mu, sg, area, wy, c2n);
   }

   // acceptance + solid angle, for both the peak fit and the window
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
   printf("\nnormalisation to FRESCO over 20-50 deg: peak fit x %.5g, window x %.5g\n", kP, kW);
   printf("\n  theta_cm     FRESCO   peakfit  ratio |   window  ratio\n");
   double chiP = 0, chiW = 0;
   int nP = 0;
   for (int b = 1; b <= nb; ++b) {
      double c = dsdo->GetBinCenter(b);
      if (c < cmMin || c > cmMax)
         continue;
      double f = fr->Eval(c);
      double d = dsdo->GetBinContent(b) * kP, w = winD->GetBinContent(b) * kW;
      if (d <= 0 && w <= 0)
         continue;
      printf("  %3.0f-%3.0f %10.4g %9.4g %6.2f | %8.4g %6.2f\n", dsdo->GetBinLowEdge(b),
             dsdo->GetBinLowEdge(b) + dsdo->GetBinWidth(b), f, d, f > 0 ? d / f : 0, w, f > 0 ? w / f : 0);
      if (f > 0 && d > 0) {
         chiP += std::pow(std::log(d / f), 2);
         chiW += std::pow(std::log(std::max(w, 1e-9) / f), 2);
         ++nP;
      }
   }
   if (nP)
      printf("\nrms of ln(data/FRESCO) over %d bins:  peak fit %.3f   fixed window %.3f\n", nP,
             std::sqrt(chiP / nP), std::sqrt(chiW / nP));

   // ---------------- figure ----------------
   TCanvas *c1 = new TCanvas("c1", "peak-fit elastic v2", 1500, 950);
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
   gMu->SetTitle(TString::Format("%s: elastic E_{x} locus and width;#theta_{cm} [deg];MeV", tag.Data()));
   gMu->GetYaxis()->SetRangeUser(-2.5, 1.2);
   gMu->Draw("AP");
   fMu.SetLineColor(kAzure + 2);
   fMu.SetLineStyle(2);
   fMu.Draw("same");
   gSg->SetMarkerStyle(24);
   gSg->SetMarkerColor(kOrange + 7);
   gSg->SetLineColor(kOrange + 7);
   gSg->Draw("P same");
   fSg.SetLineColor(kOrange + 7);
   fSg.SetLineStyle(2);
   fSg.Draw("same");
   auto *z = new TLine(cmMin, 0, cmMax, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   auto *lg1 = new TLegend(0.55, 0.73, 0.89, 0.89);
   lg1->AddEntry(gMu, "centroid (pass 2)", "p");
   lg1->AddEntry(gSg, "#sigma (pass 2)", "p");
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
   auto *lg2 = new TLegend(0.52, 0.72, 0.89, 0.88);
   lg2->AddEntry(yield, "peak area (locus-tracking)", "lp");
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
   auto *lg3 = new TLegend(0.42, 0.70, 0.89, 0.88);
   lg3->AddEntry(fr, "FRESCO KD03, E_{lab}(p)=11.58 MeV", "l");
   lg3->AddEntry(dN, "peak area, acc-corrected", "lp");
   lg3->AddEntry(wN, "fixed window, acc-corrected", "lp");
   lg3->SetTextSize(0.030);
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

   TString png = here + "/plots/elastic_yield2_" + tag + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/elastic_yield2_" + tag + ".root", "RECREATE");
   yield->Write("yield_peak");
   win->Write("yield_window");
   dsdo->Write("dsigma_dOmega");
   winD->Write("dsigma_dOmega_window");
   gMu->Write("locus");
   gSg->Write("width");
   fo.Close();
   printf("\nwrote %s and .root\n\n", png.Data());
}
